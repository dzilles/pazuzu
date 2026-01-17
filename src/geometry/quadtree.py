import warp as wp
import numpy as np
from typing import Tuple, Optional
from src.kernels.grid_kernels import compute_block_coordinates, generate_morton_codes
from src.kernels.connectivity_kernels import init_hash_map, populate_hash_map, compute_neighbors, mark_mortar_neighbors
from src.kernels.amr_kernels import mark_blocks_gradient, balance_refine_flags, prolongate_batch, restrict_batch, zero_blocks, mark_blocks_on_interface

ROOT_BOUNDS = (-1.0, -1.0, 1.0, 1.0) # x_min, y_min, x_max, y_max

class Quadtree:
    """
    Block-AMR Quadtree using Morton Encoding (Z-ordering) for implicit connectivity.
    
    Manages the grid hierarchy and mapping to the memory pool.
    """
    # --- Static Morton Encoding Helpers ---
    @staticmethod
    def part1by1(n: int) -> int:
        """Inserts a 0 bit after each of the low 16 bits of n."""
        n = (n ^ (n << 8)) & 0x00ff00ff
        n = (n ^ (n << 4)) & 0x0f0f0f0f
        n = (n ^ (n << 2)) & 0x33333333
        n = (n ^ (n << 1)) & 0x55555555
        return n
    
    @staticmethod
    def compact1by1(n: int) -> int:
        """Inverse of part1by1. Extracts bits at even positions."""
        n = n & 0x55555555
        n = (n ^ (n >> 1)) & 0x33333333
        n = (n ^ (n >> 2)) & 0x0f0f0f0f
        n = (n ^ (n >> 4)) & 0x00ff00ff
        n = (n ^ (n >> 8)) & 0x0000ffff
        return n

    @staticmethod
    def morton_encode(x: int, y: int, level: int) -> int:
        """
        Interleaves bits of x and y and adds a sentinel bit at (1 << 2*level).
        Checks for overflow to prevent collisions.
        """
        limit = 1 << level
        if x >= limit or y >= limit or x < 0 or y < 0:
            raise ValueError(f"Coordinates ({x}, {y}) out of bounds for level {level} (max {limit-1})")
        
        interleaved = Quadtree.part1by1(x) | (Quadtree.part1by1(y) << 1)
        return (1 << (2 * level)) | interleaved
    
    @staticmethod
    def morton_decode(code: int, level: int) -> Tuple[int, int]:
        """Decodes Morton code into (x, y) by removing the sentinel bit."""
        mask = (1 << (2 * level)) - 1
        interleaved = code & mask
        return Quadtree.compact1by1(interleaved), Quadtree.compact1by1(interleaved >> 1)

    def __init__(self, device: str = "cpu", max_blocks: int = 10000, root_bounds: Tuple[float, float, float, float] = (-1.0, -1.0, 1.0, 1.0), periodic_x: bool = False, periodic_y: bool = False, dtype=wp.float32, max_depth: int = 10):
        self.device = device
        self.max_blocks = max_blocks
        self.num_blocks = 0
        self.root_bounds = root_bounds
        self.periodic_x = periodic_x
        self.periodic_y = periodic_y
        self.dtype = dtype
        self.max_depth = max_depth
        
        # Morton codes for blocks resident in the pool.
        # Index i in this array corresponds to block i in SimulationState.
        self.block_morton_codes = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        self.block_levels = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        
        # Free List Management (Stack of free indices)
        # Initialize with all indices 0..max_blocks-1
        self.free_pool_indices: wp.array = wp.array(
            np.arange(max_blocks, dtype=np.int32), 
            dtype=wp.int32, 
            device=device
        )
        self.num_free = max_blocks
        
        # Hash Map for Morton -> Pool Index Lookups
        # Capacity should be larger than max_blocks to reduce collisions (e.g., 2x)
        self.map_capacity = max_blocks * 2
        self.map_keys = wp.zeros(self.map_capacity, dtype=wp.int32, device=device)
        self.map_values = wp.zeros(self.map_capacity, dtype=wp.int32, device=device)
        
        # Bounds as Warp vector
        vec4_type = wp.vec4d if dtype == wp.float64 else wp.vec4
        self.root_bounds_wp = vec4_type(root_bounds[0], root_bounds[1], root_bounds[2], root_bounds[3])
        
        # Mortar Interfaces (Phase 4)
        self.max_mortars = max_blocks * 4
        self.mortar_list = wp.zeros((self.max_mortars, 4), dtype=wp.int32, device=device)
        self.num_mortars = wp.zeros(1, dtype=wp.int32, device=device)

    def uniform_refine(self, level: int, state, basis):
        """
        Generates a uniform grid at the specified refinement level.
        
        Args:
            level (int): Refinement level (0 = root).
            state (SimulationState): The simulation state to populate.
            basis (Basis): The basis defining the nodes.
        """
        print(f"Generating uniform grid at level {level}...")
        
        grid_dim = 1 << level
        num_blocks = grid_dim * grid_dim
        
        if num_blocks > self.max_blocks:
             raise ValueError(f"Level {level} requires {num_blocks} blocks, but MAX_BLOCKS is {self.max_blocks}")

        self.num_blocks = num_blocks
        self.num_free = self.max_blocks - num_blocks # Remaining free blocks
        
        # 1. Generate Morton Codes on Device
        wp.launch(
            kernel=generate_morton_codes,
            dim=num_blocks,
            inputs=[self.block_morton_codes, grid_dim, level],
            device=self.device
        )
        
        # Set levels (all same)
        self.block_levels.fill_(level) 
        
        # 3. Update Active Block Indices in State
        # In a compact pool, the active indices are just 0, 1, ..., num_blocks-1
        active_indices = np.arange(num_blocks, dtype=np.int32)
        wp.copy(state.active_block_indices, wp.array(active_indices, dtype=wp.int32, device=self.device), count=num_blocks)
        
        # Update Free List (shift it, effectively popping the first num_blocks)
        # We just need to ensure the remaining free indices are valid.
        # Since we initialized 0..max, the remaining are num_blocks..max.
        remaining_free = np.arange(num_blocks, self.max_blocks, dtype=np.int32)
        wp.copy(self.free_pool_indices, wp.array(remaining_free, dtype=wp.int32, device=self.device), count=len(remaining_free))

        # 4. Compute Physical Coordinates
        wp.launch(
            kernel=compute_block_coordinates,
            dim=(num_blocks, basis.Np),
            inputs=[
                self.block_morton_codes,
                state.active_block_indices,
                num_blocks,
                self.block_levels,
                self.root_bounds_wp,
                basis.nodes_2d,
                state.x,
                state.y
            ],
            device=self.device
        )
        
        # 5. Build Connectivity
        self.build_connectivity(state)
        
        print(f"Created {self.num_blocks} blocks.")

    def adapt_mesh(self, state, basis, refine_threshold: float, coarsen_threshold: float, ibm_enabled: bool = False):
        """
        Adapts the mesh by refining and coarsening active blocks based on gradient thresholds.
        """
        if refine_threshold is None:
            return

        # --- Step A: Mark ---
        refine_flags = wp.zeros(self.max_blocks, dtype=wp.int32, device=self.device)
        
        wp.launch(
            kernel=mark_blocks_gradient,
            dim=self.num_blocks,
            inputs=[
                state.q,
                state.active_block_indices,
                self.num_blocks,
                self.block_levels,
                refine_flags,
                refine_threshold,
                coarsen_threshold,
                self.max_depth
            ],
            device=self.device
        )
        
        # --- Step A.2: Mark (Geometry Based - IBM) ---
        # This overrides gradient decisions to ensure the geometry is resolved.
        if ibm_enabled:
             wp.launch(
                 kernel=mark_blocks_on_interface,
                 dim=self.num_blocks,
                 inputs=[
                     state.phi,
                     state.active_block_indices,
                     self.num_blocks,
                     refine_flags
                 ],
                 device=self.device
             )
        
        # --- Step B: Balance (2:1 Constraint) ---
        # Note: Balancing only applies to refinement flags (1).
        changed = wp.zeros(1, dtype=wp.int32, device=self.device)
        h_changed = np.array([1], dtype=np.int32)
        wp_h_changed = wp.from_numpy(h_changed, dtype=wp.int32, device=self.device)
        
        while h_changed[0] > 0:
            changed.zero_()
            wp.launch(
                kernel=balance_refine_flags,
                dim=self.num_blocks,
                inputs=[
                    state.active_block_indices,
                    self.num_blocks,
                    self.block_morton_codes,
                    self.block_levels,
                    self.map_keys,
                    self.map_values,
                    self.map_capacity,
                    int(self.periodic_x),
                    int(self.periodic_y),
                    refine_flags,
                    changed
                ],
                device=self.device
            )
            wp.copy(wp_h_changed, changed)
            h_changed = wp_h_changed.numpy()

        # --- Step C: Filter and Execute Refinement ---
        h_flags = refine_flags.numpy()
        h_active = state.active_block_indices.numpy()[:self.num_blocks]
        
        blocks_to_refine = []
        for idx in h_active:
            if h_flags[idx] == 1:
                blocks_to_refine.append(idx)
                
        # Refine first
        if blocks_to_refine:
            self.refine_blocks(blocks_to_refine, state, basis)
            # Re-fetch active blocks after refinement for coarsening phase
            h_active = state.active_block_indices.numpy()[:self.num_blocks]

        # --- Step D: Execute Coarsening ---
        self.coarsen_marked_blocks(state, basis, h_flags)

    def refine_marked_blocks(self, state, basis, threshold: float):
        """
        Backward compatibility wrapper for refinement.
        """
        self.adapt_mesh(state, basis, refine_threshold=threshold, coarsen_threshold=-1.0)

    def refine_blocks(self, blocks_to_refine: list, state, basis):
        """
        Refines explicit blocks by splitting them into 4 children.
        Includes data prolongation and topology updates.
        """
        if not blocks_to_refine:
            return

        # Validate max depth
        h_levels = self.block_levels.numpy()
        valid_blocks = []
        for idx in blocks_to_refine:
            if h_levels[idx] < self.max_depth:
                valid_blocks.append(idx)
        
        blocks_to_refine = valid_blocks
        num_refine = len(blocks_to_refine)
        if num_refine == 0:
            return

        needed = num_refine * 4
        if self.num_free < needed:
            print(f"Warning: Not enough free blocks for refinement. Needed {needed}, have {self.num_free}")
            return 

        # --- Step C: Allocate & Topology Update ---
        h_morton_codes = self.block_morton_codes.numpy()
        h_free_indices = self.free_pool_indices.numpy()
        h_active = state.active_block_indices.numpy()[:self.num_blocks]
        free_ptr = 0
        
        new_active_list = []
        parents_to_remove = set(blocks_to_refine)
        
        # Build Op List for Prolongation: [parent_idx, child_idx, child_quadrant]
        prolongation_ops = [] 
        
        # Process Refinement
        for p_idx in blocks_to_refine:
            p_code = h_morton_codes[p_idx]
            p_level = h_levels[p_idx]
            
            # Get 4 children indices
            c_indices = h_free_indices[free_ptr : free_ptr + 4]
            free_ptr += 4
            
            for i in range(4):
                c_idx = c_indices[i]
                c_code = (p_code << 2) | i
                c_level = p_level + 1
                
                h_morton_codes[c_idx] = c_code
                h_levels[c_idx] = c_level
                
                # Add to Op List
                prolongation_ops.append([p_idx, c_idx, i])
        
        # --- Step D: Prolongate (Data Transfer) ---
        num_ops = len(prolongation_ops)
        ops_array: wp.array = wp.array(np.array(prolongation_ops, dtype=np.int32), dtype=wp.int32, device=self.device)
        
        wp.launch(
            kernel=prolongate_batch,
            dim=num_ops * basis.Np,
            inputs=[
                state.q,
                ops_array,
                num_ops,
                basis.P_left,
                basis.P_right
            ],
            device=self.device
        )
        
        # --- Step E: Update Active List & Recycle ---
        # Add non-refined blocks
        for idx in h_active:
            if idx not in parents_to_remove:
                new_active_list.append(idx)
        
        # Add ALL new children
        used_children = h_free_indices[0 : free_ptr]
        new_active_list.extend(used_children)
        
        # Recycle parents
        remaining_free = h_free_indices[free_ptr : self.num_free]
        freed_parents = list(parents_to_remove)
        
        new_free_count = len(remaining_free) + len(freed_parents)
        h_new_free = np.zeros(self.max_blocks, dtype=np.int32)
        h_new_free[:len(remaining_free)] = remaining_free
        h_new_free[len(remaining_free):new_free_count] = freed_parents
        
        # Upload
        wp.copy(self.free_pool_indices, wp.array(h_new_free, dtype=wp.int32, device="cpu"), count=new_free_count)
        self.num_free = new_free_count
        
        self.num_blocks = len(new_active_list)
        wp.copy(state.active_block_indices, wp.array(np.array(new_active_list, dtype=np.int32), dtype=wp.int32, device="cpu"), count=self.num_blocks)
        
        wp.copy(self.block_morton_codes, wp.array(h_morton_codes, dtype=wp.int32, device="cpu"))
        wp.copy(self.block_levels, wp.array(h_levels, dtype=wp.int32, device="cpu"))
        
        # --- Step F: Rebuild ---
        # Recompute coordinates
        wp.launch(
            kernel=compute_block_coordinates,
            dim=(self.num_blocks, basis.Np),
            inputs=[
                self.block_morton_codes,
                state.active_block_indices, 
                self.num_blocks,            
                self.block_levels,          
                self.root_bounds_wp,
                basis.nodes_2d,
                state.x,
                state.y
            ],
            device=self.device
        )
        self.build_connectivity(state)


    def coarsen_marked_blocks(self, state, basis, refine_flags: Optional[np.ndarray] = None):
        """
        Coarsens families of 4 active sibling blocks into their parent if all are marked for coarsening.
        
        Constraint: 2:1 Balance. A family cannot be coarsened if any of its members has a neighbor
        that is finer (i.e., requires the current fine level to exist). This is checked by verifying
        if any sibling is adjacent to a mortar interface where the neighbor is finer.
        """
        if refine_flags is None:
            # Default: mark all blocks for potential coarsening
            refine_flags = np.full(self.max_blocks, -1, dtype=np.int32)

        # --- Identify Coarsenable Families ---
        h_active = state.active_block_indices.numpy()[:self.num_blocks]
        h_codes = self.block_morton_codes.numpy()
        h_levels = self.block_levels.numpy()
        
        # Identify blocks that have FINER neighbors (using mortar list)
        has_finer = set()
        num_mortars_host = int(self.num_mortars.numpy()[0])
        if num_mortars_host > 0:
            h_mortars = self.mortar_list.numpy()[:num_mortars_host]
            for i in range(num_mortars_host):
                has_finer.add(h_mortars[i, 2]) # index 2 is coarse_idx
        
        candidates: dict[tuple[int, int], list[int]] = {} # parent_key -> [child_pool_idx, ...]
        
        # Group ALL active leaf nodes by their potential parent
        for idx in h_active:
            level = h_levels[idx]
            if level > 0:
                code = h_codes[idx]
                parent_code = code >> 2
                parent_level = level - 1
                key = (parent_code, parent_level)
                if key not in candidates:
                    candidates[key] = []
                candidates[key].append(idx)
        
        families_to_coarsen = []
        parents_created = [] # (parent_pool_idx, parent_code, parent_level, [child_indices])
        restriction_ops = [] # [child_idx, parent_idx, child_quadrant]
        children_set = set()

        for (p_code, p_level), members in candidates.items():
            if len(members) == 4:
                # Criteria for coarsening:
                # 1. ALL 4 siblings must be marked for coarsening (flag -1)
                # 2. No sibling can have a finer neighbor (prevents 2:1 violation with parent)
                
                can_coarsen = True
                for c_idx in members:
                    if refine_flags[c_idx] != -1:
                        can_coarsen = False
                        break
                    
                    # If any sibling has a finer neighbor, it cannot coarsen
                    if c_idx in has_finer:
                        can_coarsen = False
                        break
                
                if can_coarsen:
                    families_to_coarsen.append(((p_code, p_level), members))

        if not families_to_coarsen:
            return

        # --- Allocation & Topology ---
        h_free_indices = self.free_pool_indices.numpy()
        free_ptr = 0
        
        for (p_code, p_level), children in families_to_coarsen:
            # Pop 1 parent
            p_idx = h_free_indices[free_ptr]
            free_ptr += 1
            
            # Setup Parent
            h_codes[p_idx] = p_code
            h_levels[p_idx] = p_level
            
            parents_created.append((p_idx, children))
            
            # Determine quadrants for children
            for c_idx in children:
                c_code = h_codes[c_idx]
                quad = c_code & 3
                restriction_ops.append([c_idx, p_idx, quad])
                children_set.add(c_idx)
        
        # --- Data Transfer (Restriction) ---
        # 0. Zero Parents
        parent_indices_to_zero = np.array([p[0] for p in parents_created], dtype=np.int32)
        num_parents = len(parent_indices_to_zero)
        wp_parents_to_zero: wp.array = wp.array(parent_indices_to_zero, dtype=wp.int32, device=self.device)
        
        wp.launch(
            kernel=zero_blocks,
            dim=(num_parents, basis.Np),
            inputs=[state.q, wp_parents_to_zero, num_parents],
            device=self.device
        )

        # 1. Restrict
        num_ops = len(restriction_ops)
        ops_array: wp.array = wp.array(np.array(restriction_ops, dtype=np.int32), dtype=wp.int32, device=self.device)
        
        wp.launch(
            kernel=restrict_batch,
            dim=num_ops * basis.Np,
            inputs=[
                state.q,
                ops_array,
                num_ops,
                basis.R_left,
                basis.R_right
            ],
            device=self.device
        )
        
        # --- Update Host Arrays (Active/Free) ---
        new_active_list = []
        for idx in h_active:
            if idx not in children_set:
                new_active_list.append(idx)
        
        for p_idx, _ in parents_created:
            new_active_list.append(p_idx)
            
        self.num_blocks = len(new_active_list)
        wp.copy(state.active_block_indices, wp.array(np.array(new_active_list, dtype=np.int32), dtype=wp.int32, device="cpu"), count=self.num_blocks)
        
        remaining_free = h_free_indices[free_ptr : self.num_free]
        freed_children_list = list(children_set)
        freed_children = np.array(freed_children_list, dtype=np.int32)
        
        num_remaining = len(remaining_free)
        num_freed = len(freed_children)
        new_free_count = num_remaining + num_freed
        
        h_new_free = np.zeros(self.max_blocks, dtype=np.int32)
        h_new_free[:num_remaining] = remaining_free
        h_new_free[num_remaining : num_remaining + num_freed] = freed_children
        
        wp.copy(self.free_pool_indices, wp.array(h_new_free, dtype=wp.int32, device="cpu"), count=new_free_count)
        self.num_free = new_free_count
        
        wp.copy(self.block_morton_codes, wp.array(h_codes, dtype=wp.int32, device="cpu"))
        wp.copy(self.block_levels, wp.array(h_levels, dtype=wp.int32, device="cpu"))
        
        # --- Rebuild ---
        wp.launch(
            kernel=compute_block_coordinates,
            dim=(self.num_blocks, basis.Np),
            inputs=[
                self.block_morton_codes,
                state.active_block_indices,
                self.num_blocks,
                self.block_levels,
                self.root_bounds_wp,
                basis.nodes_2d,
                state.x,
                state.y
            ],
            device=self.device
        )
        self.build_connectivity(state)

    def build_connectivity(self, state, level: int = -1):
        """
        Rebuilds the hash map and computes neighbors for all active blocks.
        """
        # 1. Initialize Hash Map (Clear)
        wp.launch(
            kernel=init_hash_map,
            dim=self.map_capacity,
            inputs=[self.map_keys, self.map_values],
            device=self.device
        )
        
        # 2. Populate Hash Map
        wp.launch(
            kernel=populate_hash_map,
            dim=self.num_blocks,
            inputs=[
                self.block_morton_codes,
                self.block_levels,          # Added
                state.active_block_indices,
                self.num_blocks,
                self.map_keys,
                self.map_values,
                self.map_capacity
            ],
            device=self.device
        )
        
        # 3. Compute Neighbors
        self.num_mortars.zero_()
        
        wp.launch(
            kernel=compute_neighbors,
            dim=self.num_blocks,
            inputs=[
                self.block_morton_codes,
                self.block_levels,          # Added
                state.active_block_indices,
                self.num_blocks,
                self.map_keys,
                self.map_values,
                self.map_capacity,
                state.neighbors,
                self.mortar_list,
                self.num_mortars,
                self.max_mortars,
                # level removed
                int(self.periodic_x),
                int(self.periodic_y)
            ],
            device=self.device
        )
        
        # 4. Mark Coarse Mortars
        num_mortars_host = int(self.num_mortars.numpy()[0])
        if num_mortars_host > 0:
            wp.launch(
                kernel=mark_mortar_neighbors,
                dim=num_mortars_host,
                inputs=[self.mortar_list, self.num_mortars, state.neighbors],
                device=self.device
            )


    def get_bounds(self, morton_code: int) -> Tuple[float, float, float, float]:
        """
        Decodes a Morton code to get the bounding box of the block.
        """
        if morton_code == 0:
             return self.root_bounds # Root
        
        # Determine level from sentinel bit position
        # The code is (1 << 2*level) | interleaved.
        # So msb_index = 2*level. level = msb_index // 2.
        
        msb_index = morton_code.bit_length() - 1
        level = msb_index // 2
        
        ix, iy = Quadtree.morton_decode(morton_code, level)
        
        grid_dim = 1 << level
        domain_w = self.root_bounds[2] - self.root_bounds[0]
        domain_h = self.root_bounds[3] - self.root_bounds[1]
        
        dx = domain_w / grid_dim
        dy = domain_h / grid_dim
        
        x0 = self.root_bounds[0] + ix * dx
        y0 = self.root_bounds[1] + iy * dy
        
        return (x0, y0, x0 + dx, y0 + dy)

    def find_neighbors(self, morton_code: int):
        """
        Calculates neighbor Morton codes using bitwise operations.
        No explicit connectivity list is stored.
        """
        pass