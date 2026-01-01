import warp as wp
import numpy as np
from typing import Tuple, Optional
from src.kernels.grid_kernels import compute_block_coordinates, generate_morton_codes
from src.kernels.connectivity_kernels import init_hash_map, populate_hash_map, compute_neighbors

MAX_DEPTH = 10
ROOT_BOUNDS = (-1.0, -1.0, 1.0, 1.0) # x_min, y_min, x_max, y_max

# --- Python-side Morton Encoding Helpers ---
def part1by1(n: int) -> int:
    """Inserts a 0 bit after each of the low 16 bits of n."""
    n = (n ^ (n << 8)) & 0x00ff00ff
    n = (n ^ (n << 4)) & 0x0f0f0f0f
    n = (n ^ (n << 2)) & 0x33333333
    n = (n ^ (n << 1)) & 0x55555555
    return n

def morton_encode(x: int, y: int) -> int:
    """Interleaves bits of x and y (Z-order curve)."""
    return part1by1(x) | (part1by1(y) << 1)

class Quadtree:
    """
    Block-AMR Quadtree using Morton Encoding (Z-ordering) for implicit connectivity.
    
    Manages the grid hierarchy and mapping to the memory pool.
    """
    def __init__(self, device: str = "cpu", max_blocks: int = 10000, root_bounds: Tuple[float, float, float, float] = (-1.0, -1.0, 1.0, 1.0), periodic_x: bool = False, periodic_y: bool = False, dtype=wp.float32):
        self.device = device
        self.max_blocks = max_blocks
        self.num_blocks = 0
        self.root_bounds = root_bounds
        self.periodic_x = periodic_x
        self.periodic_y = periodic_y
        self.dtype = dtype
        
        # Morton codes for blocks resident in the pool.
        # Index i in this array corresponds to block i in SimulationState.
        self.block_morton_codes = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        self.block_levels = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        
        # Free List Management (Stack of free indices)
        # Initialize with all indices 0..max_blocks-1
        self.free_pool_indices = wp.array(
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

    def uniform_refine(self, level: int, state, basis):
        """
        Generates a uniform grid at the specified refinement level.
        
        Args:
            level (int): Refinement level (0 = root).
            state (SimulationState): The simulation state to populate.
            basis (Basis): The basis defining the nodes.
        """
        if level > MAX_DEPTH:
            raise ValueError(f"Level {level} exceeds MAX_DEPTH {MAX_DEPTH}")
            
        print(f"Generating uniform grid at level {level}...")
        
        grid_dim = 1 << level
        num_blocks = grid_dim * grid_dim
        
        if num_blocks > self.max_blocks:
             raise ValueError(f"Level {level} requires {num_blocks} blocks, but MAX_BLOCKS is {self.max_blocks}")

        self.num_blocks = num_blocks
        self.num_free = self.max_blocks - num_blocks # Remaining free blocks
        
        # 1. Generate Morton Codes on Device
        # Note: We simply take the first num_blocks indices for a uniform grid
        wp.launch(
            kernel=generate_morton_codes,
            dim=num_blocks,
            inputs=[self.block_morton_codes, grid_dim],
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

    def refine_blocks(self, blocks_to_refine_indices: list, state, basis):
        """
        Refines the specified active blocks by splitting them into 4 children.
        
        Args:
            blocks_to_refine_indices (list[int]): List of pool indices of blocks to refine.
            state (SimulationState): The simulation state.
            basis (Basis): The basis object.
        """
        num_refine = len(blocks_to_refine_indices)
        if num_refine == 0:
            return

        # Check free space: each split consumes 4 new slots and frees 1 slot (net +3)
        # But conceptually we consume 4 from free list and the parent becomes inactive (or returned to free list).
        # We will pop 4 for children, and push parent back later (or just swap 1 child into parent slot? 
        # Simpler to be explicit: pop 4, push 1).
        
        needed = num_refine * 4
        if self.num_free < needed:
            raise RuntimeError(f"Not enough free blocks for refinement. Needed {needed}, have {self.num_free}")

        # 1. Fetch data to Host for logic processing
        # We need parent codes and levels.
        # Ideally we'd do this on GPU, but topology changes are complex.
        
        # Get active indices from state (host copy)
        active_indices_host = state.active_block_indices.numpy()[:self.num_blocks]
        
        # Create a set for fast lookup/removal
        active_set = set(active_indices_host)
        
        # Get parent info
        # We need to read from the device arrays at specific indices
        # Optimization: Read all or copy slice? 
        # For now, let's copy the whole array to host, modify, copy back. (Slow but robust for Phase 2)
        
        h_morton_codes = self.block_morton_codes.numpy()
        h_levels = self.block_levels.numpy()
        h_free_indices = self.free_pool_indices.numpy()
        
        # Stack pointer for free list
        free_ptr = 0 # reading from 0..needed-1
        
        # New active list construction
        new_active_list = []
        
        # Mark parents for removal
        parents_to_remove = set(blocks_to_refine_indices)
        
        # Add non-refined blocks to new list
        for idx in active_indices_host:
            if idx not in parents_to_remove:
                new_active_list.append(idx)
                
        # Process Refinement
        for p_idx in blocks_to_refine_indices:
            p_code = h_morton_codes[p_idx]
            p_level = h_levels[p_idx]
            
            # Check max level
            if p_level >= MAX_DEPTH:
                # Can't refine, just keep parent? Or raise error?
                # For now, keep parent.
                new_active_list.append(p_idx)
                continue
                
            # Get 4 children indices from free pool
            c_indices = h_free_indices[free_ptr : free_ptr + 4]
            free_ptr += 4
            
            # Compute children codes
            # (2x, 2y), (2x+1, 2y), ...
            # Child codes are (p_code << 2) | i
            # i=0 (00): 2x, 2y
            # i=1 (01): 2x+1, 2y
            # i=2 (10): 2x, 2y+1
            # i=3 (11): 2x+1, 2y+1
            
            for i in range(4):
                c_idx = c_indices[i]
                c_code = (p_code << 2) | i
                c_level = p_level + 1
                
                h_morton_codes[c_idx] = c_code
                h_levels[c_idx] = c_level
                
                new_active_list.append(c_idx)
        
        # Update State on Device
        
        # Update Morton Codes and Levels
        # (We modified h_morton_codes and h_levels in place)
        self.block_morton_codes = wp.array(h_morton_codes, dtype=wp.int32, device=self.device)
        self.block_levels = wp.array(h_levels, dtype=wp.int32, device=self.device)
        
        # Update Free List
        # Shift remaining free down? Or just update the view?
        # We used `free_ptr` items.
        # Also we need to return `parents_to_remove` to the free list.
        # Update h_free_indices: shift remaining to front, append parents at end.
        
        remaining_free = h_free_indices[free_ptr : self.num_free]
        freed_parents = list(parents_to_remove)
        
        new_free_count = len(remaining_free) + len(freed_parents)
        h_new_free = np.zeros(self.max_blocks, dtype=np.int32) # or reuse buffer
        
        h_new_free[:len(remaining_free)] = remaining_free
        h_new_free[len(remaining_free):new_free_count] = freed_parents
        # Fill rest with junk or old values (doesn't matter)
        
        self.free_pool_indices = wp.array(h_new_free, dtype=wp.int32, device=self.device)
        self.num_free = new_free_count
        
        # Update Active Indices
        self.num_blocks = len(new_active_list)
        new_active_array = np.array(new_active_list, dtype=np.int32)
        wp.copy(state.active_block_indices, wp.array(new_active_array, dtype=wp.int32, device=self.device), count=self.num_blocks)
        
        # Recompute Coordinates for NEW blocks
        # Actually, we need to recompute for ALL blocks or just new ones?
        # Simpler to recompute all active blocks to be safe, or optimize to only compute children.
        # Optimization: Only compute for new children.
        # But `compute_block_coordinates` takes indices.
        # Let's just run it for all active blocks for now to ensure consistency. 
        # (Coordinate computation is cheap).
        
        wp.launch(
            kernel=compute_block_coordinates,
            dim=(self.num_blocks, basis.Np),
            inputs=[
                self.block_morton_codes,
                state.active_block_indices, # Use active indices!
                self.num_blocks,            # count
                self.block_levels,          # Now passed as array! Wait, kernel needs update?
                self.root_bounds_wp,
                basis.nodes_2d,
                state.x,
                state.y
            ],
            device=self.device
        )
        
        # Rebuild Connectivity
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
                # level removed
                int(self.periodic_x),
                int(self.periodic_y)
            ],
            device=self.device
        )


    def get_bounds(self, morton_code: int) -> Tuple[float, float, float, float]:

        """
        Decodes a Morton code to get the bounding box of the block.
        """
        # Placeholder
        return self.root_bounds

    def find_neighbors(self, morton_code: int):
        """
        Calculates neighbor Morton codes using bitwise operations.
        No explicit connectivity list is stored.
        """
        pass