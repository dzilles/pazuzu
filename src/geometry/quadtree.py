import warp as wp
import numpy as np
from typing import Tuple, Optional
from src.kernels.grid_kernels import compute_block_coordinates
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
    def __init__(self, device: str = "cpu", max_blocks: int = 10000, root_bounds: Tuple[float, float, float, float] = (-1.0, -1.0, 1.0, 1.0), periodic_x: bool = False, periodic_y: bool = False):
        self.device = device
        self.max_blocks = max_blocks
        self.num_blocks = 0
        self.root_bounds = root_bounds
        self.periodic_x = periodic_x
        self.periodic_y = periodic_y
        
        # Morton codes for blocks resident in the pool.
        # Index i in this array corresponds to block i in SimulationState.
        self.block_morton_codes = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        self.block_levels = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        
        # Hash Map for Morton -> Pool Index Lookups
        # Capacity should be larger than max_blocks to reduce collisions (e.g., 2x)
        self.map_capacity = max_blocks * 2
        self.map_keys = wp.zeros(self.map_capacity, dtype=wp.int32, device=device)
        self.map_values = wp.zeros(self.map_capacity, dtype=wp.int32, device=device)
        
        # Bounds as Warp vector
        self.root_bounds_wp = wp.vec4(root_bounds[0], root_bounds[1], root_bounds[2], root_bounds[3])

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
        
        # 1. Generate Morton Codes on Host
        host_codes = np.zeros(num_blocks, dtype=np.int32)
        idx = 0
        # Z-order iteration
        for iy in range(grid_dim):
            for ix in range(grid_dim):
                code = morton_encode(ix, iy)
                host_codes[idx] = code
                idx += 1
        
        # 2. Upload to Device
        # We fill the first num_blocks slots of the pool
        wp.copy(self.block_morton_codes, wp.array(host_codes, dtype=wp.int32, device=self.device), count=num_blocks)
        
        # Set levels (all same)
        self.block_levels.fill_(level) 
        
        # 3. Update Active Block Indices in State
        # In a compact pool, the active indices are just 0, 1, ..., num_blocks-1
        active_indices = np.arange(num_blocks, dtype=np.int32)
        wp.copy(state.active_block_indices, wp.array(active_indices, dtype=wp.int32, device=self.device), count=num_blocks)
        
        # 4. Compute Physical Coordinates
        wp.launch(
            kernel=compute_block_coordinates,
            dim=(num_blocks, basis.Np),
            inputs=[
                self.block_morton_codes,
                level,
                self.root_bounds_wp,
                basis.nodes_2d,
                state.x,
                state.y
            ],
            device=self.device
        )
        
        # 5. Build Connectivity
        self.build_connectivity(state, level)
        
        print(f"Created {self.num_blocks} blocks.")

    def build_connectivity(self, state, level: int):
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
                state.active_block_indices,
                self.num_blocks,
                self.map_keys,
                self.map_values,
                self.map_capacity
            ],
            device=self.device
        )
        
        # 3. Compute Neighbors
        # Reset neighbors first? Not strictly necessary if we overwrite all, 
        # but good practice if logic is partial. Here we overwrite.
        
        wp.launch(
            kernel=compute_neighbors,
            dim=self.num_blocks,
            inputs=[
                self.block_morton_codes,
                state.active_block_indices,
                self.num_blocks,
                self.map_keys,
                self.map_values,
                self.map_capacity,
                state.neighbors,
                level,
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