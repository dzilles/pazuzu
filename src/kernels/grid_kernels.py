import warp as wp
from typing import Any
from src.kernels import boundary_conditions as bc

@wp.func
def part1by1(n: int):
    """Inserts a 0 bit after each of the low 16 bits of n."""
    n = (n ^ (n << 8)) & 0x00ff00ff
    n = (n ^ (n << 4)) & 0x0f0f0f0f
    n = (n ^ (n << 2)) & 0x33333333
    n = (n ^ (n << 1)) & 0x55555555
    return n

@wp.func
def morton_encode(x: int, y: int):
    """Interleaves bits of x and y (Z-order curve)."""
    return part1by1(x) | (part1by1(y) << 1)

@wp.func
def compact1by1(n: int):
    """Inverse of part1by1."""
    n = n & 0x55555555
    n = (n ^ (n >> 1)) & 0x33333333
    n = (n ^ (n >> 2)) & 0x0f0f0f0f
    n = (n ^ (n >> 4)) & 0x00ff00ff
    n = (n ^ (n >> 8)) & 0x0000ffff
    return n

@wp.func
def morton_decode(code: int):
    """Decodes Morton code into (x, y)."""
    return compact1by1(code), compact1by1(code >> 1)

@wp.kernel
def generate_morton_codes(
    codes: Any,
    grid_dim: int
):
    tid = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    iy = tid // grid_dim  # type: ignore # Warp type inference
    ix = tid % grid_dim  # type: ignore # Warp type inference
    codes[tid] = morton_encode(ix, iy)  # type: ignore # Warp type inference

@wp.kernel
def compute_block_coordinates(
    morton_codes: Any,
    active_indices: Any,
    num_active: int,
    block_levels: Any,
    root_bounds: Any, 
    nodes_2d: Any,
    out_x: Any,
    out_y: Any
):
    """
    Computes physical coordinates for all nodes in active blocks.
    
    Args:
        morton_codes: Array of Morton codes (indexed by pool_idx).
        active_indices: Array of active block pool indices.
        num_active: Number of active blocks.
        block_levels: Array of refinement levels (indexed by pool_idx).
        root_bounds: Domain boundaries.
        nodes_2d: Reference element nodes in [-1, 1].
        out_x: Output X coordinates (max_blocks, Np).
        out_y: Output Y coordinates (max_blocks, Np).
    """
    tid_block, tid_node = wp.tid()
    
    if tid_block >= num_active:
        return
        
    pool_idx = active_indices[tid_block]
    code = morton_codes[pool_idx]
    level = block_levels[pool_idx]
    
    ix, iy = morton_decode(code)
    
    # Grid dimensions at this level
    grid_dim = 1 << level
    
    # Map reference node [-1, 1] to physical block
    ref_node = nodes_2d[tid_node]
    r = ref_node[0]
    s = ref_node[1]

    # Domain size
    domain_w = root_bounds[2] - root_bounds[0]
    domain_h = root_bounds[3] - root_bounds[1]
    
    # Block size
    f_grid_dim = bc.get_any_generic(r, grid_dim)
    dx = domain_w / f_grid_dim
    dy = domain_h / f_grid_dim
    
    # Block origin (bottom-left)
    x0 = root_bounds[0] + bc.get_any_generic(r, ix) * dx
    y0 = root_bounds[1] + bc.get_any_generic(r, iy) * dy
    
    # x = x0 + (r + 1)/2 * dx
    # y = y0 + (s + 1)/2 * dy
    one = bc.get_one_generic(r)
    half = bc.get_half_generic(r)
    phys_x = x0 + (r + one) * half * dx
    phys_y = y0 + (s + one) * half * dy
    
    out_x[pool_idx, tid_node] = phys_x
    out_y[pool_idx, tid_node] = phys_y
