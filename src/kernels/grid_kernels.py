import warp as wp

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
def compute_block_coordinates(
    morton_codes: wp.array(dtype=wp.int32),
    level: int,
    root_bounds: wp.vec4, # min_x, min_y, max_x, max_y
    nodes_2d: wp.array(dtype=wp.vec2),
    out_x: wp.array(dtype=float, ndim=2),
    out_y: wp.array(dtype=float, ndim=2)
):
    """
    Computes physical coordinates for all nodes in active blocks.
    
    Args:
        morton_codes: Array of Morton codes for active blocks.
        level: Refinement level of the blocks (uniform grid assumption for this kernel).
        root_bounds: Domain boundaries.
        nodes_2d: Reference element nodes in [-1, 1].
        out_x: Output X coordinates (num_blocks, Np).
        out_y: Output Y coordinates (num_blocks, Np).
    """
    block_idx, node_idx = wp.tid()
    
    code = morton_codes[block_idx]
    ix, iy = morton_decode(code)
    
    # Grid dimensions at this level
    grid_dim = 1 << level
    
    # Domain size
    domain_w = root_bounds[2] - root_bounds[0]
    domain_h = root_bounds[3] - root_bounds[1]
    
    # Block size
    dx = domain_w / float(grid_dim)
    dy = domain_h / float(grid_dim)
    
    # Block origin (bottom-left)
    x0 = root_bounds[0] + float(ix) * dx
    y0 = root_bounds[1] + float(iy) * dy
    
    # Map reference node [-1, 1] to physical block
    ref_node = nodes_2d[node_idx]
    r = ref_node[0]
    s = ref_node[1]
    
    # x = x0 + (r + 1)/2 * dx
    # y = y0 + (s + 1)/2 * dy
    phys_x = x0 + (r + 1.0) * 0.5 * dx
    phys_y = y0 + (s + 1.0) * 0.5 * dy
    
    out_x[block_idx, node_idx] = phys_x
    out_y[block_idx, node_idx] = phys_y
