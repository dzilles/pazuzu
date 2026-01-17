import warp as wp
from src.kernels.grid_kernels import morton_encode
from src.kernels.connectivity_kernels import map_lookup

@wp.kernel
def generate_sdf_cylinder(
    x: wp.array(dtype=float, ndim=2),
    y: wp.array(dtype=float, ndim=2),
    phi_out: wp.array(dtype=float, ndim=2),
    active_block_indices: wp.array(dtype=int),
    center_x: float,
    center_y: float,
    radius: float
):
    """
    Computes the Analytical Signed Distance Field (SDF) for a 2D cylinder.
    phi = sqrt((x-cx)^2 + (y-cy)^2) - r
    """
    block_id, node_id = wp.tid()
    pool_idx = active_block_indices[block_id]
    
    px = x[pool_idx, node_id]
    py = y[pool_idx, node_id]
    
    dx = px - center_x
    dy = py - center_y
    
    dist = wp.sqrt(dx*dx + dy*dy)
    phi_out[pool_idx, node_id] = dist - radius

@wp.kernel
def generate_sdf_from_mesh(
    x: wp.array(dtype=float, ndim=2),
    y: wp.array(dtype=float, ndim=2),
    phi_out: wp.array(dtype=float, ndim=2),
    active_block_indices: wp.array(dtype=int),
    mesh_id: wp.uint64,
    max_dist: float,
    invert: int
):
    """
    Computes the Signed Distance Field (SDF) using a Warp Mesh (BVH) query.
    """
    block_id, node_id = wp.tid()
    pool_idx = active_block_indices[block_id]
    
    px = x[pool_idx, node_id]
    py = y[pool_idx, node_id]
    
    # Construct 3D query point (assuming 2D simulation is on Z=0 plane)
    p = wp.vec3(px, py, 0.0)
    
    # Variables to store query results
    face_index = int(0)
    face_u = float(0.0)
    face_v = float(0.0)
    sign = float(0.0)
    
    if wp.mesh_query_point_sign_normal(
        mesh_id,
        p,
        max_dist,
        sign,
        face_index,
        face_u,
        face_v
    ):
        p_mesh = wp.mesh_eval_position(mesh_id, face_index, face_u, face_v)
        dist = wp.length(p - p_mesh)
        
        # Final signed distance
        res = dist * sign
        
        if invert == 1:
            res = -res
            
        phi_out[pool_idx, node_id] = res
    else:
        phi_out[pool_idx, node_id] = max_dist

@wp.func
def sample_state_at_point(
    pool_idx: int,
    x_p: float,
    y_p: float,
    q: wp.array(dtype=wp.vec4, ndim=2),
    x: wp.array(dtype=float, ndim=2),
    y: wp.array(dtype=float, ndim=2),
    nodes_1d: wp.array(dtype=float),
    N1: int # Number of nodes in 1D (N+1)
):
    """
    Interpolates the state q at point (x_p, y_p) within the block specified by pool_idx.
    Uses bilinear interpolation on the GLL sub-grid.
    """
    # 1. Determine Block Bounds
    # Assuming tensor product ordering: 0 is bottom-left, Np-1 is top-right
    Np = N1 * N1
    
    # Bounds of the block
    x_min = x[pool_idx, 0]
    x_max = x[pool_idx, Np - 1]
    y_min = y[pool_idx, 0]
    y_max = y[pool_idx, Np - 1]
    
    # 2. Clamp to Block Boundaries (Nearest Neighbor Fallback for ghost cells outside)
    xp_c = wp.clamp(x_p, x_min, x_max)
    yp_c = wp.clamp(y_p, y_min, y_max)
    
    # 3. Map to Reference Coordinates [-1, 1]
    # xi = 2 * (x - xmin) / (xmax - xmin) - 1
    xi = 2.0 * (xp_c - x_min) / (x_max - x_min) - 1.0
    eta = 2.0 * (yp_c - y_min) / (y_max - y_min) - 1.0
    
    # 4. Find the GLL cell containing (xi, eta)
    # nodes_1d is sorted [-1, ... 1]
    # We find i such that nodes_1d[i] <= xi <= nodes_1d[i+1]
    
    idx_i = int(0)
    idx_j = int(0)
    
    # Linear scan for i (x-direction)
    # Range is 0 to N1-2 (since we check i and i+1)
    for k in range(N1 - 1):
        if xi >= nodes_1d[k] and xi <= nodes_1d[k+1]:
            idx_i = k
            break
            
    # Linear scan for j (y-direction)
    for k in range(N1 - 1):
        if eta >= nodes_1d[k] and eta <= nodes_1d[k+1]:
            idx_j = k
            break
            
    # 5. Bilinear Interpolation Weights
    # u_local in [0, 1] within the interval
    xi_0 = nodes_1d[idx_i]
    xi_1 = nodes_1d[idx_i+1]
    eta_0 = nodes_1d[idx_j]
    eta_1 = nodes_1d[idx_j+1]
    
    u_w = (xi - xi_0) / (xi_1 - xi_0)
    v_w = (eta - eta_0) / (eta_1 - eta_0)
    
    # 6. Fetch 4 corner values
    # Indices in the flattened block array
    # idx = j * N1 + i
    node_00 = idx_j * N1 + idx_i
    node_10 = idx_j * N1 + (idx_i + 1)
    node_01 = (idx_j + 1) * N1 + idx_i
    node_11 = (idx_j + 1) * N1 + (idx_i + 1)
    
    q00 = q[pool_idx, node_00]
    q10 = q[pool_idx, node_10]
    q01 = q[pool_idx, node_01]
    q11 = q[pool_idx, node_11]
    
    # Interpolate
    # (1-u)(1-v) * q00 + u(1-v) * q10 + (1-u)v * q01 + uv * q11
    q_val = (1.0 - u_w) * (1.0 - v_w) * q00 + \
            u_w * (1.0 - v_w) * q10 + \
            (1.0 - u_w) * v_w * q01 + \
            u_w * v_w * q11
            
    return q_val

@wp.func
def find_block_id(
    x: float,
    y: float,
    level: int,
    map_keys: wp.array(dtype=int),
    map_values: wp.array(dtype=int),
    capacity: int,
    root_bounds: wp.vec4
):
    # Unwrap bounds
    x_min = root_bounds[0]
    y_min = root_bounds[1]
    x_max = root_bounds[2]
    y_max = root_bounds[3]
    
    grid_dim = 1 << level
    
    # Calculate cell size
    domain_w = x_max - x_min
    domain_h = y_max - y_min
    
    dx = domain_w / float(grid_dim)
    dy = domain_h / float(grid_dim)
    
    # Calculate integer coordinates
    # We add a small epsilon to handle boundary cases
    ix = int((x - x_min) / dx)
    iy = int((y - y_min) / dy)
    
    # Boundary checks
    if ix < 0 or ix >= grid_dim or iy < 0 or iy >= grid_dim:
        return -1
        
    code = morton_encode(ix, iy, level)
    
    return map_lookup(map_keys, map_values, capacity, code)

@wp.kernel
def apply_ibm_forcing(
    q: wp.array(dtype=wp.vec4, ndim=2),       # (max_blocks, Np)
    x: wp.array(dtype=float, ndim=2),         # (max_blocks, Np)
    y: wp.array(dtype=float, ndim=2),         # (max_blocks, Np)
    phi: wp.array(dtype=float, ndim=2),       # (max_blocks, Np)
    active_block_indices: wp.array(dtype=int),# (num_active_blocks)
    mesh_id: wp.uint64,
    nodes_1d: wp.array(dtype=float),          # (N1)
    N1: int,
    cutoff_dist: float,                        # e.g., 2.0 * dx_min
    # New Arguments
    map_keys: wp.array(dtype=int),
    map_values: wp.array(dtype=int),
    map_capacity: int,
    root_bounds: wp.vec4,
    block_levels: wp.array(dtype=int),
    invert: int,
    boundary_type_id: int # 0: Slip, 1: No-Slip
):
    """
    Applies Ghost-Cell forcing for IBM.
    Reflects the state at 'image point' to the 'ghost node' (solid).
    """
    block_id, node_id = wp.tid()
    
    # Map to global pool index
    pool_idx = active_block_indices[block_id]
    
    # 1. Check Ghost Node Condition
    # Node is Solid (phi < 0)
    phi_val = phi[pool_idx, node_id]
    
    # Force ALL solid nodes to prevent polynomial ringing
    if phi_val < 0.0:
        
        # 2. Geometry Query
        px = x[pool_idx, node_id]
        py = y[pool_idx, node_id]
        p_vec = wp.vec3(px, py, 0.0)
        
        face_index = int(0)
        face_u = float(0.0)
        face_v = float(0.0)
        sign = float(0.0)
        
        # Use the cutoff_dist passed from solver (which should be large)
        if wp.mesh_query_point_sign_normal(mesh_id, p_vec, cutoff_dist, sign, face_index, face_u, face_v):
            
            p_surf = wp.mesh_eval_position(mesh_id, face_index, face_u, face_v)
            n_surf = wp.mesh_eval_face_normal(mesh_id, face_index)
            
            # Distance d = |p - p_surf|
            diff = p_vec - p_surf
            d = wp.length(diff)
            
            # Normalize normal
            nx = n_surf[0]
            ny = n_surf[1]
            
            if invert == 1:
                nx = -nx
                ny = -ny

            inv_len = 1.0 / wp.sqrt(nx*nx + ny*ny + 1e-9)
            nx *= inv_len
            ny *= inv_len
            
            # 3. Image Point
            x_img = px + 2.0 * d * nx
            y_img = py + 2.0 * d * ny
            
            # 4. Lookup Block ID for Image Point
            current_level = block_levels[pool_idx]
            target_block_idx = find_block_id(x_img, y_img, current_level, map_keys, map_values, map_capacity, root_bounds)
            
            # Fallback to current block if not found
            if target_block_idx == -1:
                target_block_idx = pool_idx
            
            # 5. Interpolate State at Image Point
            q_img = sample_state_at_point(target_block_idx, x_img, y_img, q, x, y, nodes_1d, N1)
            
            # 6. Apply Boundary Condition (Slip Wall)
            rho_img = q_img[0]
            rhou_img = q_img[1]
            rhov_img = q_img[2]
            E_img = q_img[3]
            
            u_img = rhou_img / rho_img
            v_img = rhov_img / rho_img
            
            # Velocity Reflection
            if boundary_type_id == 1: # No-Slip
                # V_ghost = -V_image
                u_ghost = -u_img
                v_ghost = -v_img
            else: # Slip (Default)
                v_dot_n = u_img * nx + v_img * ny
                vn_x = v_dot_n * nx
                vn_y = v_dot_n * ny
                
                vt_x = u_img - vn_x
                vt_y = v_img - vn_y
                
                # V_ghost = Vt - Vn
                u_ghost = vt_x - vn_x
                v_ghost = vt_y - vn_y
            
            # Construct Ghost State
            q_ghost = wp.vec4(
                rho_img,
                rho_img * u_ghost,
                rho_img * v_ghost,
                E_img
            )
            
            # 7. Write Back
            q[pool_idx, node_id] = q_ghost