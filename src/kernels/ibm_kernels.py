import warp as wp
from src.kernels.grid_kernels import morton_encode
from src.kernels.connectivity_kernels import map_lookup
from typing import Any

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
    phi: wp.array(dtype=float, ndim=2),
    x: wp.array(dtype=float, ndim=2),
    y: wp.array(dtype=float, ndim=2),
    nodes_1d: wp.array(dtype=float),
    N1: int,
    params: Any
):
    """
    Interpolates the state q at point (x_p, y_p) within the block specified by pool_idx.
    Uses Renormalized Bilinear Interpolation to exclude solid nodes (Implicit Feedback Loop).
    """
    # 1. Determine Block Bounds via Extrapolation
    
    # Reference coordinates
    xi_start = nodes_1d[0]
    xi_end = nodes_1d[N1 - 1]
    
    # X-direction (using bottom row: 0 to N1-1)
    x_start = x[pool_idx, 0]
    x_end = x[pool_idx, N1 - 1]
    
    # Jacobian dX/dxi
    d_xi = xi_end - xi_start
    if d_xi < 1e-12:
         d_xi = 2.0
    
    J_x = (x_end - x_start) / d_xi
    
    true_x_min = x_start - J_x * (xi_start - (-1.0))
    true_x_max = x_end + J_x * (1.0 - xi_end)
    
    # Y-direction (using left column: 0 to (N1-1)*N1)
    idx_top_left = (N1 - 1) * N1
    y_start = y[pool_idx, 0]
    y_end = y[pool_idx, idx_top_left]
    
    J_y = (y_end - y_start) / d_xi
    
    true_y_min = y_start - J_y * (xi_start - (-1.0))
    true_y_max = y_end + J_y * (1.0 - xi_end)
    
    # 2. Clamp to True Block Boundaries
    xp_c = wp.clamp(x_p, true_x_min, true_x_max)
    yp_c = wp.clamp(y_p, true_y_min, true_y_max)
    
    # 3. Map to Reference Coordinates [-1, 1]
    xi = 2.0 * (xp_c - true_x_min) / (true_x_max - true_x_min) - 1.0
    eta = 2.0 * (yp_c - true_y_min) / (true_y_max - true_y_min) - 1.0
    
    # 4. Find the GLL cell containing (xi, eta)
    idx_i = int(0)
    idx_j = int(0)
    
    # Linear scan for i (x-direction)
    for k in range(N1 - 1):
        if xi >= nodes_1d[k] and xi <= nodes_1d[k+1]:
            idx_i = k
            break
            
    # Linear scan for j (y-direction)
    for k in range(N1 - 1):
        if eta >= nodes_1d[k] and eta <= nodes_1d[k+1]:
            idx_j = k
            break
            
    # 5. Bilinear Interpolation Weights (Base)
    xi_0 = nodes_1d[idx_i]
    xi_1 = nodes_1d[idx_i+1]
    eta_0 = nodes_1d[idx_j]
    eta_1 = nodes_1d[idx_j+1]
    
    u_w = (xi - xi_0) / (xi_1 - xi_0)
    v_w = (eta - eta_0) / (eta_1 - eta_0)
    
    w00 = (1.0 - u_w) * (1.0 - v_w)
    w10 = u_w * (1.0 - v_w)
    w01 = (1.0 - u_w) * v_w
    w11 = u_w * v_w
    
    # 6. Fetch Data and Apply Renormalization Mask
    # Indices in the flattened block array
    node_00 = idx_j * N1 + idx_i
    node_10 = idx_j * N1 + (idx_i + 1)
    node_01 = (idx_j + 1) * N1 + idx_i
    node_11 = (idx_j + 1) * N1 + (idx_i + 1)
    
    # Fetch State
    q00 = q[pool_idx, node_00]
    q10 = q[pool_idx, node_10]
    q01 = q[pool_idx, node_01]
    q11 = q[pool_idx, node_11]
    
    # Fetch SDF (phi)
    phi00 = phi[pool_idx, node_00]
    phi10 = phi[pool_idx, node_10]
    phi01 = phi[pool_idx, node_01]
    phi11 = phi[pool_idx, node_11]
    
    # --- STRICT MASKING to Break Feedback Loop ---
    # Nodes with phi <= 0.0 are FORCED. We must NOT use them.
    # We use a positive epsilon to be safe.
    safe_epsilon = 1e-9
    
    m00 = wp.where(phi00 > safe_epsilon, 1.0, 0.0)
    m10 = wp.where(phi10 > safe_epsilon, 1.0, 0.0)
    m01 = wp.where(phi01 > safe_epsilon, 1.0, 0.0)
    m11 = wp.where(phi11 > safe_epsilon, 1.0, 0.0)
    
    # Effective Weights
    ew00 = w00 * m00
    ew10 = w10 * m10
    ew01 = w01 * m01
    ew11 = w11 * m11
    
    # Renormalize
    w_total = ew00 + ew10 + ew01 + ew11
    
    q_final = wp.vec4(0.0, 0.0, 0.0, 0.0)
    
    # --- SAFE FALLBACK ---
    if w_total > 1e-9:
        inv_w = 1.0 / w_total
        q_final = (q00 * ew00 + q10 * ew10 + q01 * ew01 + q11 * ew11) * inv_w
    else:
        # Fallback to Freestream (Dirichlet Anchor)
        rho_inf = params.rho_inf
        u_inf = params.u_inf
        v_inf = params.v_inf
        p_inf = params.p_inf
        gamma = params.gamma
        
        v2 = u_inf*u_inf + v_inf*v_inf
        E_inf = p_inf / (gamma - 1.0) + 0.5 * rho_inf * v2
        
        q_final = wp.vec4(rho_inf, rho_inf * u_inf, rho_inf * v_inf, E_inf)
            
    return q_final

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
    ix = int((x - x_min) / dx)
    iy = int((y - y_min) / dy)
    
    # Boundary checks
    if ix < 0 or ix >= grid_dim or iy < 0 or iy >= grid_dim:
        return -1
        
    code = morton_encode(ix, iy, level)
    
    return map_lookup(map_keys, map_values, capacity, code)

@wp.kernel
def apply_ibm_forcing(
    q: wp.array(dtype=wp.vec4, ndim=2),
    x: wp.array(dtype=float, ndim=2),
    y: wp.array(dtype=float, ndim=2),
    phi: wp.array(dtype=float, ndim=2),
    active_block_indices: wp.array(dtype=int),
    mesh_id: wp.uint64,
    nodes_1d: wp.array(dtype=float),
    N1: int,
    cutoff_dist: float,
    map_keys: wp.array(dtype=int),
    map_values: wp.array(dtype=int),
    map_capacity: int,
    root_bounds: wp.vec4,
    block_levels: wp.array(dtype=int),
    invert: int,
    boundary_type_id: int,
    params: Any
):
    """
    Applies Ghost-Cell forcing for IBM.
    Reflects the state at 'image point' to the 'ghost node' (solid).
    Includes safety clamping to prevent vacuum/velocity explosions.
    """
    block_id, node_id = wp.tid()
    pool_idx = active_block_indices[block_id]
    
    phi_val = phi[pool_idx, node_id]
    
    # Only force nodes that are strictly solid/interface
    if phi_val < 0.0:
        px = x[pool_idx, node_id]
        py = y[pool_idx, node_id]
        p_vec = wp.vec3(px, py, 0.0)
        
        face_index = int(0)
        face_u = float(0.0)
        face_v = float(0.0)
        sign = float(0.0)
        
        if wp.mesh_query_point_sign_normal(mesh_id, p_vec, cutoff_dist, sign, face_index, face_u, face_v):
            p_surf = wp.mesh_eval_position(mesh_id, face_index, face_u, face_v)
            n_surf = wp.mesh_eval_face_normal(mesh_id, face_index)
            
            diff = p_vec - p_surf
            d = wp.length(diff)
            
            nx = n_surf[0]
            ny = n_surf[1]
            if invert == 1:
                nx = -nx
                ny = -ny

            inv_len = 1.0 / wp.sqrt(nx*nx + ny*ny + 1e-9)
            nx *= inv_len
            ny *= inv_len
            
            x_img = px + 2.0 * d * nx
            y_img = py + 2.0 * d * ny
            
            current_level = block_levels[pool_idx]
            target_block_idx = find_block_id(x_img, y_img, current_level, map_keys, map_values, map_capacity, root_bounds)
            
            if target_block_idx == -1:
                target_block_idx = pool_idx
            
            # Interpolate
            q_img = sample_state_at_point(target_block_idx, x_img, y_img, q, phi, x, y, nodes_1d, N1, params)
            
            rho_img = q_img[0]
            rhou_img = q_img[1]
            rhov_img = q_img[2]
            E_img = q_img[3]
            
            # --- START FIX: Sanitize Image State ---
            
            # 1. Clamp Density to Physics Floor
            rho_safe = wp.max(rho_img, params.rho_floor)
            
            u_img = float(0.0)
            v_img = float(0.0)
            
            # 2. Extract Velocity Safely
            # If density is dangerously close to the floor (vacuum), 
            # we zero out the velocity to prevent numerical explosion.
            # 2.0x safety margin is arbitrary but robust.
            if rho_img > 2.0 * params.rho_floor:
                inv_rho = 1.0 / rho_img
                u_img = rhou_img * inv_rho
                v_img = rhov_img * inv_rho
            
            # 3. Calculate Image Pressure (safely)
            ke_img = 0.5 * rho_safe * (u_img*u_img + v_img*v_img)
            p_img = (params.gamma - 1.0) * (E_img - ke_img)
            
            # 4. Clamp Pressure
            p_safe = wp.max(p_img, params.p_floor)
            
            # --- End FIX ---

            # Apply Boundary Condition (Reflection)
            if boundary_type_id == 1: # No-Slip
                u_ghost = -u_img
                v_ghost = -v_img
            else: # Slip
                v_dot_n = u_img * nx + v_img * ny
                vn_x = v_dot_n * nx
                vn_y = v_dot_n * ny
                vt_x = u_img - vn_x
                vt_y = v_img - vn_y
                u_ghost = vt_x - vn_x
                v_ghost = vt_y - vn_y
            
            # Reconstruct Ghost State
            # We copy density/pressure from image (Neumann-ish for P, Dirichlet for rho)
            # but use the Reflected velocity.
            
            rho_ghost = rho_safe 
            ke_ghost = 0.5 * rho_ghost * (u_ghost*u_ghost + v_ghost*v_ghost)
            E_ghost = p_safe / (params.gamma - 1.0) + ke_ghost
            
            q_ghost = wp.vec4(rho_ghost, rho_ghost * u_ghost, rho_ghost * v_ghost, E_ghost)
            q[pool_idx, node_id] = q_ghost