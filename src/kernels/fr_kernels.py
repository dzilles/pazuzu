import warp as wp
from src.physics.laws import euler
from src.kernels import boundary_conditions as bc

@wp.kernel
def compute_fr_update(
    q: wp.array(dtype=wp.vec4, ndim=2),
    active_indices: wp.array(dtype=int),
    neighbors: wp.array(dtype=int, ndim=2),
    num_active: int,
    rhs: wp.array(dtype=wp.vec4, ndim=2),
    # Geometry
    nodes_1d: wp.array(dtype=float),
    D1D: wp.array(dtype=float, ndim=2),
    dg_L: wp.array(dtype=float),
    dg_R: wp.array(dtype=float),
    root_bounds: wp.vec4,
    level: int,
    # Physics
    params: euler.EquationParams32, 
    t: float
):
    # 1D Thread Index -> Block + Node
    tid = wp.tid()
    N1 = nodes_1d.shape[0] # N+1 points
    Np = N1 * N1
    
    if tid >= num_active * Np:
        return
        
    # Decompose Index
    block_offset = tid // Np
    node_local = tid % Np
    pool_idx = active_indices[block_offset]
    
    # 2D Node Index (i, j)
    j = node_local // N1
    i = node_local % N1
    
    # Grid Spacing (Uniform Cartesian)
    grid_dim = 1 << level
    domain_w = root_bounds[2] - root_bounds[0]
    domain_h = root_bounds[3] - root_bounds[1]
    dx = domain_w / float(grid_dim)
    dy = domain_h / float(grid_dim)
    
    # Geometric Factors: d/dx = (2/dx) * d/dr
    inv_J_x = 2.0 / dx
    inv_J_y = 2.0 / dy
    
    # --- 1. Volume Gradient (Divergence) ---
    val_div = wp.vec4(0.0, 0.0, 0.0, 0.0)
    
    # Loop over k (1D line)
    for k in range(N1):
        # dF/dx terms (along i-line)
        idx_k_x = j * N1 + k
        q_k_x = q[pool_idx, idx_k_x]
        F_k = euler.flux_x(q_k_x, params)
        val_div += F_k * D1D[i, k] * inv_J_x
        
        # dG/dy terms (along j-line)
        idx_k_y = k * N1 + i
        q_k_y = q[pool_idx, idx_k_y]
        G_k = euler.flux_y(q_k_y, params)
        val_div += G_k * D1D[j, k] * inv_J_y
        
    # --- 2. Interface Corrections (X) ---
    # Left Interface (i=0)
    idx_L = j * N1 + 0 
    q_L_internal = q[pool_idx, idx_L]
    f_L_internal = euler.flux_x(q_L_internal, params)
    
    neigh_L = neighbors[pool_idx, 0]
    q_L_ghost = q_L_internal # Default (Transmissive)
    if neigh_L != -1:
        # Neighbor's Right boundary (i=N1-1)
        idx_neigh = j * N1 + (N1 - 1)
        q_L_ghost = q[neigh_L, idx_neigh]
        
    F_star_L = euler.rusanov_flux(q_L_ghost, q_L_internal, 1.0, 0.0, params)
    corr_x = (F_star_L - f_L_internal) * dg_L[i]
    
    # Right Interface (i=N1-1)
    idx_R = j * N1 + (N1 - 1)
    q_R_internal = q[pool_idx, idx_R]
    f_R_internal = euler.flux_x(q_R_internal, params)
    
    neigh_R = neighbors[pool_idx, 1]
    q_R_ghost = q_R_internal
    if neigh_R != -1:
        # Neighbor's Left boundary (i=0)
        idx_neigh = j * N1 + 0
        q_R_ghost = q[neigh_R, idx_neigh]
        
    F_star_R = euler.rusanov_flux(q_R_internal, q_R_ghost, 1.0, 0.0, params)
    corr_x += (F_star_R - f_R_internal) * dg_R[i]
    
    corr_x *= inv_J_x
    
    # --- 3. Interface Corrections (Y) ---
    # Bottom Interface (j=0)
    idx_B = 0 * N1 + i
    q_B_internal = q[pool_idx, idx_B]
    g_B_internal = euler.flux_y(q_B_internal, params)
    
    neigh_B = neighbors[pool_idx, 2]
    q_B_ghost = q_B_internal
    if neigh_B != -1:
        # Neighbor's Top boundary (j=N1-1)
        idx_neigh = (N1 - 1) * N1 + i
        q_B_ghost = q[neigh_B, idx_neigh]
        
    # Flux in Y direction, Normal=(0,1)
    G_star_B = euler.rusanov_flux(q_B_ghost, q_B_internal, 0.0, 1.0, params)
    corr_y = (G_star_B - g_B_internal) * dg_L[j] # Uses Left correction poly for Bottom (-1)
    
    # Top Interface (j=N1-1)
    idx_T = (N1 - 1) * N1 + i
    q_T_internal = q[pool_idx, idx_T]
    g_T_internal = euler.flux_y(q_T_internal, params)
    
    neigh_T = neighbors[pool_idx, 3]
    q_T_ghost = q_T_internal
    if neigh_T != -1:
        # Neighbor's Bottom boundary (j=0)
        idx_neigh = 0 * N1 + i
        q_T_ghost = q[neigh_T, idx_neigh]
        
    G_star_T = euler.rusanov_flux(q_T_internal, q_T_ghost, 0.0, 1.0, params)
    corr_y += (G_star_T - g_T_internal) * dg_R[j] # Uses Right correction poly for Top (+1)
    
    corr_y *= inv_J_y
    
    # Total Update
    rhs[pool_idx, node_local] = -(val_div + corr_x + corr_y)
