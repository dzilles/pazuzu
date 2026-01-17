import warp as wp
from src.physics.laws import euler
from src.numerics import riemann_solvers as rs
from src.kernels import utils as u
from src.kernels import boundary_conditions as bc

from typing import Any

# Face Indices (B-R-T-L)
FACE_BOTTOM = 0
FACE_RIGHT = 1
FACE_TOP = 2
FACE_LEFT = 3
MORTAR_FLAG = -2

@wp.func
def compute_interface_flux(q_L: Any, q_R: Any, nx: Any, ny: Any, params: Any):
    """
    Computes the numerical flux at an interface using the selected Riemann solver.
    """
    if params.flux_type == 1:
        return rs.hllc_flux(q_L, q_R, nx, ny, params)
    return rs.rusanov_flux(q_L, q_R, nx, ny, params)

@wp.kernel
def compute_fr_update(
    q: Any,
    active_indices: Any,
    neighbors: Any,
    solver_mode: Any,
    # Boundary Data
    bc_mask: Any,
    bc_data: Any,
    x: Any,
    y: Any,
    # Grid info
    num_active: int,
    rhs: Any,
    # Geometry
    nodes_1d: Any,
    D1D: Any,
    dg_L: Any,
    dg_R: Any,
    root_bounds: Any,
    block_levels: Any,
    # Physics
    params: Any, 
    t: Any
):
    """
    Computes the Flux Reconstruction (FR) update for the Euler equations.
    Includes Taylor Series Extrapolation for Geometric Consistency at boundaries.
    """
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
    
    # Check Solver Mode (Skip if FV)
    if solver_mode[pool_idx] == 1:
        return
    
    # 2D Node Index (i, j)
    j = node_local // N1
    i = node_local % N1
    
    # Grid Spacing (Local AMR)
    my_level = block_levels[pool_idx]
    grid_dim = 1 << my_level
    domain_w = root_bounds[2] - root_bounds[0]
    domain_h = root_bounds[3] - root_bounds[1]
    
    template = q[0, 0][0]
    f_grid_dim = u.get_any_generic(template, grid_dim)
    dx = domain_w / f_grid_dim
    dy = domain_h / f_grid_dim
    
    # Geometric Factors
    one = u.get_one_generic(template)
    two = one + one
    inv_J_x = two / dx
    inv_J_y = two / dy
    zero = one - one
    ramp_time = params.ramp_up_time
    
    # Extrapolation Distances (Reference Space)
    # dist_L = -1.0 - nodes_1d[0]
    # dist_R = 1.0 - nodes_1d[N1-1]
    dist_L = -one - nodes_1d[0]
    dist_R = one - nodes_1d[N1-1]
    
    # --- 1. Volume Gradient (Divergence) + Extrapolation Gradients ---
    val_div = u.make_vec4_generic(zero, zero, zero, zero)
    
    # Accumulators for Face Extrapolation (Reference gradients dQ/dxi)
    dq_dxi_L = u.make_vec4_generic(zero, zero, zero, zero) # For Left Face (i=0)
    dq_dxi_R = u.make_vec4_generic(zero, zero, zero, zero) # For Right Face (i=N1-1)
    
    dq_deta_B = u.make_vec4_generic(zero, zero, zero, zero) # For Bottom Face (j=0)
    dq_deta_T = u.make_vec4_generic(zero, zero, zero, zero) # For Top Face (j=N1-1)

    for k in range(N1):
        # --- X-Direction Integration ---
        # dF/dx terms for the node
        idx_k_x = j * N1 + k
        q_k_x = q[pool_idx, idx_k_x]
        F_k = euler.flux_x(q_k_x, params)
        val_div += F_k * D1D[i, k] * inv_J_x
        
        # Extrapolation: Compute dQ/dxi for the Left (k=0) and Right (k=N1-1) boundaries of THIS row (j)
        # We need dQ/dxi at the boundary node to extrapolate to the face.
        # dQ/dxi @ node 0 = sum(D1D[0, k] * Q[k])
        dq_dxi_L += q_k_x * D1D[0, k]
        dq_dxi_R += q_k_x * D1D[N1-1, k]
        
        # --- Y-Direction Integration ---
        # dG/dy terms for the node
        idx_k_y = k * N1 + i
        q_k_y = q[pool_idx, idx_k_y]
        G_k = euler.flux_y(q_k_y, params)
        val_div += G_k * D1D[j, k] * inv_J_y
        
        # Extrapolation: Compute dQ/deta for Bottom and Top boundaries of THIS col (i)
        dq_deta_B += q_k_y * D1D[0, k]
        dq_deta_T += q_k_y * D1D[N1-1, k]

    # --- 2. Interface Corrections (X) ---
    # Left Interface (i=0)
    idx_L = j * N1 + 0 
    q_L_internal = q[pool_idx, idx_L]
    f_L_internal = euler.flux_x(q_L_internal, params)
    
    # EXTRAPOLATE to Face: q_face = q_node + dq_dxi * dist
    q_L_face = q_L_internal + dq_dxi_L * dist_L
    
    neigh_L = neighbors[pool_idx, FACE_LEFT]
    if neigh_L >= 0:
        idx_neigh = j * N1 + (N1 - 1)
        q_L_ghost = q[neigh_L, idx_neigh] # Neighbor is used as-is (0th order mismatch accepted for neighbor)
    else:
        # Generic Boundary: Use Extrapolated State
        bc_idx = bc_mask[pool_idx, FACE_LEFT]
        x_val = x[pool_idx, idx_L]
        y_val = y[pool_idx, idx_L]
        nx = -one
        ny = zero
        # Apply BC using the correct FACE value, not internal value
        q_L_ghost = bc.apply_boundary_condition(
            bc_idx, bc_data, q_L_face, 
            nx, ny, x_val, y_val, t, ramp_time, params
        )
        
    F_star_L = compute_interface_flux(q_L_ghost, q_L_face, one, zero, params)
    
    corr_x = u.make_vec4_generic(zero, zero, zero, zero)
    if neigh_L != MORTAR_FLAG:
        corr_x = (F_star_L - f_L_internal) * dg_L[i]
    
    # Right Interface (i=N1-1)
    idx_R = j * N1 + (N1 - 1)
    q_R_internal = q[pool_idx, idx_R]
    f_R_internal = euler.flux_x(q_R_internal, params)
    
    # EXTRAPOLATE
    q_R_face = q_R_internal + dq_dxi_R * dist_R

    neigh_R = neighbors[pool_idx, FACE_RIGHT]
    if neigh_R >= 0:
        idx_neigh = j * N1 + 0
        q_R_ghost = q[neigh_R, idx_neigh]
    else:
        # Generic Boundary
        bc_idx = bc_mask[pool_idx, FACE_RIGHT]
        x_val = x[pool_idx, idx_R]
        y_val = y[pool_idx, idx_R]
        nx = one
        ny = zero
        q_R_ghost = bc.apply_boundary_condition(
            bc_idx, bc_data, q_R_face, 
            nx, ny, x_val, y_val, t, ramp_time, params
        )
        
    F_star_R = compute_interface_flux(q_R_face, q_R_ghost, one, zero, params)
    if neigh_R != MORTAR_FLAG:
        corr_x += (F_star_R - f_R_internal) * dg_R[i]
    
    corr_x *= inv_J_x
    
    # --- 3. Interface Corrections (Y) ---
    # Bottom Interface (j=0)
    idx_B = 0 * N1 + i
    q_B_internal = q[pool_idx, idx_B]
    g_B_internal = euler.flux_y(q_B_internal, params)
    
    # EXTRAPOLATE
    q_B_face = q_B_internal + dq_deta_B * dist_L

    neigh_B = neighbors[pool_idx, FACE_BOTTOM]
    if neigh_B >= 0:
        idx_neigh = (N1 - 1) * N1 + i
        q_B_ghost = q[neigh_B, idx_neigh]
    else:
        # Generic Boundary
        bc_idx = bc_mask[pool_idx, FACE_BOTTOM]
        x_val = x[pool_idx, idx_B]
        y_val = y[pool_idx, idx_B]
        nx = zero
        ny = -one
        q_B_ghost = bc.apply_boundary_condition(
            bc_idx, bc_data, q_B_face, 
            nx, ny, x_val, y_val, t, ramp_time, params
        )
        
    G_star_B = compute_interface_flux(q_B_ghost, q_B_face, zero, one, params)
    corr_y = u.make_vec4_generic(zero, zero, zero, zero)
    if neigh_B != MORTAR_FLAG:
        corr_y = (G_star_B - g_B_internal) * dg_L[j]
    
    # Top Interface (j=N1-1)
    idx_T = (N1 - 1) * N1 + i
    q_T_internal = q[pool_idx, idx_T]
    g_T_internal = euler.flux_y(q_T_internal, params)
    
    # EXTRAPOLATE
    q_T_face = q_T_internal + dq_deta_T * dist_R

    neigh_T = neighbors[pool_idx, FACE_TOP]
    if neigh_T >= 0:
        idx_neigh = 0 * N1 + i
        q_T_ghost = q[neigh_T, idx_neigh]
    else:
        # Generic Boundary
        bc_idx = bc_mask[pool_idx, FACE_TOP]
        x_val = x[pool_idx, idx_T]
        y_val = y[pool_idx, idx_T]
        nx = zero
        ny = one
        q_T_ghost = bc.apply_boundary_condition(
            bc_idx, bc_data, q_T_face, 
            nx, ny, x_val, y_val, t, ramp_time, params
        )
        
    G_star_T = compute_interface_flux(q_T_face, q_T_ghost, zero, one, params)
    if neigh_T != MORTAR_FLAG:
        corr_y += (G_star_T - g_T_internal) * dg_R[j]
    
    corr_y *= inv_J_y
    
    # Total Update
    rhs[pool_idx, node_local] += -(val_div + corr_x + corr_y)