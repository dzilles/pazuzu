import warp as wp
from src.kernels import utils as u
from src.numerics import riemann_solvers as rs
from typing import Any
from src.kernels import boundary_conditions as bc # Make sure this is imported

# Face Indices (Same as FR kernels)
FACE_BOTTOM = 0
FACE_RIGHT = 1
FACE_TOP = 2
FACE_LEFT = 3

@wp.func
def get_subcell_geometry(
    weights_1d: Any,
    grid_dim: int,
    root_bounds: Any,
    i: int,
    j: int
):
    """
    Computes the physical area of the sub-cell (i, j).
    Reference area = w_i * w_j.
    Physical area = Reference area * (Lx/2 * Ly/2).
    """
    # 1. Reference Weights
    w_i = weights_1d[i]
    w_j = weights_1d[j]
    
    # 2. Physical scaling
    # Block size
    domain_w = root_bounds[2] - root_bounds[0]
    domain_h = root_bounds[3] - root_bounds[1]
    
    # One block's physical size
    # Cast grid_dim to the same type as domain_w using get_any_generic
    f_grid_dim = u.get_any_generic(domain_w, grid_dim)
    
    block_w = domain_w / f_grid_dim
    block_h = domain_h / f_grid_dim
    
    # Jacobian Det = (block_w/2) * (block_h/2)
    # Use generic helper to ensure 0.5 matches the type of block_w/block_h
    half_w = block_w * u.get_half_generic(block_w)
    half_h = block_h * u.get_half_generic(block_h)
    
    jacobian = half_w * half_h
    
    area = w_i * w_j * jacobian
    return area, block_w, block_h

@wp.kernel
def compute_fv_update(
    q: Any,
    rhs: Any,
    active_indices: Any,
    neighbors: Any,
    solver_mode: Any,
    # --- NEW: Boundary Data ---
    bc_mask: Any,
    bc_data: Any,
    x: Any,
    y: Any,
    # --------------------------
    weights_1d: Any,
    root_bounds: Any,
    block_levels: Any,
    params: Any,
    t: Any # Needed for time-dependent BCs (inlet)
):
    tid = wp.tid()
    N1 = weights_1d.shape[0]
    Np = N1 * N1
    num_active = active_indices.shape[0]
    
    if tid >= num_active * Np:
        return
        
    block_offset = tid // Np
    node_local = tid % Np
    pool_idx = active_indices[block_offset]
    
    if solver_mode[pool_idx] != 1:
        return
        
    j = node_local // N1
    i = node_local % N1
    
    my_level = block_levels[pool_idx]
    grid_dim = 1 << my_level
    vol, block_w, block_h = get_subcell_geometry(weights_1d, grid_dim, root_bounds, i, j)
    inv_vol = u.get_one_generic(vol) / vol
    
    one = u.get_one_generic(vol)
    zero = u.get_any_generic(vol, 0.0)
    ramp_time = params.ramp_up_time
    
    q_curr = q[pool_idx, node_local]
    
    # --- X-Direction ---
    dy_phys = weights_1d[j] * (block_h * u.get_half_generic(block_h))
    
    # 1. Left Interface
    q_L_ghost = q_curr
    if i > 0:
        idx_L = j * N1 + (i - 1)
        q_L_ghost = q[pool_idx, idx_L]
    else:
        neigh_L = neighbors[pool_idx, FACE_LEFT]
        if neigh_L >= 0:
            idx_neigh = j * N1 + (N1 - 1)
            q_L_ghost = q[neigh_L, idx_neigh]
        else:
            # --- FIX: Apply Boundary Condition ---
            bc_idx = bc_mask[pool_idx, FACE_LEFT]
            if bc_idx >= 0:
                nx = -one
                ny = zero
                q_L_ghost = bc.apply_boundary_condition(
                    bc_idx, bc_data, q_curr, nx, ny, 
                    x[pool_idx, node_local], y[pool_idx, node_local], 
                    t, ramp_time, params
                )
    
    F_L = rs.rusanov_flux(q_L_ghost, q_curr, one, zero, params)
    
    # 2. Right Interface
    q_R_ghost = q_curr
    if i < N1 - 1:
        idx_R = j * N1 + (i + 1)
        q_R_ghost = q[pool_idx, idx_R]
    else:
        neigh_R = neighbors[pool_idx, FACE_RIGHT]
        if neigh_R >= 0:
            idx_neigh = j * N1 + 0
            q_R_ghost = q[neigh_R, idx_neigh]
        else:
            # --- FIX: Apply Boundary Condition ---
            bc_idx = bc_mask[pool_idx, FACE_RIGHT]
            if bc_idx >= 0:
                nx = one
                ny = zero
                q_R_ghost = bc.apply_boundary_condition(
                    bc_idx, bc_data, q_curr, nx, ny, 
                    x[pool_idx, node_local], y[pool_idx, node_local], 
                    t, ramp_time, params
                )

    F_R = rs.rusanov_flux(q_curr, q_R_ghost, one, zero, params)
    flux_balance_x = (F_R - F_L) * dy_phys
    
    # --- Y-Direction ---
    dx_phys = weights_1d[i] * (block_w * u.get_half_generic(block_w))
    
    # 3. Bottom Interface
    q_B_ghost = q_curr
    if j > 0:
        idx_B = (j - 1) * N1 + i
        q_B_ghost = q[pool_idx, idx_B]
    else:
        neigh_B = neighbors[pool_idx, FACE_BOTTOM]
        if neigh_B >= 0:
            idx_neigh = (N1 - 1) * N1 + i
            q_B_ghost = q[neigh_B, idx_neigh]
        else:
             # --- FIX: Apply Boundary Condition ---
            bc_idx = bc_mask[pool_idx, FACE_BOTTOM]
            if bc_idx >= 0:
                nx = zero
                ny = -one
                q_B_ghost = bc.apply_boundary_condition(
                    bc_idx, bc_data, q_curr, nx, ny, 
                    x[pool_idx, node_local], y[pool_idx, node_local], 
                    t, ramp_time, params
                )

    G_B = rs.rusanov_flux(q_B_ghost, q_curr, zero, one, params)
    
    # 4. Top Interface
    q_T_ghost = q_curr
    if j < N1 - 1:
        idx_T = (j + 1) * N1 + i
        q_T_ghost = q[pool_idx, idx_T]
    else:
        neigh_T = neighbors[pool_idx, FACE_TOP]
        if neigh_T >= 0:
            idx_neigh = 0 * N1 + i
            q_T_ghost = q[neigh_T, idx_neigh]
        else:
             # --- FIX: Apply Boundary Condition ---
            bc_idx = bc_mask[pool_idx, FACE_TOP]
            if bc_idx >= 0:
                nx = zero
                ny = one
                q_T_ghost = bc.apply_boundary_condition(
                    bc_idx, bc_data, q_curr, nx, ny, 
                    x[pool_idx, node_local], y[pool_idx, node_local], 
                    t, ramp_time, params
                )

    G_T = rs.rusanov_flux(q_curr, q_T_ghost, zero, one, params)
    flux_balance_y = (G_T - G_B) * dx_phys
    
    total_div = (flux_balance_x + flux_balance_y) * inv_vol
    rhs[pool_idx, node_local] += -total_div
