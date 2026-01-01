import warp as wp
from src.physics.laws import euler
from src.numerics import riemann_solvers as rs
from src.kernels import utils as u

from typing import Any

@wp.func
def compute_interface_flux(q_L: Any, q_R: Any, nx: Any, ny: Any, params: Any):
    if params.flux_type == 1:
        return rs.hllc_flux(q_L, q_R, nx, ny, params)
    return rs.rusanov_flux(q_L, q_R, nx, ny, params)

@wp.kernel
def compute_fr_update(
    q: Any,
    active_indices: Any,
    neighbors: Any,
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
    
    # Grid Spacing (Local AMR)
    my_level = block_levels[pool_idx]
    grid_dim = 1 << my_level
    domain_w = root_bounds[2] - root_bounds[0]
    domain_h = root_bounds[3] - root_bounds[1]
    
    template = q[0, 0][0]
    f_grid_dim = u.get_any_generic(template, grid_dim)
    dx = domain_w / f_grid_dim
    dy = domain_h / f_grid_dim
    
    # Geometric Factors: d/dx = (2/dx) * d/dr
    one = u.get_one_generic(template)
    two = one + one
    inv_J_x = two / dx
    inv_J_y = two / dy
    
    zero = one - one
    
    # --- 1. Volume Gradient (Divergence) ---
    val_div = q[0, 0] - q[0, 0] # Generic zero vector of correct precision
    
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
    if neigh_L >= 0:
        # Neighbor's Right boundary (i=N1-1)
        idx_neigh = j * N1 + (N1 - 1)
        q_L_ghost = q[neigh_L, idx_neigh]
    else:
        # --- Boundary Condition (neigh_L == -1) ---
        # TODO: Implement Physical BCs (Inflow/Outflow/Wall).
        # Currently defaults to Transmissive (q_ghost = q_internal).
        # For Periodic setups, neigh_L should never be -1.
        pass
        
    F_star_L = compute_interface_flux(q_L_ghost, q_L_internal, one, zero, params)
    
    # Correction: Skip if mortar (mortar kernel handles it)
    corr_x = q[0, 0] - q[0, 0]
    if neigh_L != -2:
        corr_x = (F_star_L - f_L_internal) * dg_L[i]
    
    # Right Interface (i=N1-1)
    idx_R = j * N1 + (N1 - 1)
    q_R_internal = q[pool_idx, idx_R]
    f_R_internal = euler.flux_x(q_R_internal, params)
    
    neigh_R = neighbors[pool_idx, 1]
    q_R_ghost = q_R_internal
    if neigh_R >= 0:
        # Neighbor's Left boundary (i=0)
        idx_neigh = j * N1 + 0
        q_R_ghost = q[neigh_R, idx_neigh]
    else:
        # --- Boundary Condition (neigh_R == -1) ---
        # TODO: Implement Physical BCs. Defaults to Transmissive.
        pass
        
    F_star_R = compute_interface_flux(q_R_internal, q_R_ghost, one, zero, params)
    if neigh_R != -2:
        corr_x += (F_star_R - f_R_internal) * dg_R[i]
    
    corr_x *= inv_J_x
    
    # --- 3. Interface Corrections (Y) ---
    # Bottom Interface (j=0)
    idx_B = 0 * N1 + i
    q_B_internal = q[pool_idx, idx_B]
    g_B_internal = euler.flux_y(q_B_internal, params)
    
    neigh_B = neighbors[pool_idx, 2]
    q_B_ghost = q_B_internal
    if neigh_B >= 0:
        # Neighbor's Top boundary (j=N1-1)
        idx_neigh = (N1 - 1) * N1 + i
        q_B_ghost = q[neigh_B, idx_neigh]
    else:
        # --- Boundary Condition (neigh_B == -1) ---
        # TODO: Implement Physical BCs. Defaults to Transmissive.
        pass
        
    # Flux in Y direction, Normal=(0,1)
    G_star_B = compute_interface_flux(q_B_ghost, q_B_internal, zero, one, params)
    corr_y = q[0, 0] - q[0, 0]
    if neigh_B != -2:
        corr_y = (G_star_B - g_B_internal) * dg_L[j]
    
    # Top Interface (j=N1-1)
    idx_T = (N1 - 1) * N1 + i
    q_T_internal = q[pool_idx, idx_T]
    g_T_internal = euler.flux_y(q_T_internal, params)
    
    neigh_T = neighbors[pool_idx, 3]
    q_T_ghost = q_T_internal
    if neigh_T >= 0:
        # Neighbor's Bottom boundary (j=0)
        idx_neigh = 0 * N1 + i
        q_T_ghost = q[neigh_T, idx_neigh]
    else:
        # --- Boundary Condition (neigh_T == -1) ---
        # TODO: Implement Physical BCs. Defaults to Transmissive.
        pass
        
    G_star_T = compute_interface_flux(q_T_internal, q_T_ghost, zero, one, params)
    if neigh_T != -2:
        corr_y += (G_star_T - g_T_internal) * dg_R[j]
    
    corr_y *= inv_J_y
    
    # Total Update
    rhs[pool_idx, node_local] = -(val_div + corr_x + corr_y)
