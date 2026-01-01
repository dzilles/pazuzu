import warp as wp
from typing import Any
from src.numerics.riemann_solvers import rusanov_flux
from src.physics.laws import euler
from src.kernels import utils as u
from src.kernels.structs import EquationParams32

@wp.kernel
def compute_mortar_fluxes(
    q: Any,                 # (MAX, Np) vec4
    rhs: Any,               # (MAX, Np) vec4 - Flux accumulator
    mortar_list: Any,       # (MAX_MORTARS, 4) int32: [fine_idx, face, coarse_idx, subface]
    num_mortars: Any,       # (1) int32
    # Geometry
    nodes_1d: Any,
    face_nodes: Any,        # (4, Nfp) int32
    dg_L: Any,              # (N1) float
    dg_R: Any,              # (N1) float
    P_left: Any,            # (N1, N1) float
    P_right: Any,           # (N1, N1) float
    R_left: Any,            # (N1, N1) float
    R_right: Any,           # (N1, N1) float
    root_bounds: Any,
    block_levels: Any,
    # Physics
    params: Any
):
    tid = wp.tid()
    count = num_mortars[0]
    if tid >= count:
        return
        
    # 1. Load Mortar Info
    fine_idx = mortar_list[tid, 0]
    face = mortar_list[tid, 1]
    coarse_idx = mortar_list[tid, 2]
    subface = mortar_list[tid, 3] # 0 or 1
    
    # 2. Geometric Setup
    N1 = nodes_1d.shape[0]
    Np = N1 * N1
    domain_w = root_bounds[2] - root_bounds[0]
    domain_h = root_bounds[3] - root_bounds[1]
    
    template = q[0, 0][0]
    zero_s = template - template
    zero_v = u.make_vec4_generic(zero_s, zero_s, zero_s, zero_s)
    one_s = u.get_one_generic(zero_s)
    half_s = u.get_half_generic(zero_s)
    two_s = one_s + one_s
    
    # Normal and Opposing Face
    # Opposing face in B-R-T-L ordering: 0<->2 (B<->T), 1<->3 (R<->L)
    coarse_face = (face + 2) % 4
        
    # Jacobians
    level_f = block_levels[fine_idx]
    dim_f = 1 << level_f
    dx_f = domain_w / u.get_any_generic(zero_s, dim_f)
    dy_f = domain_h / u.get_any_generic(zero_s, dim_f)
    
    # Selection of coordinate spacing based on interface normal
    # Face 0,2 are Y-interfaces (Normal=(0,1)), Face 1,3 are X-interfaces (Normal=(1,0))
    inv_J_f = two_s / dy_f
    if face == 1 or face == 3: inv_J_f = two_s / dx_f
    
    level_c = block_levels[coarse_idx]
    dim_c = 1 << level_c
    dx_c = domain_w / u.get_any_generic(zero_s, dim_c)
    dy_c = domain_h / u.get_any_generic(zero_s, dim_c)
    inv_J_c = two_s / dy_c
    if face == 1 or face == 3: inv_J_c = two_s / dx_c

    # 3. Compute Interface Fluxes (F*)
    # All fluxes along POSITIVE x or y
    for k in range(N1):
        idx_f = face_nodes[face, k]
        q_f = q[fine_idx, idx_f]
        
        # Internal flux component along POSITIVE axis
        f_int_f = zero_v
        if face == 1 or face == 3: f_int_f = euler.flux_x(q_f, params)
        else:                      f_int_f = euler.flux_y(q_f, params)
        
        # Projected Coarse state
        q_c_proj = zero_v
        for m in range(N1):
            idx_c = face_nodes[coarse_face, m]
            weight = P_left[k, m]
            if subface == 1: weight = P_right[k, m]
            q_c_proj += q[coarse_idx, idx_c] * weight
            
        # Numerical Flux (F*) along POSITIVE axis
        # Use (1,0) or (0,1)
        F_star = zero_v
        if face == 1 or face == 3: 
            if face == 3: # Fine is Right, Coarse is Left
                F_star = rusanov_flux(q_c_proj, q_f, one_s, zero_s, params)
            else: # Fine is Left, Coarse is Right
                F_star = rusanov_flux(q_f, q_c_proj, one_s, zero_s, params)
        else:
            if face == 0: # Fine is Top, Coarse is Bottom
                F_star = rusanov_flux(q_c_proj, q_f, zero_s, one_s, params)
            else: # Fine is Bottom, Coarse is Top
                F_star = rusanov_flux(q_f, q_c_proj, zero_s, one_s, params)
        
        # 4. Apply Correction to Fine Block (Volume)
        # jump = F*_pos - f_int_pos
        jump_f = F_star - f_int_f
        
        dg_val_f = dg_L # Default for faces 0 and 3 (-y and -x)
        if (face == 1) or (face == 2): dg_val_f = dg_R # Faces 1 and 2 (+x and +y)
        
        for ii in range(N1):
            xi_idx = ii
            # For Face 1,3 (X), correction varies in i. For Face 0,2 (Y), correction varies in j.
            if face == 1 or face == 3: node_idx = k * N1 + ii
            else:                      node_idx = ii * N1 + k
            
            w_corr = dg_val_f[xi_idx]
            wp.atomic_add(rhs, fine_idx, node_idx, -jump_f * w_corr * inv_J_f)
            
        # 5. Apply Correction to Coarse Block (Volume)
        dg_val_c = dg_L
        if (coarse_face == 1) or (coarse_face == 2): dg_val_c = dg_R
        
        for m in range(N1):
            weight_r = R_left[m, k]
            if subface == 1: weight_r = R_right[m, k]
            
            f_star_c_contrib = F_star * weight_r
            
            for ii in range(N1):
                xi_idx_c = ii
                if coarse_face == 1 or coarse_face == 3: node_idx_c = m * N1 + ii
                else:                                    node_idx_c = ii * N1 + m
                
                w_corr_c = dg_val_c[xi_idx_c]
                # RHS_c -= F_star part
                wp.atomic_add(rhs, coarse_idx, node_idx_c, -f_star_c_contrib * w_corr_c * inv_J_c)
                
    # 6. Finalize Coarse Internal Correction
    dg_val_c = dg_L
    if (coarse_face == 1) or (coarse_face == 2): dg_val_c = dg_R

    for m in range(N1):
        idx_c_m = face_nodes[coarse_face, m]
        q_c_m = q[coarse_idx, idx_c_m]
        
        f_int_c_m = zero_v
        if face == 1 or face == 3: f_int_c_m = euler.flux_x(q_c_m, params)
        else:                      f_int_c_m = euler.flux_y(q_c_m, params)
        
        for ii in range(N1):
            xi_idx_c = ii
            if coarse_face == 1 or coarse_face == 3: node_idx_c = m * N1 + ii
            else:                                    node_idx_c = ii * N1 + m
            
            w_corr_c = dg_val_c[xi_idx_c]
            # RHS_c += f_int_c part (weight 0.5 because coarse covers 2 fine blocks)
            wp.atomic_add(rhs, coarse_idx, node_idx_c, f_int_c_m * half_s * w_corr_c * inv_J_c)