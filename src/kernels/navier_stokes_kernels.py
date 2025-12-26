import warp as wp
from src.kernels import boundary_conditions as bc
from src.physics.models.navier_stokes_laws import temperature, viscous_flux_x, viscous_flux_y
from typing import Any

@wp.kernel
def compute_primitive_gradients_volume(
    q: wp.array(dtype=Any, ndim=2),      # (NumElems, Np)
    grad_u: wp.array(dtype=Any, ndim=2),
    grad_v: wp.array(dtype=Any, ndim=2),
    grad_T: wp.array(dtype=Any, ndim=2),
    Dr: wp.array(dtype=Any, ndim=2),
    Ds: wp.array(dtype=Any, ndim=2),
    rx: wp.array(dtype=Any, ndim=2),
    ry: wp.array(dtype=Any, ndim=2),
    sx: wp.array(dtype=Any, ndim=2),
    sy: wp.array(dtype=Any, ndim=2),
    Np: wp.int32,
    params: Any
):
    """
    Pass 1: Computes volume part of gradients of primitive variables.
    grad(phi) = (dphi/dr * dr/dx + dphi/ds * ds/dx, dphi/dr * dr/dy + dphi/ds * ds/dy)
    """
    e, i = wp.tid()
    
    dr_dx = rx[e, i]
    dr_dy = ry[e, i]
    ds_dx = sx[e, i]
    ds_dy = sy[e, i]
    
    du_dr = params.one - params.one
    du_ds = params.one - params.one
    dv_dr = params.one - params.one
    dv_ds = params.one - params.one
    dT_dr = params.one - params.one
    dT_ds = params.one - params.one
    
    for j in range(Np):
        q_j = q[e, j]
        rho_j = wp.max(q_j[0], params.rho_floor)
        u_j = q_j[1] / rho_j
        v_j = q_j[2] / rho_j
        T_j = temperature(q_j, params)
        
        dr = Dr[i, j]
        ds = Ds[i, j]
        
        du_dr += u_j * dr
        du_ds += u_j * ds
        dv_dr += v_j * dr
        dv_ds += v_j * ds
        dT_dr += T_j * dr
        dT_ds += T_j * ds
        
    grad_u[e, i] = bc.make_vec2_generic(du_dr * dr_dx + du_ds * ds_dx, du_dr * dr_dy + du_ds * ds_dy)
    grad_v[e, i] = bc.make_vec2_generic(dv_dr * dr_dx + dv_ds * ds_dx, dv_dr * dr_dy + dv_ds * ds_dy)
    grad_T[e, i] = bc.make_vec2_generic(dT_dr * dr_dx + dT_ds * ds_dx, dT_dr * dr_dy + dT_ds * ds_dy)

@wp.kernel
def compute_primitive_gradients_surface(
    q: wp.array(dtype=Any, ndim=2),
    grad_u: wp.array(dtype=Any, ndim=2),
    grad_v: wp.array(dtype=Any, ndim=2),
    grad_T: wp.array(dtype=Any, ndim=2),
    connectivity: wp.array(dtype=wp.int32, ndim=2),
    neighbor_face_indices: wp.array(dtype=wp.int32, ndim=2),
    face_map: wp.array(dtype=wp.int32, ndim=2),
    LIFT: wp.array(dtype=Any, ndim=2),
    face_geo_factors: wp.array(dtype=Any, ndim=3), 
    J: wp.array(dtype=Any, ndim=2),     
    bc_mask: wp.array(dtype=wp.int32, ndim=2),
    bc_data: wp.array(dtype=Any, ndim=1),
    coord_x: wp.array(dtype=Any, ndim=2),
    coord_y: wp.array(dtype=Any, ndim=2),
    Nfp: wp.int32,
    t: Any,
    ramp_time: Any,
    params: Any
):
    """
    Pass 2: Lifts surface jumps of primitive variables to the gradients (BR1 method).
    grad(phi) += LIFT * (phi_star - phi_inner) * n / J
    where phi_star = 0.5 * (phi_inner + phi_outer)
    """
    e = wp.tid()
    
    for face_idx in range(4):
        nx = face_geo_factors[e, face_idx, 0]
        ny = face_geo_factors[e, face_idx, 1]
        surf_J = face_geo_factors[e, face_idx, 2]
        
        neighbor_e = connectivity[e, face_idx]
        
        for k in range(Nfp):
            node_idx_local = face_map[face_idx, k]
            q_i = q[e, node_idx_local]
            rho_i = wp.max(q_i[0], params.rho_floor)
            u_i = q_i[1] / rho_i
            v_i = q_i[2] / rho_i
            T_i = temperature(q_i, params)
            
            u_o = u_i
            v_o = v_i
            T_o = T_i
            
            if neighbor_e >= 0:
                neighbor_face = neighbor_face_indices[e, face_idx]
                neighbor_node_idx = face_map[neighbor_face, k]
                q_o = q[neighbor_e, neighbor_node_idx]
                rho_o = wp.max(q_o[0], params.rho_floor)
                u_o = q_o[1] / rho_o
                v_o = q_o[2] / rho_o
                T_o = temperature(q_o, params)
            else:
                bc_index = bc_mask[e, face_idx]
                x = coord_x[e, node_idx_local]
                y = coord_y[e, node_idx_local]
                q_bc = bc.apply_boundary_condition(bc_index, bc_data, q_i, nx, ny, x, y, t, ramp_time, params)
                rho_bc = wp.max(q_bc[0], params.rho_floor)
                u_o = q_bc[1] / rho_bc
                v_o = q_bc[2] / rho_bc
                T_o = temperature(q_bc, params)
                
            # Central flux for gradients
            u_star = params.half * (u_i + u_o)
            v_star = params.half * (v_i + v_o)
            T_star = params.half * (T_i + T_o)
            
            # Jump contributions
            jump_u = (u_star - u_i) * surf_J
            jump_v = (v_star - v_i) * surf_J
            jump_T = (T_star - T_i) * surf_J
            
            lift_col = face_idx * Nfp + k
            for i in range(q.shape[1]):
                lift_val = LIFT[i, lift_col]
                inv_vol_J = params.one / J[e, i]
                
                common = lift_val * inv_vol_J
                grad_u[e, i] += bc.make_vec2_generic(jump_u * nx * common, jump_u * ny * common)
                grad_v[e, i] += bc.make_vec2_generic(jump_v * nx * common, jump_v * ny * common)
                grad_T[e, i] += bc.make_vec2_generic(jump_T * nx * common, jump_T * ny * common)

@wp.kernel
def compute_viscous_volume_term(
    q: wp.array(dtype=Any, ndim=2),
    grad_u: wp.array(dtype=Any, ndim=2),
    grad_v: wp.array(dtype=Any, ndim=2),
    grad_T: wp.array(dtype=Any, ndim=2),
    rhs: wp.array(dtype=Any, ndim=2),
    Dr: wp.array(dtype=Any, ndim=2),
    Ds: wp.array(dtype=Any, ndim=2),
    rx: wp.array(dtype=Any, ndim=2),
    ry: wp.array(dtype=Any, ndim=2),
    sx: wp.array(dtype=Any, ndim=2),
    sy: wp.array(dtype=Any, ndim=2),
    Np: wp.int32,
    params: Any
):
    """
    Computes div(Fv) and adds it to RHS.
    """
    e, i = wp.tid()
    
    dr_dx = rx[e, i]
    dr_dy = ry[e, i]
    ds_dx = sx[e, i]
    ds_dy = sy[e, i]
    
    dFv_dr = q[e, i] - q[e, i]
    dFv_ds = q[e, i] - q[e, i]
    dGv_dr = q[e, i] - q[e, i]
    dGv_ds = q[e, i] - q[e, i]
    
    for j in range(Np):
        q_j = q[e, j]
        gu_j = grad_u[e, j]
        gv_j = grad_v[e, j]
        gT_j = grad_T[e, j]
        
        Fv = viscous_flux_x(q_j, gu_j, gv_j, gT_j, params)
        Gv = viscous_flux_y(q_j, gu_j, gv_j, gT_j, params)
        
        dr = Dr[i, j]
        ds = Ds[i, j]
        
        dFv_dr += Fv * dr
        dFv_ds += Fv * ds
        dGv_dr += Gv * dr
        dGv_ds += Gv * ds
        
    div_Fv = (dFv_dr * dr_dx + dFv_ds * ds_dx) + (dGv_dr * dr_dy + dGv_ds * ds_dy)
    rhs[e, i] += div_Fv

@wp.kernel
def compute_viscous_surface_term(
    q: wp.array(dtype=Any, ndim=2),
    grad_u: wp.array(dtype=Any, ndim=2),
    grad_v: wp.array(dtype=Any, ndim=2),
    grad_T: wp.array(dtype=Any, ndim=2),
    rhs: wp.array(dtype=Any, ndim=2),
    connectivity: wp.array(dtype=wp.int32, ndim=2),
    neighbor_face_indices: wp.array(dtype=wp.int32, ndim=2),
    face_map: wp.array(dtype=wp.int32, ndim=2),
    LIFT: wp.array(dtype=Any, ndim=2),
    face_geo_factors: wp.array(dtype=Any, ndim=3), 
    J: wp.array(dtype=Any, ndim=2),     
    bc_mask: wp.array(dtype=wp.int32, ndim=2),
    bc_data: wp.array(dtype=Any, ndim=1),
    coord_x: wp.array(dtype=Any, ndim=2),
    coord_y: wp.array(dtype=Any, ndim=2),
    Nfp: wp.int32,
    t: Any,
    ramp_time: Any,
    params: Any
):
    """
    Computes viscous numerical flux surface integral: LIFT * (Fv_star . n - Fv_inner . n) / J
    """
    e = wp.tid()
    
    for face_idx in range(4):
        nx = face_geo_factors[e, face_idx, 0]
        ny = face_geo_factors[e, face_idx, 1]
        surf_J = face_geo_factors[e, face_idx, 2]
        
        neighbor_e = connectivity[e, face_idx]
        
        for k in range(Nfp):
            node_idx_local = face_map[face_idx, k]
            
            qi = q[e, node_idx_local]
            gui = grad_u[e, node_idx_local]
            gvi = grad_v[e, node_idx_local]
            gTi = grad_T[e, node_idx_local]
            
            Fvi = viscous_flux_x(qi, gui, gvi, gTi, params) * nx + viscous_flux_y(qi, gui, gvi, gTi, params) * ny
            
            qo = qi
            guo = gui
            gvo = gvi
            gTo = gTi
            
            if neighbor_e >= 0:
                neighbor_face = neighbor_face_indices[e, face_idx]
                neighbor_node_idx = face_map[neighbor_face, k]
                qo = q[neighbor_e, neighbor_node_idx]
                guo = grad_u[neighbor_e, neighbor_node_idx]
                gvo = grad_v[neighbor_e, neighbor_node_idx]
                gTo = grad_T[neighbor_e, neighbor_node_idx]
            else:
                bc_index = bc_mask[e, face_idx]
                x = coord_x[e, node_idx_local]
                y = coord_y[e, node_idx_local]
                qo = bc.apply_boundary_condition(bc_index, bc_data, qi, nx, ny, x, y, t, ramp_time, params)
                guo = gui
                gvo = gvi
                gTo = gTi
                
            Fvo = viscous_flux_x(qo, guo, gvo, gTo, params) * nx + viscous_flux_y(qo, guo, gvo, gTo, params) * ny
            
            # Central viscous flux star
            Fv_star = params.half * (Fvi + Fvo)
            
            # RHS contribution (strong form)
            jump = (Fv_star - Fvi) * surf_J
            
            lift_col = face_idx * Nfp + k
            for i in range(q.shape[1]):
                lift_val = LIFT[i, lift_col]
                inv_vol_J = params.one / J[e, i]
                rhs[e, i] += (lift_val * inv_vol_J) * jump