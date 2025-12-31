import warp as wp
from src.kernels import boundary_conditions as bc
from src.physics.laws.navier_stokes import temperature, viscous_flux_x, viscous_flux_y
from typing import Any

@wp.kernel
def compute_primitive_gradients_volume(
    q: Any,      # (NumElems, Np)
    grad_u: Any,
    grad_v: Any,
    grad_T: Any,
    Dr: Any,
    Ds: Any,
    rx: Any,
    ry: Any,
    sx: Any,
    sy: Any,
    Np: wp.int32,
    params: Any
):
    """
    Pass 1: Computes volume part of gradients of primitive variables.
    grad(phi) = (dphi/dr * dr/dx + dphi/ds * ds/dx, dphi/dr * dr/dy + dphi/ds * ds/dy)
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    dr_dx = rx[e, i]  # type: ignore # Warp type inference
    dr_dy = ry[e, i]  # type: ignore # Warp type inference
    ds_dx = sx[e, i]  # type: ignore # Warp type inference
    ds_dy = sy[e, i]  # type: ignore # Warp type inference
    
    du_dr = params.one - params.one
    du_ds = params.one - params.one
    dv_dr = params.one - params.one
    dv_ds = params.one - params.one
    dT_dr = params.one - params.one
    dT_ds = params.one - params.one
    
    for j in range(Np):
        q_j = q[e, j]  # type: ignore # Warp type inference
        rho_j = wp.max(q_j[0], params.rho_floor)
        u_j = q_j[1] / rho_j
        v_j = q_j[2] / rho_j
        T_j = temperature(q_j, params)
        
        dr = Dr[i, j]  # type: ignore # Warp type inference
        ds = Ds[i, j]  # type: ignore # Warp type inference
        
        du_dr += u_j * dr
        du_ds += u_j * ds
        dv_dr += v_j * dr
        dv_ds += v_j * ds
        dT_dr += T_j * dr
        dT_ds += T_j * ds
        
    grad_u[e, i] = bc.make_vec2_generic(du_dr * dr_dx + du_ds * ds_dx, du_dr * dr_dy + du_ds * ds_dy)  # type: ignore # Warp type inference
    grad_v[e, i] = bc.make_vec2_generic(dv_dr * dr_dx + dv_ds * ds_dx, dv_dr * dr_dy + dv_ds * ds_dy)  # type: ignore # Warp type inference
    grad_T[e, i] = bc.make_vec2_generic(dT_dr * dr_dx + dT_ds * ds_dx, dT_dr * dr_dy + dT_ds * ds_dy)  # type: ignore # Warp type inference

@wp.kernel
def compute_primitive_gradients_surface(
    q: Any,
    grad_u: Any,
    grad_v: Any,
    grad_T: Any,
    connectivity: Any,
    neighbor_face_indices: Any,
    face_map: Any,
    LIFT: Any,
    face_geo_factors: Any, 
    J: Any,     
    bc_mask: Any,
    bc_data: Any,
    coord_x: Any,
    coord_y: Any,
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
    e = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    for face_idx in range(4):
        nx = face_geo_factors[e, face_idx, 0]  # type: ignore # Warp type inference
        ny = face_geo_factors[e, face_idx, 1]  # type: ignore # Warp type inference
        surf_J = face_geo_factors[e, face_idx, 2]  # type: ignore # Warp type inference
        
        neighbor_e = connectivity[e, face_idx]  # type: ignore # Warp type inference
        
        for k in range(Nfp):
            node_idx_local = face_map[face_idx, k]  # type: ignore # Warp type inference
            q_i = q[e, node_idx_local]  # type: ignore # Warp type inference
            rho_i = wp.max(q_i[0], params.rho_floor)
            u_i = q_i[1] / rho_i
            v_i = q_i[2] / rho_i
            T_i = temperature(q_i, params)
            
            u_o = u_i
            v_o = v_i
            T_o = T_i
            
            if neighbor_e >= 0:
                neighbor_face = neighbor_face_indices[e, face_idx]  # type: ignore # Warp type inference
                neighbor_node_idx = face_map[neighbor_face, k]  # type: ignore # Warp type inference
                q_o = q[neighbor_e, neighbor_node_idx]  # type: ignore # Warp type inference
                rho_o = wp.max(q_o[0], params.rho_floor)
                u_o = q_o[1] / rho_o
                v_o = q_o[2] / rho_o
                T_o = temperature(q_o, params)
            else:
                bc_index = bc_mask[e, face_idx]  # type: ignore # Warp type inference
                x = coord_x[e, node_idx_local]  # type: ignore # Warp type inference
                y = coord_y[e, node_idx_local]  # type: ignore # Warp type inference
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

            if neighbor_e < 0:
                bc_type = bc_data[bc_index].type
                if bc_type == bc.BC_NO_SLIP_WALL or bc_type == bc.BC_ISOTHERMAL_WALL:
                    # For No-Slip walls, the interface value of velocity is exactly zero.
                    u_star_wall = params.one - params.one
                    v_star_wall = params.one - params.one
                    jump_u = (u_star_wall - u_i) * surf_J
                    jump_v = (v_star_wall - v_i) * surf_J
            
            lift_col = face_idx * Nfp + k
            for i in range(q.shape[1]):
                lift_val = LIFT[i, lift_col]  # type: ignore # Warp type inference
                inv_vol_J = params.one / J[e, i]  # type: ignore # Warp type inference
                
                common = lift_val * inv_vol_J
                grad_u[e, i] += bc.make_vec2_generic(jump_u * nx * common, jump_u * ny * common)  # type: ignore # Warp type inference
                grad_v[e, i] += bc.make_vec2_generic(jump_v * nx * common, jump_v * ny * common)  # type: ignore # Warp type inference
                grad_T[e, i] += bc.make_vec2_generic(jump_T * nx * common, jump_T * ny * common)  # type: ignore # Warp type inference

@wp.kernel
def compute_viscous_volume_term(
    q: Any,
    grad_u: Any,
    grad_v: Any,
    grad_T: Any,
    rhs: Any,
    Dr: Any,
    Ds: Any,
    rx: Any,
    ry: Any,
    sx: Any,
    sy: Any,
    Np: wp.int32,
    params: Any
):
    """
    Computes div(Fv) and adds it to RHS.
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    dr_dx = rx[e, i]  # type: ignore # Warp type inference
    dr_dy = ry[e, i]  # type: ignore # Warp type inference
    ds_dx = sx[e, i]  # type: ignore # Warp type inference
    ds_dy = sy[e, i]  # type: ignore # Warp type inference
    
    dFv_dr = q[e, i] - q[e, i]  # type: ignore # Warp type inference
    dFv_ds = q[e, i] - q[e, i]  # type: ignore # Warp type inference
    dGv_dr = q[e, i] - q[e, i]  # type: ignore # Warp type inference
    dGv_ds = q[e, i] - q[e, i]  # type: ignore # Warp type inference
    
    for j in range(Np):
        q_j = q[e, j]  # type: ignore # Warp type inference
        gu_j = grad_u[e, j]  # type: ignore # Warp type inference
        gv_j = grad_v[e, j]  # type: ignore # Warp type inference
        gT_j = grad_T[e, j]  # type: ignore # Warp type inference
        
        Fv = viscous_flux_x(q_j, gu_j, gv_j, gT_j, params)
        Gv = viscous_flux_y(q_j, gu_j, gv_j, gT_j, params)
        
        dr = Dr[i, j]  # type: ignore # Warp type inference
        ds = Ds[i, j]  # type: ignore # Warp type inference
        
        dFv_dr += Fv * dr
        dFv_ds += Fv * ds
        dGv_dr += Gv * dr
        dGv_ds += Gv * ds
        
    div_Fv = (dFv_dr * dr_dx + dFv_ds * ds_dx) + (dGv_dr * dr_dy + dGv_ds * ds_dy)
    rhs[e, i] += div_Fv  # type: ignore # Warp type inference

@wp.kernel
def compute_viscous_surface_term_kernel(
    q: Any,
    grad_u: Any,
    grad_v: Any,
    grad_T: Any,
    rhs: Any,
    connectivity: Any,
    neighbor_face_indices: Any,
    face_map: Any,
    LIFT: Any,
    face_geo_factors: Any, 
    J: Any,     
    bc_mask: Any,
    bc_data: Any,
    coord_x: Any,
    coord_y: Any,
    Nfp: wp.int32,
    t: Any,
    ramp_time: Any,
    params: Any
):
    e = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    compute_viscous_surface_term(
        e, q, grad_u, grad_v, grad_T, rhs, connectivity, neighbor_face_indices, face_map,
        LIFT, face_geo_factors, J, bc_mask, bc_data, coord_x, coord_y, Nfp, t, ramp_time, params
    )

@wp.func
def compute_viscous_surface_term(
    e: int,
    q: Any,
    grad_u: Any,
    grad_v: Any,
    grad_T: Any,
    rhs: Any,
    connectivity: Any,
    neighbor_face_indices: Any,
    face_map: Any,
    LIFT: Any,
    face_geo_factors: Any, 
    J: Any,     
    bc_mask: Any,
    bc_data: Any,
    coord_x: Any,
    coord_y: Any,
    Nfp: wp.int32,
    t: Any,
    ramp_time: Any,
    params: Any
):
    """
    Computes viscous numerical flux surface integral: LIFT * (Fv_star . n - Fv_inner . n) / J
    """
    
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
                
                # For No-Slip Wall (Adiabatic or Isothermal), we must ensure consistent viscous stresses.
                bc_type = bc_data[bc_index].type
                if bc_type == bc.BC_NO_SLIP_WALL or bc_type == bc.BC_ISOTHERMAL_WALL:
                    # 1. Temperature Gradient
                    if bc_type == bc.BC_NO_SLIP_WALL:
                        # Adiabatic Wall (Neumann=0): Heat Flux must be zero.
                        # Reflection: Flip the NORMAL component of grad T.
                        dot_T = gTi[0] * nx + gTi[1] * ny
                        gTo = bc.make_vec2_generic(gTi[0] - params.one * (dot_T + dot_T) * nx, gTi[1] - params.one * (dot_T + dot_T) * ny)
                    else:
                        # Isothermal Wall: T is fixed, grad T remains as extrapolated from interior.
                        gTo = gTi

                    # 2. No-Slip Wall (Dirichlet=0): Velocity must be zero at the interface.
                    # Reflection: Preserve NORMAL component, Flip TANGENTIAL component of grad U/V.
                    # Formula: g_out = 2*(g_in . n) * n - g_in
                    dot_u = gui[0] * nx + gui[1] * ny
                    guo = bc.make_vec2_generic(params.one * (dot_u + dot_u) * nx - gui[0], params.one * (dot_u + dot_u) * ny - gui[1])
                    
                    dot_v = gvi[0] * nx + gvi[1] * ny
                    gvo = bc.make_vec2_generic(params.one * (dot_v + dot_v) * nx - gvi[0], params.one * (dot_v + dot_v) * ny - gvi[1])
                else:
                    gTo = gTi
                    guo = gui
                    gvo = gvi
                
            Fvo = viscous_flux_x(qo, guo, gvo, gTo, params) * nx + viscous_flux_y(qo, guo, gvo, gTo, params) * ny
            
            # Central viscous flux star
            Fv_star = params.half * (Fvi + Fvo)
            
            # Enforce exact zero energy flux for Adiabatic No-Slip walls.
            # This removes spurious viscous work terms (u.tau) that don't cancel perfectly.
            # We use a check on neighbor_e < 0 to ensure we are on a boundary.
            if neighbor_e < 0:
                 bc_type_check = bc_data[bc_index].type
                 if bc_type_check == bc.BC_NO_SLIP_WALL:
                     # Set Energy flux (index 3) to 0.0
                     zero = params.one - params.one
                     Fv_star = bc.set_vec4_generic(Fv_star, 3, zero)
            
            # RHS contribution (strong form)
            jump = (Fv_star - Fvi) * surf_J
            
            lift_col = face_idx * Nfp + k
            for i in range(q.shape[1]):
                lift_val = LIFT[i, lift_col]
                inv_vol_J = params.one / J[e, i]
                rhs[e, i] += (lift_val * inv_vol_J) * jump