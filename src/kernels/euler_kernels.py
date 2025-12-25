import warp as wp
from src.kernels import boundary_conditions as bc
from src.physics.models.euler_laws import pressure, flux_x, flux_y
from typing import Any

# Use constants from BC module
BC_WALL = bc.BC_WALL
BC_FARFIELD = bc.BC_FARFIELD
BC_INLET = bc.BC_INLET
BC_OUTLET = bc.BC_OUTLET

# Flux Types
FLUX_RUSANOV = 0
FLUX_HLLC = 1

# --- Equation Kernels (Helper functions used by other kernels) ---

@wp.func
def get_max_wave_speed(q: Any, nx: Any, ny: Any, params: Any):
    """
    Calculates the acoustic wave speed |u_n| + c in the direction normal to a face.
    """
    zero = q[0] - q[0]
    rho = wp.max(q[0], zero + params.rho_floor)
    p = pressure(q, params)
    c = wp.sqrt(params.gamma * p / rho)
    
    # Velocity normal to the face
    u_n = (q[1] * nx + q[2] * ny) / rho
    return wp.abs(u_n) + c

@wp.func
def rusanov_flux(q_l: Any, q_r: Any, nx: Any, ny: Any, params: Any):
    """
    Computes the Lax-Friedrichs / Rusanov numerical flux across an interface.
    """
    # Fluxes projected onto normal
    F_l = flux_x(q_l, params) * nx + flux_y(q_l, params) * ny
    F_r = flux_x(q_r, params) * nx + flux_y(q_r, params) * ny
    
    # Wave speeds
    lambda_l = get_max_wave_speed(q_l, nx, ny, params)
    lambda_r = get_max_wave_speed(q_r, nx, ny, params)
    alpha = wp.max(lambda_l, lambda_r)
    
    return params.half * (F_l + F_r) - params.half * alpha * (q_r - q_l)

@wp.func
def hllc_flux(q_l: Any, q_r: Any, nx: Any, ny: Any, params: Any):
    """
    Computes the HLLC (Harten-Lax-van Leer-Contact) numerical flux.
    """
    zero = q_l[0] - q_l[0]
    
    # 1. Primitives and Normal Velocities
    rho_l = wp.max(q_l[0], zero + params.rho_floor)
    u_l = q_l[1] / rho_l
    v_l = q_l[2] / rho_l
    p_l = pressure(q_l, params)
    
    rho_r = wp.max(q_r[0], zero + params.rho_floor)
    u_r = q_r[1] / rho_r
    v_r = q_r[2] / rho_r
    p_r = pressure(q_r, params)

    # Normal velocities
    un_l = u_l * nx + v_l * ny
    un_r = u_r * nx + v_r * ny

    # Sound speeds
    c_l = wp.sqrt(params.gamma * p_l / rho_l)
    c_r = wp.sqrt(params.gamma * p_r / rho_r)

    # 2. Wave Speed Estimates (Davis estimate)
    s_l = wp.min(un_l - c_l, un_r - c_r)
    s_r = wp.max(un_l + c_l, un_r + c_r)

    # 3. Check for Inter-wave states
    if s_l >= 0.0:
        # Supersonic flow to the right -> Flux is F_L
        return flux_x(q_l, params) * nx + flux_y(q_l, params) * ny
    elif s_r <= 0.0:
        # Supersonic flow to the left -> Flux is F_R
        return flux_x(q_r, params) * nx + flux_y(q_r, params) * ny
    else:
        # Subsonic / Contact wave involved
        
        numer = p_r - p_l + rho_l * un_l * (s_l - un_l) - rho_r * un_r * (s_r - un_r)
        denom = rho_l * (s_l - un_l) - rho_r * (s_r - un_r)
        
        # Safety for denom
        if wp.abs(denom) < 1.0e-8:
            s_star = params.half * (un_l + un_r) # Fallback to average
        else:
            s_star = numer / denom

        # Select Side (Left or Right of Contact Discontinuity)
        if s_star >= 0.0:
            # We are in the Star-Left region
            D_l = (s_l - un_l) / (s_l - s_star)
            
            rho_star_l = rho_l * D_l
            u_star_l = u_l + (s_star - un_l) * nx
            v_star_l = v_l + (s_star - un_l) * ny
            
            E_term = (q_l[3] / rho_l) + (s_star - un_l) * (s_star + p_l / (rho_l * (s_l - un_l)))
            E_star_tot_l = rho_star_l * E_term
            
            U_star_l = bc.make_vec4_generic(rho_star_l, rho_star_l * u_star_l, rho_star_l * v_star_l, E_star_tot_l)
            
            # Flux F_l projected
            F_l_n = flux_x(q_l, params) * nx + flux_y(q_l, params) * ny
            
            return F_l_n + s_l * (U_star_l - q_l)

        else:
            # We are in the Star-Right region
            D_r = (s_r - un_r) / (s_r - s_star)
            
            rho_star_r = rho_r * D_r
            u_star_r = u_r + (s_star - un_r) * nx
            v_star_r = v_r + (s_star - un_r) * ny
            
            E_term = (q_r[3] / rho_r) + (s_star - un_r) * (s_star + p_r / (rho_r * (s_r - un_r)))
            E_star_tot_r = rho_star_r * E_term
            
            U_star_r = bc.make_vec4_generic(rho_star_r, rho_star_r * u_star_r, rho_star_r * v_star_r, E_star_tot_r)
            
            F_r_n = flux_x(q_r, params) * nx + flux_y(q_r, params) * ny
            
            return F_r_n + s_r * (U_star_r - q_r)

# --- Solver Kernels ---

@wp.kernel
def compute_volume_term(
    q: wp.array(dtype=Any, ndim=2),      # (NumElems, Np)
    rhs: wp.array(dtype=Any, ndim=2),    # Output RHS
    Dr: wp.array(dtype=Any, ndim=2),  # Differentiation Matrix r
    Ds: wp.array(dtype=Any, ndim=2),  # Differentiation Matrix s
    rx: wp.array(dtype=Any, ndim=2),  # Metric dr/dx (NumElems, Np)
    ry: wp.array(dtype=Any, ndim=2),  # Metric dr/dy (NumElems, Np)
    sx: wp.array(dtype=Any, ndim=2),  # Metric ds/dx (NumElems, Np)
    sy: wp.array(dtype=Any, ndim=2),  # Metric ds/dy (NumElems, Np)
    Np: wp.int32,                             # Number of points per element
    params: Any
):
    """
    Computes the divergence of the flux (volume integral).
    """
    e, i = wp.tid() # Element e, Node i

    # Load Metrics for this element and node
    dr_dx = rx[e, i]
    dr_dy = ry[e, i]
    ds_dx = sx[e, i]
    ds_dy = sy[e, i]

    # Initialize accumulation vectors with zeros of correct type
    zero_vec = q[e, i] - q[e, i]
    
    dF_dr = zero_vec
    dF_ds = zero_vec
    dG_dr = zero_vec
    dG_ds = zero_vec

    # Matrix-Vector Multiplication
    for j in range(Np):
        q_val = q[e, j]
        
        # Fluxes at node j
        F_val = flux_x(q_val, params)
        G_val = flux_y(q_val, params)
        
        # Accumulate derivatives
        dr = Dr[i, j]
        ds = Ds[i, j]
        
        dF_dr += F_val * dr
        dF_ds += F_val * ds
        
        dG_dr += G_val * dr
        dG_ds += G_val * ds

    # Apply Chain Rule
    dF_dx = dF_dr * dr_dx + dF_ds * ds_dx
    dG_dy = dG_dr * dr_dy + dG_ds * ds_dy
    
    # RHS update
    rhs[e, i] = -(dF_dx + dG_dy) 


@wp.kernel
def compute_surface_term(
    q: wp.array(dtype=Any, ndim=2),        # State
    rhs: wp.array(dtype=Any, ndim=2),      # RHS to accumulate into
    connectivity: wp.array(dtype=wp.int32, ndim=2),
    neighbor_face_indices: wp.array(dtype=wp.int32, ndim=2),
    face_map: wp.array(dtype=wp.int32, ndim=2),
    LIFT: wp.array(dtype=Any, ndim=2),
    face_geo_factors: wp.array(dtype=Any, ndim=3), 
    J: wp.array(dtype=Any, ndim=2),     
    bc_mask: wp.array(dtype=wp.int32, ndim=2),
    coord_x: wp.array(dtype=Any, ndim=2),
    coord_y: wp.array(dtype=Any, ndim=2),
    Nfp: wp.int32,
    t: Any,
    ramp_time: Any,
    flux_type: wp.int32,
    params: Any
):
    """
    Computes the surface integral (flux jump) and lifts it to the volume nodes.
    """
    e = wp.tid() # One thread per element

    zero_vec = q[e, 0] - q[e, 0]

    # Loop over all 4 faces of the quadrilateral
    for face_idx in range(4):
        
        # Load geometric factors for this face
        nx = face_geo_factors[e, face_idx, 0]
        ny = face_geo_factors[e, face_idx, 1]
        surf_J = face_geo_factors[e, face_idx, 2]

        # Get Neighbor Info
        neighbor_e = connectivity[e, face_idx] 
        
        # Loop over nodes on this face
        for k in range(Nfp):
            # Node Index on current element
            node_idx_local = face_map[face_idx, k]
            
            q_inner = q[e, node_idx_local]
            q_outer = q_inner # Default initialization
            
            if neighbor_e >= 0:
                # --- Interior Face ---
                neighbor_face = neighbor_face_indices[e, face_idx]
                neighbor_k = k 
                neighbor_node_idx = face_map[neighbor_face, neighbor_k]
                
                q_outer = q[neighbor_e, neighbor_node_idx]
            else:
                # --- Boundary Face ---
                bc_type = bc_mask[e, face_idx]
                x = coord_x[e, node_idx_local]
                y = coord_y[e, node_idx_local]
                q_outer = bc.apply_boundary_condition(bc_type, q_inner, nx, ny, x, y, t, ramp_time, params)
            
            # 1. Numerical Flux (F*)
            f_star = zero_vec
            if flux_type == FLUX_HLLC:
                f_star = hllc_flux(q_inner, q_outer, nx, ny, params)
            else:
                f_star = rusanov_flux(q_inner, q_outer, nx, ny, params)
            
            # 2. Normal Flux from interior (F_n)
            f_n = flux_x(q_inner, params) * nx + flux_y(q_inner, params) * ny
            
            # 3. Flux Jump (F* - F_n) scaled by Surface Jacobian
            flux_jump = (f_n - f_star) * surf_J
            
            # 4. LIFTing: Add contribution to ALL volume nodes
            lift_col = face_idx * Nfp + k
            
            for i in range(q.shape[1]): # Iterate over all volume nodes (Np)
                lift_val = LIFT[i, lift_col]
                vol_J = J[e, i]
                val = lift_val * flux_jump / vol_J
                rhs[e, i] += val

@wp.kernel
def compute_max_wave_speed(
    q: wp.array(dtype=Any, ndim=2),      # Shape: (num_elements, Np)
    max_speed: wp.array(dtype=Any, ndim=1), # Shape: (1,)
    params: Any
):
    """
    Computes the maximum wave speed in the entire domain for CFL calculation.
    """
    e, i = wp.tid()
    
    # Load state
    val = q[e, i]
    zero = val[0] - val[0]
    rho = wp.max(val[0], zero + params.rho_floor)
    
    # Primitive variables
    u = val[1] / rho
    v = val[2] / rho
    p = pressure(val, params)
    
    # Sound speed c
    c = wp.sqrt(params.gamma * p / rho)
    
    # Velocity magnitude |u|
    vel_mag = wp.sqrt(u*u + v*v)
    
    # Wave speed lambda = |u| + c
    wave_speed = vel_mag + c
    
    # Update global maximum atomically
    wp.atomic_max(max_speed, 0, wave_speed)
