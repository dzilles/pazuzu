import warp as wp
from src.kernels import boundary_conditions as bc
from src.physics.laws.euler import pressure, flux_x, flux_y
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
def rusanov_flux(
    q_l: Any, q_r: Any, 
    F_l_n: Any, F_r_n: Any, # <--- NEW: Pre-calculated Normal Fluxes
    nx: Any, ny: Any, params: Any
):
    """
    Computes the Lax-Friedrichs / Rusanov numerical flux using PRE-CALCULATED fluxes.
    """
    # Wave speeds still depend on state q
    lambda_l = get_max_wave_speed(q_l, nx, ny, params)
    lambda_r = get_max_wave_speed(q_r, nx, ny, params)
    alpha = wp.max(lambda_l, lambda_r)
    
    # Use the Projected Fluxes for the average part
    # F* = 0.5 * (F_L + F_R) - 0.5 * alpha * (q_R - q_L)
    return params.half * (F_l_n + F_r_n) - params.half * alpha * (q_r - q_l)

@wp.func
def hllc_flux(
    q_l: Any, q_r: Any, 
    F_l_n: Any, F_r_n: Any, # <--- NEW: Pre-calculated Normal Fluxes
    nx: Any, ny: Any, params: Any
):
    """
    Computes HLLC flux using PRE-CALCULATED fluxes.
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

    un_l = u_l * nx + v_l * ny
    un_r = u_r * nx + v_r * ny

    # Sound speeds
    c_l = wp.sqrt(params.gamma * p_l / rho_l)
    c_r = wp.sqrt(params.gamma * p_r / rho_r)

    # 2. Wave Speed Estimates (Davis estimate)
    s_l = wp.min(un_l - c_l, un_r - c_r)
    s_r = wp.max(un_l + c_l, un_r + c_r)

    # 3. Logic with Pre-calculated Fluxes
    if s_l >= 0.0:
        return F_l_n
    elif s_r <= 0.0:
        return F_r_n
    else:
        # Subsonic / Contact wave involved
        
        numer = p_r - p_l + rho_l * un_l * (s_l - un_l) - rho_r * un_r * (s_r - un_r)
        denom = rho_l * (s_l - un_l) - rho_r * (s_r - un_r)
        
        # Safety for denom
        if wp.abs(denom) < params.epsilon:
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
            
            # HLLC modification: F*_l = F_l + S_l (U*_l - U_l)
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
            
            return F_r_n + s_r * (U_star_r - q_r)

# --- Solver Kernels ---

@wp.kernel
def compute_nodal_fluxes(
    q: Any,      # (NumElems, Np)
    f_x: Any,    # Output
    f_y: Any,    # Output
    params: Any
):
    """ Computes fluxes at basis nodes. """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    q_val = q[e, i]
    f_x[e, i] = flux_x(q_val, params)
    f_y[e, i] = flux_y(q_val, params)

@wp.kernel
def compute_volume_term(
    f_x: Any,    # (NumElems, Np) - Flux X
    f_y: Any,    # (NumElems, Np) - Flux Y
    rhs: Any,    # Output RHS
    Dr: Any,     # Differentiation Matrix r
    Ds: Any,     # Differentiation Matrix s
    rx: Any,     # Metric dr/dx (NumElems, Np)
    ry: Any,     # Metric dr/dy (NumElems, Np)
    sx: Any,     # Metric ds/dx (NumElems, Np)
    sy: Any,     # Metric ds/dy (NumElems, Np)
    Np: wp.int32,                        # Number of points per element
    params: Any
):
    """
    Computes the divergence of the flux (volume integral) using precomputed nodal fluxes.
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime # Element e, Node i

    # Load Metrics for this element and node
    dr_dx = rx[e, i]
    dr_dy = ry[e, i]
    ds_dx = sx[e, i]
    ds_dy = sy[e, i]

    zero_vec = f_x[e, i] - f_x[e, i]
    dF_dr = zero_vec
    dF_ds = zero_vec
    dG_dr = zero_vec
    dG_ds = zero_vec

    # Matrix-Vector Multiplication for derivatives
    for j in range(Np):
        F_val = f_x[e, j]
        G_val = f_y[e, j]
        
        dr = Dr[i, j]
        ds = Ds[i, j]
        
        dF_dr += F_val * dr
        dF_ds += F_val * ds
        dG_dr += G_val * dr
        dG_ds += G_val * ds

    # Apply Chain Rule
    dF_dx = dF_dr * dr_dx + dF_ds * ds_dx
    dG_dy = dG_dr * dr_dy + dG_ds * ds_dy
    
    rhs[e, i] = -(dF_dx + dG_dy) 

# --- Over-Integration Kernels ---

@wp.kernel
def interpolate_to_quadrature(
    q: Any,      # (NumElems, Np)
    q_q: Any,    # (NumElems, Nq)
    Interp: Any, # (Nq, Np)
    Np: wp.int32
):
    """ Interpolates nodal values to quadrature points. """
    e, iq = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    val = q[e, 0] - q[e, 0]
    for j in range(Np):
        val += q[e, j] * Interp[iq, j]
    q_q[e, iq] = val

@wp.kernel
def compute_projected_fluxes(
    q_q: Any,    # (NumElems, Nq)
    f_x_n: Any,  # (NumElems, Np) - Output Projected Flux X
    f_y_n: Any,  # (NumElems, Np) - Output Projected Flux Y
    Proj: Any,   # (Np, Nq)
    Nq: wp.int32,
    params: Any
):
    """
    Computes non-linear fluxes at quadrature points and projects them back to GLL nodes.
    This effectively performs a L2 projection of the non-linear flux into the polynomial space,
    filtering out high-frequency aliasing modes.
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime # Element e, Node i

    zero_vec = q_q[e, 0] - q_q[e, 0]
    fx_acc = zero_vec
    fy_acc = zero_vec

    for iq in range(Nq):
        q_val = q_q[e, iq]
        
        # Compute non-linear flux at quadrature point
        fx_q = flux_x(q_val, params)
        fy_q = flux_y(q_val, params)
        
        # Accumulate projection: sum_iq (Proj[i, iq] * F_q[iq])
        p_val = Proj[i, iq]
        fx_acc += fx_q * p_val
        fy_acc += fy_q * p_val

    f_x_n[e, i] = fx_acc
    f_y_n[e, i] = fy_acc

@wp.kernel
def accumulate_interface_fluxes(
    q: Any,        
    f_x: Any,      # Projected Flux X
    f_y: Any,      # Projected Flux Y
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
    flux_type: wp.int32,
    params: Any
):
    """
    Computes the surface integral by parallelizing over each face node.
    Exploits the sparsity of the LIFT matrix for GLL nodes.
    """
    e, face_idx, k = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    # 1. Geometry and Connectivity
    nx = face_geo_factors[e, face_idx, 0]
    ny = face_geo_factors[e, face_idx, 1]
    surf_J = face_geo_factors[e, face_idx, 2]
    neighbor_e = connectivity[e, face_idx] 
    
    vol_idx = face_map[face_idx, k]
    
    # --- 2. Load Data (States AND Fluxes) ---
    q_inner = q[e, vol_idx]
    
    # Projected Flux from Interior
    fx_in = f_x[e, vol_idx]
    fy_in = f_y[e, vol_idx]
    f_n_inner = fx_in * nx + fy_in * ny
    
    q_outer = q_inner # Default
    f_n_outer = f_n_inner # Default
    
    if neighbor_e >= 0:
        # Neighbor exists: Load its State AND its Projected Flux
        neighbor_face = neighbor_face_indices[e, face_idx]
        neighbor_node_idx = face_map[neighbor_face, k]
        
        q_outer = q[neighbor_e, neighbor_node_idx]
        
        fx_nb = f_x[neighbor_e, neighbor_node_idx]
        fy_nb = f_y[neighbor_e, neighbor_node_idx]
        f_n_outer = fx_nb * nx + fy_nb * ny
        
    else:
        # Boundary: Compute Boundary State
        bc_index = bc_mask[e, face_idx]
        x = coord_x[e, vol_idx]
        y = coord_y[e, vol_idx]
        q_outer = bc.apply_boundary_condition(bc_index, bc_data, q_inner, nx, ny, x, y, t, ramp_time, params)
        
        # Boundary Flux: We use the analytical flux of the boundary state.
        fx_bc = flux_x(q_outer, params)
        fy_bc = flux_y(q_outer, params)
        f_n_outer = fx_bc * nx + fy_bc * ny
    
    # --- 3. Calculate Riemann Flux using Projected Flux inputs ---
    zero_vec = q_inner - q_inner
    f_star = zero_vec
    if flux_type == FLUX_HLLC:
        f_star = hllc_flux(q_inner, q_outer, f_n_inner, f_n_outer, nx, ny, params)
    else:
        f_star = rusanov_flux(q_inner, q_outer, f_n_inner, f_n_outer, nx, ny, params)
    
    # --- FIX START ---
    # Explicitly enforce No-Penetration Condition on the Numerical Flux.
    # For any solid wall (Slip or No-Slip), u_n = 0.
    # Therefore, Convective Mass Flux (index 0) and Energy Flux (index 3) MUST be zero.
    if neighbor_e < 0:
        bc_index = bc_mask[e, face_idx]
        bc_type = bc_data[bc_index].type
        
        # Check for any Wall type (Slip, No-Slip, Isothermal, Cylinder)
        is_wall = (bc_type == bc.BC_WALL) or \
                  (bc_type == bc.BC_NO_SLIP_WALL) or \
                  (bc_type == bc.BC_ISOTHERMAL_WALL) or \
                  (bc_type == bc.BC_CYLINDER_WALL)
        
        if is_wall:
            zero = params.one - params.one
            # Zero out Mass Flux (rho * u_n)
            f_star = bc.set_vec4_generic(f_star, 0, zero)
            # Zero out Energy Flux ( (rho*E + P) * u_n )
            f_star = bc.set_vec4_generic(f_star, 3, zero)
            # Note: Momentum flux (indices 1, 2) remains as computed (Pressure forces)
    # --- FIX END ---

    # --- 4. Flux Jump ---
    flux_jump = (f_n_inner - f_star) * surf_J
    
    # --- 5. LIFT Operator (Sparse Application) ---
    # For GLL Nodal basis, each face node maps to exactly one volume node.
    lift_col = face_idx * Nfp + k
    lift_val = LIFT[vol_idx, lift_col]
    
    # Accumulate contribution to the volume node
    # Since corner nodes are shared by two faces, we use atomic_add
    vol_J = J[e, vol_idx]
    val = (lift_val * flux_jump) / vol_J
    wp.atomic_add(rhs, e, vol_idx, val)

@wp.kernel
def compute_max_wave_speed(
    q: Any,      # Shape: (num_elements, Np)
    max_speed: Any, # Shape: (1,)
    params: Any
):
    """
    Computes the maximum wave speed in the entire domain for CFL calculation.
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
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

# --- Limiter Kernels ---

@wp.kernel
def compute_cell_averages(
    q: Any,
    q_avg: Any,
    weights: Any,
    J: Any,
    Np: wp.int32,
    params: Any
):
    """
    Computes the cell-average state for each element.
    Avg = (Sum Q_j * w_j * J_j) / (Sum w_j * J_j)
    """
    e = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    # Use template to get correct zero vector and zero scalar of the same precision
    zero_vec = q[e, 0] - q[e, 0]
    zero_scalar = J[e, 0] - J[e, 0]
    
    total_q = zero_vec
    total_vol = zero_scalar
    
    for j in range(Np):
        w_j = weights[j]
        J_j = J[e, j]
        vol_j = w_j * J_j
        
        total_q += q[e, j] * vol_j
        total_vol += vol_j
        
    q_avg[e] = total_q / total_vol

@wp.kernel
def compute_neighbor_min_max(
    q_avg: Any,
    q_min: Any,
    q_max: Any,
    connectivity: Any,
    params: Any
):
    """
    Finds the minimum and maximum cell averages among an element and its neighbors.
    """
    e = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    avg_e = q_avg[e]
    
    min_val = avg_e
    max_val = avg_e
    
    for face_idx in range(4):
        neighbor_e = connectivity[e, face_idx]
        if neighbor_e >= 0:
            avg_nb = q_avg[neighbor_e]
            
            # Warp vector min/max
            # Note: For vec4, we want element-wise min/max
            for c in range(4):
                if avg_nb[c] < min_val[c]:
                    min_val = bc.set_vec4_generic(min_val, c, avg_nb[c])
                if avg_nb[c] > max_val[c]:
                    max_val = bc.set_vec4_generic(max_val, c, avg_nb[c])
                    
    q_min[e] = min_val
    q_max[e] = max_val

@wp.func
def minmod(a: Any, b: Any, c: Any):
    """Standard 3-argument minmod function."""
    zero = a - a
    res = zero
    if a > zero and b > zero and c > zero:
        res = wp.min(a, wp.min(b, c))
    elif a < zero and b < zero and c < zero:
        res = wp.max(a, wp.max(b, c))
    return res

@wp.kernel
def compute_gradients_green_gauss(
    q_avg: Any,
    connectivity: Any,
    face_geo_factors: Any,
    vol: Any,
    grad_x: Any,
    grad_y: Any,
    params: Any
):
    """
    Estimates the cell-center gradient using Green-Gauss theorem.
    grad(q) = (1/Vol) * sum_faces (q_face * n * area)
    """
    e = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    avg_e = q_avg[e]
    v_e = vol[e]
    
    zero_vec = avg_e - avg_e
    gx = zero_vec
    gy = zero_vec
    
    for face_idx in range(4):
        nx = face_geo_factors[e, face_idx, 0]
        ny = face_geo_factors[e, face_idx, 1]
        area = face_geo_factors[e, face_idx, 2]
        
        neighbor_e = connectivity[e, face_idx]
        
        q_nb = avg_e # Default for boundary (simple extrapolation)
        if neighbor_e >= 0:
            q_nb = q_avg[neighbor_e]
        
        # Arithmetic average at face
        q_face = params.half * (avg_e + q_nb)
        
        # Accumulate: q_face * n * area
        gx += q_face * nx * area
        gy += q_face * ny * area
        
    grad_x[e] = gx / v_e
    grad_y[e] = gy / v_e

@wp.kernel
def apply_minmod_limiter(
    q: Any,
    q_avg: Any,
    q_min: Any,
    q_max: Any,
    grad_x: Any,
    grad_y: Any,
    centroid: Any,
    coord_x: Any,
    coord_y: Any,
    Np: wp.int32,
    params: Any
):
    """
    Applies a gradient-based Minmod slope limiter.
    Ensures reconstructed nodal values are within neighbor min/max.
    """
    e = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    avg_e = q_avg[e]
    gx = grad_x[e]
    gy = grad_y[e]
    c_e = centroid[e]
    
    min_e = q_min[e]
    max_e = q_max[e]
    
    # Use template to get correct zero/one scalars of the same precision
    zero_scalar = avg_e[0] - avg_e[0]
    one_scalar = params.one
    
    # Standard Minmod-like slope limiter phi
    phi = one_scalar
    eps = params.rho_floor * params.rho_floor # Small tolerance
    
    for j in range(Np):
        dx = coord_x[e, j] - c_e[0]
        dy = coord_y[e, j] - c_e[1]
        
        # Reconstruction: q_j = q_avg + phi * (grad_q dot delta_x)
        dq = gx * dx + gy * dy
        
        for c in range(4):
            if wp.abs(dq[c]) > eps:
                if dq[c] > zero_scalar:
                    ratio = (max_e[c] - avg_e[c]) / dq[c]
                    phi = wp.min(phi, ratio)
                else:
                    ratio = (min_e[c] - avg_e[c]) / dq[c]
                    phi = wp.min(phi, ratio)
                    
    phi = wp.clamp(phi, zero_scalar, one_scalar)
    
    # Apply limited gradient reconstruction
    if phi < one_scalar:
        for j in range(Np):
            dx = coord_x[e, j] - c_e[0]
            dy = coord_y[e, j] - c_e[1]
            q[e, j] = avg_e + phi * (gx * dx + gy * dy)
    else:
        # Even if phi=1, we still reconstruct to ensure linear consistency 
        # (or we could just leave high-order DG nodes as is, but that's risky for shocks)
        # Actually, for DG, if phi=1 we usually keep the original DG polynomial.
        # But this is a slope limiter which reduces DG to P1-limited.
        for j in range(Np):
            dx = coord_x[e, j] - c_e[0]
            dy = coord_y[e, j] - c_e[1]
            q[e, j] = avg_e + gx * dx + gy * dy

@wp.kernel
def apply_barth_jespersen_limiter(
    q: Any,
    q_avg: Any,
    q_min: Any,
    q_max: Any,
    Np: wp.int32,
    params: Any
):
    """
    Applies the Barth-Jespersen limiter to the nodal values.
    Computes a scaling factor alpha_e such that the limited values are within [q_min, q_max].
    """
    e = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    avg_e = q_avg[e]
    min_e = q_min[e]
    max_e = q_max[e]
    
    # Use template to get correct zero/one scalars of the same precision
    zero_scalar = q_avg[e][0] - q_avg[e][0]
    one_scalar = params.one
    
    # Initialize alpha to 1.0 (unlimited)
    alpha = one_scalar
    eps = params.rho_floor * params.rho_floor
    
    for j in range(Np):
        q_j = q[e, j]
        diff = q_j - avg_e
        
        for c in range(4):
            if diff[c] > eps:
                # Potential overshoot
                ratio = (max_e[c] - avg_e[c]) / diff[c]
                alpha = wp.min(alpha, ratio)
            elif diff[c] < -eps:
                # Potential undershoot
                ratio = (min_e[c] - avg_e[c]) / diff[c]
                alpha = wp.min(alpha, ratio)
                
    # Safeguard alpha
    alpha = wp.clamp(alpha, zero_scalar, one_scalar)
    
    # Apply limiting
    if alpha < one_scalar:
        for j in range(Np):
            q[e, j] = avg_e + alpha * (q[e, j] - avg_e)