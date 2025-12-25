import warp as wp
from typing import Any

# --- Boundary Condition Type Constants ---
BC_INTERNAL = 0
BC_WALL = 1      # Slip Wall (Euler)
BC_FARFIELD = 2  # Farfield (Freestream / Characteristic)
BC_INLET = 3     # Inlet (Freestream)
BC_OUTLET = 4    # Subsonic Outlet
BC_EXTRAPOLATION = 5 # Zero-Gradient / Outflow
BC_CYLINDER_WALL = 6 # Slip Wall with analytical normals for cylinder

# --- Helper Functions ---

@wp.func
def make_vec4_generic(x: wp.float32, y: wp.float32, z: wp.float32, w: wp.float32):
    return wp.vec4(x, y, z, w)

@wp.func
def make_vec4_generic(x: wp.float64, y: wp.float64, z: wp.float64, w: wp.float64):
    return wp.vec4d(x, y, z, w)

@wp.func
def set_vec4_generic(v: wp.vec4, i: int, val: wp.float32):
    res = v
    if i == 0: res = wp.vec4(val, v[1], v[2], v[3])
    elif i == 1: res = wp.vec4(v[0], val, v[2], v[3])
    elif i == 2: res = wp.vec4(v[0], v[1], val, v[3])
    elif i == 3: res = wp.vec4(v[0], v[1], v[2], val)
    return res

@wp.func
def set_vec4_generic(v: wp.vec4d, i: int, val: wp.float64):
    res = v
    if i == 0: res = wp.vec4d(val, v[1], v[2], v[3])
    elif i == 1: res = wp.vec4d(v[0], val, v[2], v[3])
    elif i == 2: res = wp.vec4d(v[0], v[1], val, v[3])
    elif i == 3: res = wp.vec4d(v[0], v[1], v[2], val)
    return res

@wp.func
def get_freestream_state(t: Any, ramp_time: Any, template_q: Any, params: Any):
    """Returns the freestream state with ramping (rho=1, u=ramped, v=0, p=1)."""
    zero = template_q[0] - template_q[0]
    
    rho = zero + params.one
    
    # Ramp u from 0 to 1 over ramp_time
    target_u = zero + params.one
    factor = zero + params.one
    
    if t < ramp_time:
        factor = t / ramp_time
        
    u = target_u * factor
    v = zero
    p = zero + params.one
    
    rho_u = rho * u
    rho_v = rho * v
    kinetic_energy = params.half * rho * (u*u + v*v)
    E = p / (params.gamma - params.one) + kinetic_energy
    
    return make_vec4_generic(rho, rho_u, rho_v, E)

@wp.func
def compute_characteristic_state(q_inner: Any, q_target: Any, nx: Any, ny: Any, params: Any):
    """
    Computes the boundary state using Riemann invariants (Characteristic Boundary Conditions).
    """
    # --- Inner Primitives ---
    rho_i = wp.max(params.rho_floor, q_inner[0])
    inv_rho_i = params.one / rho_i
    u_i = q_inner[1] * inv_rho_i
    v_i = q_inner[2] * inv_rho_i
    
    # Pressure and Sound Speed
    kin_i = params.half * (q_inner[1]*u_i + q_inner[2]*v_i)
    p_i = wp.max(params.p_floor, (params.gamma - params.one) * (q_inner[3] - kin_i))
    c_i = wp.sqrt(params.gamma * p_i * inv_rho_i)
    
    un_i = u_i * nx + v_i * ny
    ut_i = -u_i * ny + v_i * nx # Tangential velocity
    
    # --- Target Primitives ---
    rho_t = wp.max(params.rho_floor, q_target[0])
    inv_rho_t = params.one / rho_t
    u_t = q_target[1] * inv_rho_t
    v_t = q_target[2] * inv_rho_t
    
    # Pressure and Sound Speed
    kin_t = params.half * (q_target[1]*u_t + q_target[2]*v_t)
    p_t = wp.max(params.p_floor, (params.gamma - params.one) * (q_target[3] - kin_t))
    c_t = wp.sqrt(params.gamma * p_t * inv_rho_t)
    
    un_t = u_t * nx + v_t * ny
    ut_t = -u_t * ny + v_t * nx
    
    # --- Riemann Invariants ---
    gm1 = params.gamma - params.one
    inv_gm1 = params.one / gm1
    two = params.one + params.one
    
    # J+ (outgoing from inner), J- (incoming from target)
    j_inner = un_i + two * c_i * inv_gm1
    j_target = un_t - two * c_t * inv_gm1
    
    # --- Flow Regime ---
    zero = params.one - params.one
    
    if un_i > c_i: # Supersonic Outflow
        return q_inner
    elif un_i < -c_i: # Supersonic Inflow
        return q_target
    else:
        # Subsonic Outflow or Inflow
        un_b = params.half * (j_inner + j_target)
        c_b = params.half * params.half * gm1 * (j_inner - j_target)
        
        ut_b = zero
        s_b = zero
        
        if un_i >= zero: # Subsonic Outflow (one characteristic enters, three leave)
            # Entropy and tangential velocity from inner
            s_i = p_i / wp.pow(rho_i, params.gamma)
            ut_b = ut_i
            s_b = s_i
        else: # Subsonic Inflow (three characteristics enter, one leaves)
            # Entropy and tangential velocity from target
            s_t = p_t / wp.pow(rho_t, params.gamma)
            ut_b = ut_t
            s_b = s_t
            
        # Reconstruct boundary state Qb
        # rho = (c^2 / (gamma * s))^(1/(gamma-1))
        rho_b = wp.pow((c_b * c_b) / (params.gamma * s_b), inv_gm1)
        p_b = s_b * wp.pow(rho_b, params.gamma)
        
        # Velocity components from un_b and ut_b
        u_b = un_b * nx - ut_b * ny
        v_b = un_b * ny + ut_b * nx
        
        # Conservative variables
        rho_u_b = rho_b * u_b
        rho_v_b = rho_b * v_b
        kin_b = params.half * rho_b * (u_b*u_b + v_b*v_b)
        E_b = p_b * inv_gm1 + kin_b
        
        return make_vec4_generic(rho_b, rho_u_b, rho_v_b, E_b)

# --- Boundary Condition Functions ---

@wp.func
def apply_slip_wall(q_inner: Any, nx: Any, ny: Any, params: Any):
    """
    Applies Slip Wall boundary condition.
    """
    rho = q_inner[0]
    rhou = q_inner[1]
    rhov = q_inner[2]
    E = q_inner[3]
    
    # Momentum dot Normal
    mom_dot_n = rhou * nx + rhov * ny
    
    # Reflected momentum: v_ghost = v - 2 * (v . n) * n
    # 2.0 = one + one
    two = params.one + params.one
    rhou_ghost = rhou - two * mom_dot_n * nx
    rhov_ghost = rhov - two * mom_dot_n * ny
    
    return make_vec4_generic(rho, rhou_ghost, rhov_ghost, E)

@wp.func
def apply_cylinder_wall(q_inner: Any, x: Any, y: Any, params: Any):
    """
    Applies Slip Wall boundary condition using analytical normals for a cylinder centered at (0,0).
    """
    # Calculate analytical normal for cylinder centered at (0,0)
    radius = wp.sqrt(x*x + y*y)
    
    # Avoid division by zero
    nx_analytisch = x - x + params.one
    ny_analytisch = x - x
    if radius > 1e-8:
        nx_analytisch = x / radius
        ny_analytisch = y / radius
        
    return apply_slip_wall(q_inner, nx_analytisch, ny_analytisch, params)

@wp.func
def apply_extrapolation(q_inner: Any):
    return q_inner

@wp.func
def apply_farfield(q_inner: Any, nx: Any, ny: Any, t: Any, ramp_time: Any, params: Any):
    """
    Applies Characteristic Farfield boundary condition.
    """
    if q_inner[0] < params.rho_floor:
        return get_freestream_state(t, ramp_time, q_inner, params)
        
    q_freestream = get_freestream_state(t, ramp_time, q_inner, params)
    return compute_characteristic_state(q_inner, q_freestream, nx, ny, params)

@wp.func
def apply_inlet(q_inner: Any, nx: Any, ny: Any, t: Any, ramp_time: Any, params: Any):
    """
    Applies Inlet boundary condition using characteristics.
    """
    q_inlet = get_freestream_state(t, ramp_time, q_inner, params)
    return compute_characteristic_state(q_inner, q_inlet, nx, ny, params)

@wp.func
def apply_outlet(q_inner: Any, nx: Any, ny: Any, params: Any):
    """
    Applies Subsonic Outlet boundary condition using characteristics.
    """
    # Target state: rho, u, v from inner, but p = p_back (usually 1.0)
    rho = wp.max(params.rho_floor, q_inner[0])
    inv_rho = params.one / rho
    u = q_inner[1] * inv_rho
    v = q_inner[2] * inv_rho
    
    p_back = params.one
    
    kin = params.half * rho * (u*u + v*v)
    E_target = p_back / (params.gamma - params.one) + kin
    
    q_target = make_vec4_generic(rho, q_inner[1], q_inner[2], E_target)
    
    return compute_characteristic_state(q_inner, q_target, nx, ny, params)

@wp.func
def apply_boundary_condition(
    bc_type: wp.int32, 
    q_inner: Any, 
    nx: Any, 
    ny: Any,
    x: Any,
    y: Any,
    t: Any,
    ramp_time: Any,
    params: Any
):
    """
    Dispatcher for boundary conditions.
    """
    q_outer = q_inner
    
    if bc_type == BC_WALL:
        q_outer = apply_slip_wall(q_inner, nx, ny, params)
    elif bc_type == BC_CYLINDER_WALL:
        q_outer = apply_cylinder_wall(q_inner, x, y, params)
    elif bc_type == BC_FARFIELD:
        q_outer = apply_farfield(q_inner, nx, ny, t, ramp_time, params)
    elif bc_type == BC_INLET:
        q_outer = apply_inlet(q_inner, nx, ny, t, ramp_time, params)
    elif bc_type == BC_OUTLET:
        q_outer = apply_outlet(q_inner, nx, ny, params)
    elif bc_type == BC_EXTRAPOLATION:
        q_outer = apply_extrapolation(q_inner)
        
    return q_outer
