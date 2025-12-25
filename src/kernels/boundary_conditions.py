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
BC_DOUBLE_MACH_EXACT = 7 # Exact shock motion for DMR boundaries

# --- Helper Functions ---

@wp.func
def make_vec4_generic(x: wp.float32, y: wp.float32, z: wp.float32, w: wp.float32):
    return wp.vec4(x, y, z, w)

@wp.func
def make_vec4_generic(x: wp.float64, y: wp.float64, z: wp.float64, w: wp.float64):
    return wp.vec4d(x, y, z, w)

@wp.func
def get_half_generic(template: wp.float32):
    return wp.float32(0.5)

@wp.func
def get_half_generic(template: wp.float64):
    return wp.float64(0.5)

@wp.func
def get_one_generic(template: wp.float32):
    return wp.float32(1.0)

@wp.func
def get_one_generic(template: wp.float64):
    return wp.float64(1.0)

@wp.func
def get_any_generic(template: wp.float32, val: wp.float32):
    return wp.float32(val)

@wp.func
def get_any_generic(template: wp.float32, val: wp.float64):
    return wp.float32(val)

@wp.func
def get_any_generic(template: wp.float64, val: wp.float32):
    return wp.float64(val)

@wp.func
def get_any_generic(template: wp.float64, val: wp.float64):
    return wp.float64(val)

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
    """Returns the freestream state with ramping."""
    zero = template_q[0] - template_q[0]
    
    rho = zero + params.rho_inf
    
    # Ramp velocity if needed
    factor = zero + params.one
    if t < ramp_time:
        factor = t / ramp_time
        
    u = params.u_inf * factor
    v = params.v_inf * factor
    p = zero + params.p_inf
    
    rho_u = rho * u
    rho_v = rho * v
    kinetic_energy = params.half * rho * (u*u + v*v)
    E = p / (params.gamma - params.one) + kinetic_energy
    
    return make_vec4_generic(rho, rho_u, rho_v, E)

@wp.func
def compute_characteristic_state(q_inner: Any, q_target: Any, nx: Any, ny: Any, params: Any):
    """
    Computes the boundary state using Riemann invariants (Characteristic Boundary Conditions).
    A simplified version that is more robust against reflections.
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
    
    # --- Riemann Invariants ---
    gm1 = params.gamma - params.one
    inv_gm1 = params.one / gm1
    two = params.one + params.one
    half = get_half_generic(nx)
    
    # J+ (outgoing from inner), J- (incoming from target)
    j_plus = un_i + two * c_i * inv_gm1
    j_minus = un_t - two * c_t * inv_gm1
    
    # --- Flow Regime ---
    zero = params.one - params.one
    
    if un_i > c_i: # Supersonic Outflow
        return q_inner
    elif un_i < -c_i: # Supersonic Inflow
        return q_target
    else:
        # Subsonic Outflow or Inflow
        un_b = half * (j_plus + j_minus)
        c_b = half * half * gm1 * (j_plus - j_minus)
        
        # Determine if we use entropy/tangential from inner or target
        if un_i >= zero: # Outflow
            # Tangential velocity and entropy from inner
            ut_i = -u_i * ny + v_i * nx
            s_i = p_i / wp.pow(rho_i, params.gamma)
            
            rho_b = wp.pow((c_b * c_b) / (params.gamma * s_i), inv_gm1)
            p_b = s_i * wp.pow(rho_b, params.gamma)
            
            u_b = un_b * nx - ut_i * ny
            v_b = un_b * ny + ut_i * nx
        else: # Inflow
            # Tangential velocity and entropy from target
            ut_t = -u_t * ny + v_t * nx
            s_t = p_t / wp.pow(rho_t, params.gamma)
            
            rho_b = wp.pow((c_b * c_b) / (params.gamma * s_t), inv_gm1)
            p_b = s_t * wp.pow(rho_b, params.gamma)
            
            u_b = un_b * nx - ut_t * ny
            v_b = un_b * ny + ut_t * nx
            
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
def apply_farfield(q_inner: Any, nx: Any, ny: Any, bc_state: Any, params: Any):
    """
    Applies Characteristic Farfield boundary condition using state from bc_state.
    """
    zero = params.one - params.one
    
    # Calculate energy for target state
    kinetic_energy = params.half * bc_state.v0 * (bc_state.v1*bc_state.v1 + bc_state.v2*bc_state.v2)
    E_target = bc_state.v3 / (params.gamma - params.one) + kinetic_energy
    q_target = make_vec4_generic(bc_state.v0, bc_state.v0 * bc_state.v1, bc_state.v0 * bc_state.v2, E_target)

    return compute_characteristic_state(q_inner, q_target, nx, ny, params)

@wp.func
def apply_inlet(q_inner: Any, nx: Any, ny: Any, t: Any, ramp_time: Any, bc_state: Any, params: Any):
    """
    Applies Inlet boundary condition using characteristics and state from bc_state.
    """
    # Create target state from bc_state parameters with ramping
    rho = bc_state.v0
    
    factor = nx - nx + params.one
    if t < ramp_time:
        factor = t / ramp_time
        
    u = bc_state.v1 * factor
    v = bc_state.v2 * factor
    p = bc_state.v3
    
    rho_u = rho * u
    rho_v = rho * v
    kinetic_energy = params.half * rho * (u*u + v*v)
    E = p / (params.gamma - params.one) + kinetic_energy
    
    q_inlet = make_vec4_generic(rho, rho_u, rho_v, E)
    return compute_characteristic_state(q_inner, q_inlet, nx, ny, params)

@wp.func
def apply_outlet(q_inner: Any, nx: Any, ny: Any, bc_state: Any, params: Any):
    """
    Applies Subsonic Outlet boundary condition using characteristics and p_back from bc_state.
    """
    # Target state: rho, u, v from inner, but p = p_back from bc_state.v0
    rho = wp.max(params.rho_floor, q_inner[0])
    inv_rho = params.one / rho
    u = q_inner[1] * inv_rho
    v = q_inner[2] * inv_rho
    
    p_back = bc_state.v0
    
    kin = params.half * rho * (u*u + v*v)
    E_target = p_back / (params.gamma - params.one) + kin
    
    q_target = make_vec4_generic(rho, q_inner[1], q_inner[2], E_target)
    
    return compute_characteristic_state(q_inner, q_target, nx, ny, params)

@wp.func
def apply_double_mach_exact(q_inner: Any, nx: Any, ny: Any, x: Any, y: Any, t: Any, params: Any):
    """
    Exact shock motion for DMR boundaries.
    """
    # Use template for correct precision constants
    template = nx
    half = get_half_generic(template)
    one = get_one_generic(template)
    
    # Shock line: x = 1/6 + y/tan(60) + (10/sin(60))*t
    # sin(60) = sqrt(3)/2, tan(60) = sqrt(3)
    three = get_any_generic(template, 3.0)
    sqrt3 = wp.sqrt(three)
    sin60 = sqrt3 * half
    
    six = get_any_generic(template, 6.0)
    ten = get_any_generic(template, 10.0)
    shock_x = one / six + y / sqrt3 + (ten / sin60) * t
    
    if x < shock_x:
        # Post-shock (Left) State: rho=8, p=116.5, u=8.25*cos(30), v=-8.25*sin(30)
        rho = get_any_generic(template, 8.0)
        u825 = get_any_generic(template, 8.25)
        # cos(30) = sin(60)
        u = u825 * sin60
        # sin(30) = 0.5
        v = -u825 * half
        p = get_any_generic(template, 116.5)
    else:
        # Pre-shock (Right) State: rho=1.4, u=0, v=0, p=1.0
        rho = get_any_generic(template, 1.4)
        u = get_any_generic(template, 0.0)
        v = get_any_generic(template, 0.0)
        p = one
        
    rho_u = rho * u
    rho_v = rho * v
    kin = half * rho * (u*u + v*v)
    E = p / (params.gamma - params.one) + kin
    
    q_target = make_vec4_generic(rho, rho_u, rho_v, E)
    return compute_characteristic_state(q_inner, q_target, nx, ny, params)

@wp.func
def apply_boundary_condition(
    bc_index: wp.int32,
    bc_data: Any,
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
    bc_state = bc_data[bc_index]
    bc_type = bc_state.type
    q_outer = q_inner
    
    if bc_type == BC_WALL:
        q_outer = apply_slip_wall(q_inner, nx, ny, params)
    elif bc_type == BC_CYLINDER_WALL:
        q_outer = apply_cylinder_wall(q_inner, x, y, params)
    elif bc_type == BC_FARFIELD:
        q_outer = apply_farfield(q_inner, nx, ny, bc_state, params)
    elif bc_type == BC_INLET:
        q_outer = apply_inlet(q_inner, nx, ny, t, ramp_time, bc_state, params)
    elif bc_type == BC_OUTLET:
        q_outer = apply_outlet(q_inner, nx, ny, bc_state, params)
    elif bc_type == BC_EXTRAPOLATION:
        q_outer = apply_extrapolation(q_inner)
    elif bc_type == BC_DOUBLE_MACH_EXACT:
        q_outer = apply_double_mach_exact(q_inner, nx, ny, x, y, t, params)
        
    return q_outer
