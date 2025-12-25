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
    rho = q_inner[0]
    # Protect against vacuum
    if rho < 1e-6:
        return get_freestream_state(t, ramp_time, q_inner, params)
        
    u = q_inner[1] / rho
    v = q_inner[2] / rho
    
    vn = u * nx + v * ny
    
    if vn > 0.0:
        return q_inner
    else:
        return get_freestream_state(t, ramp_time, q_inner, params)

@wp.func
def apply_inlet(t: Any, ramp_time: Any, template_q: Any, params: Any):
    return get_freestream_state(t, ramp_time, template_q, params)

@wp.func
def apply_outlet(q_inner: Any, params: Any):
    """
    Applies Subsonic Outlet boundary condition.
    """
    rho = q_inner[0]
    rhou = q_inner[1]
    rhov = q_inner[2]
    
    p_back = rho - rho + params.one
    
    kin_energy = params.half * (rhou*rhou + rhov*rhov) / rho
    E_outer = p_back / (params.gamma - params.one) + kin_energy
    
    return make_vec4_generic(rho, rhou, rhov, E_outer)

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
        q_outer = apply_inlet(t, ramp_time, q_inner, params)
    elif bc_type == BC_OUTLET:
        q_outer = apply_outlet(q_inner, params)
    elif bc_type == BC_EXTRAPOLATION:
        q_outer = apply_extrapolation(q_inner)
        
    return q_outer
