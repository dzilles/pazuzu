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
def get_freestream_state(t: Any, ramp_time: Any, template_q: Any, gamma: Any, half: Any, one: Any):
    """Returns the freestream state with ramping (rho=1, u=ramped, v=0, p=1)."""
    zero = template_q[0] - template_q[0]
    
    rho = zero + one
    
    # Ramp u from 0 to 1 over ramp_time
    target_u = zero + one
    factor = zero + one
    
    if t < ramp_time:
        factor = t / ramp_time
        
    u = target_u * factor
    v = zero
    p = zero + one
    
    rho_u = rho * u
    rho_v = rho * v
    kinetic_energy = half * rho * (u*u + v*v)
    E = p / (gamma - one) + kinetic_energy
    
    return make_vec4_generic(rho, rho_u, rho_v, E)

# --- Boundary Condition Functions ---

@wp.func
def apply_slip_wall(q_inner: Any, nx: Any, ny: Any, one: Any):
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
    two = one + one
    rhou_ghost = rhou - two * mom_dot_n * nx
    rhov_ghost = rhov - two * mom_dot_n * ny
    
    return make_vec4_generic(rho, rhou_ghost, rhov_ghost, E)

@wp.func
def apply_cylinder_wall(q_inner: Any, x: Any, y: Any, one: Any):
    """
    Applies Slip Wall boundary condition using analytical normals for a cylinder centered at (0,0).
    """
    # Calculate analytical normal for cylinder centered at (0,0)
    radius = wp.sqrt(x*x + y*y)
    
    # Avoid division by zero
    nx_analytisch = x - x + one
    ny_analytisch = x - x
    if radius > 1e-8:
        nx_analytisch = x / radius
        ny_analytisch = y / radius
        
    return apply_slip_wall(q_inner, nx_analytisch, ny_analytisch, one)

@wp.func
def apply_extrapolation(q_inner: Any):
    return q_inner

@wp.func
def apply_farfield(q_inner: Any, nx: Any, ny: Any, t: Any, ramp_time: Any, gamma: Any, half: Any, one: Any):
    """
    Applies Characteristic Farfield boundary condition.
    """
    rho = q_inner[0]
    # Protect against vacuum
    if rho < 1e-6:
        return get_freestream_state(t, ramp_time, q_inner, gamma, half, one)
        
    u = q_inner[1] / rho
    v = q_inner[2] / rho
    
    vn = u * nx + v * ny
    
    if vn > 0.0:
        return q_inner
    else:
        return get_freestream_state(t, ramp_time, q_inner, gamma, half, one)

@wp.func
def apply_inlet(t: Any, ramp_time: Any, template_q: Any, gamma: Any, half: Any, one: Any):
    return get_freestream_state(t, ramp_time, template_q, gamma, half, one)

@wp.func
def apply_outlet(q_inner: Any, gamma: Any, half: Any, one: Any):
    """
    Applies Subsonic Outlet boundary condition.
    """
    rho = q_inner[0]
    rhou = q_inner[1]
    rhov = q_inner[2]
    
    p_back = rho - rho + one
    
    kin_energy = half * (rhou*rhou + rhov*rhov) / rho
    E_outer = p_back / (gamma - one) + kin_energy
    
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
    gamma: Any,
    half: Any,
    one: Any
):
    """
    Dispatcher for boundary conditions.
    """
    q_outer = q_inner
    
    if bc_type == BC_WALL:
        q_outer = apply_slip_wall(q_inner, nx, ny, one)
    elif bc_type == BC_CYLINDER_WALL:
        q_outer = apply_cylinder_wall(q_inner, x, y, one)
    elif bc_type == BC_FARFIELD:
        q_outer = apply_farfield(q_inner, nx, ny, t, ramp_time, gamma, half, one)
    elif bc_type == BC_INLET:
        q_outer = apply_inlet(t, ramp_time, q_inner, gamma, half, one)
    elif bc_type == BC_OUTLET:
        q_outer = apply_outlet(q_inner, gamma, half, one)
    elif bc_type == BC_EXTRAPOLATION:
        q_outer = apply_extrapolation(q_inner)
        
    return q_outer
