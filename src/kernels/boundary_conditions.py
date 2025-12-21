import warp as wp

# --- Constants ---
# These should ideally match physics/equations.py, but for kernels we define them here
gamma = 1.4

# --- Boundary Condition Type Constants ---
BC_INTERNAL = 0
BC_WALL = 1      # Slip Wall (Euler)
BC_FARFIELD = 2  # Farfield / Symmetry
BC_INLET = 3     # Inlet (Freestream)
BC_OUTLET = 4    # Subsonic Outlet

# --- Helper Functions ---

@wp.func
def get_freestream_state(t: wp.float32, ramp_time: wp.float32) -> wp.vec4:
    """Returns the freestream state with ramping (rho=1, u=ramped, v=0, p=1)."""
    rho = 1.0
    
    # Ramp u from 0 to 1 over ramp_time
    target_u = 1.0
    factor = 1.0
    if t < ramp_time:
        factor = t / ramp_time
        
    u = target_u * factor
    v = 0.0
    p = 1.0
    
    rho_u = rho * u
    rho_v = rho * v
    kinetic_energy = 0.5 * rho * (u*u + v*v)
    E = p / (gamma - 1.0) + kinetic_energy
    
    return wp.vec4(rho, rho_u, rho_v, E)

# --- Boundary Condition Functions ---

@wp.func
def apply_slip_wall(q_inner: wp.vec4, nx: wp.float32, ny: wp.float32) -> wp.vec4:
    """
    Applies Slip Wall boundary condition.
    Mirrors the velocity vector across the wall surface (removes normal component).
    
    Args:
        q_inner (wp.vec4): State inside the domain.
        nx, ny (float): Normal vector pointing OUT of the domain.
        
    Returns:
        wp.vec4: The ghost state q_outer.
    """
    rho = q_inner[0]
    rhou = q_inner[1]
    rhov = q_inner[2]
    E = q_inner[3]
    
    # Momentum dot Normal
    mom_dot_n = rhou * nx + rhov * ny
    
    # Reflected momentum: v_ghost = v - 2 * (v . n) * n
    rhou_ghost = rhou - 2.0 * mom_dot_n * nx
    rhov_ghost = rhov - 2.0 * mom_dot_n * ny
    
    return wp.vec4(rho, rhou_ghost, rhov_ghost, E)

@wp.func
def apply_farfield(q_inner: wp.vec4, nx: wp.float32, ny: wp.float32) -> wp.vec4:
    """
    Applies Farfield / Outflow boundary condition.
    
    Current implementation treats Farfield as Slip Wall (Symmetry) for stability 
    on parallel boundaries (preventing grazing flow instabilities).
    
    Args:
        q_inner (wp.vec4): State inside.
        nx, ny (float): Normal vector.
        
    Returns:
        wp.vec4: Ghost state.
    """
    # Use Slip Wall logic for stability
    return apply_slip_wall(q_inner, nx, ny)

@wp.func
def apply_inlet(t: wp.float32, ramp_time: wp.float32) -> wp.vec4:
    """
    Applies Inlet boundary condition (Dirichlet Freestream).
    """
    return get_freestream_state(t, ramp_time)

@wp.func
def apply_outlet(q_inner: wp.vec4) -> wp.vec4:
    """
    Applies Subsonic Outlet boundary condition.
    Fixes pressure to p_back=1.0, extrapolates density and velocity.
    """
    # p_back = 1.0
    rho = q_inner[0]
    rhou = q_inner[1]
    rhov = q_inner[2]
    
    # Compute inner velocities to reconstruct Energy with new Pressure
    p_back = float(1.0)
    
    # Kinetic Energy = 0.5 * (rhou^2 + rhov^2) / rho
    # Note: If rho is very small, this could be unstable, but q_inner should be valid.
    kin_energy = 0.5 * (rhou*rhou + rhov*rhov) / rho
    E_outer = p_back / (gamma - 1.0) + kin_energy
    
    return wp.vec4(rho, rhou, rhov, E_outer)

@wp.func
def apply_boundary_condition(
    bc_type: wp.int32, 
    q_inner: wp.vec4, 
    nx: wp.float32, 
    ny: wp.float32,
    t: wp.float32,
    ramp_time: wp.float32
) -> wp.vec4:
    """
    Dispatcher for boundary conditions.
    
    Args:
        bc_type (int): The boundary condition ID.
        q_inner (wp.vec4): State inside.
        nx, ny (float): Normal vector.
        t (float): Current simulation time.
        ramp_time (float): Ramping parameter.
        
    Returns:
        wp.vec4: The ghost state q_outer.
    """
    # Default to inner state (transmissive/extrapolation)
    q_outer = q_inner
    
    if bc_type == BC_WALL:
        q_outer = apply_slip_wall(q_inner, nx, ny)
    elif bc_type == BC_FARFIELD:
        q_outer = apply_farfield(q_inner, nx, ny)
    elif bc_type == BC_INLET:
        q_outer = apply_inlet(t, ramp_time)
    elif bc_type == BC_OUTLET:
        q_outer = apply_outlet(q_inner)
        
    return q_outer