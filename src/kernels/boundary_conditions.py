import warp as wp

# --- Boundary Condition Type Constants ---
# These must match the IDs used in the Solver configuration
BC_INTERNAL = 0
BC_WALL = 1      # Slip Wall (Euler)
BC_FARFIELD = 2  # Non-reflecting / Outflow

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
def apply_farfield(q_inner: wp.vec4) -> wp.vec4:
    """
    Applies Farfield / Outflow boundary condition.
    Uses zeroth-order extrapolation (q_outer = q_inner).
    
    Args:
        q_inner (wp.vec4): State inside the domain.
        
    Returns:
        wp.vec4: The ghost state q_outer.
    """
    return q_inner

@wp.func
def apply_boundary_condition(bc_type: wp.int32, q_inner: wp.vec4, nx: wp.float32, ny: wp.float32) -> wp.vec4:
    """
    Dispatcher for boundary conditions.
    
    Args:
        bc_type (int): The boundary condition ID.
        q_inner (wp.vec4): State inside.
        nx, ny (float): Normal vector.
        
    Returns:
        wp.vec4: The ghost state q_outer.
    """
    # Default to inner state (transmissive)
    q_outer = q_inner
    
    if bc_type == BC_WALL:
        q_outer = apply_slip_wall(q_inner, nx, ny)
    elif bc_type == BC_FARFIELD:
        q_outer = apply_farfield(q_inner)
        
    return q_outer
