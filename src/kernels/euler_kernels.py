import warp as wp

# --- Constants ---
gamma = 1.4
rho_floor = 1.0e-5
p_floor = 1.0e-5

BC_WALL = 1
BC_FARFIELD = 2
BC_INLET = 3
BC_OUTLET = 4

# --- Equation Kernels (Helper functions used by other kernels) ---

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

@wp.func
def pressure(q: wp.vec4) -> wp.float32:
    """
    Calculates the pressure from conservative variables with safety checks.
    
    p = (gamma - 1) * (E - 0.5 * rho * |u|^2)

    Args:
        q (wp.vec4): Conservative state vector [rho, rho*u, rho*v, E].

    Returns:
        float: The pressure, clamped to p_floor.
    """
    rho = wp.max(q[0], rho_floor)
    rho_u = q[1]
    rho_v = q[2]
    E = q[3]
    
    # Kinetic Energy: 0.5 * rho * (u^2 + v^2) = 0.5 * (rho_u^2 + rho_v^2) / rho
    kin_energy = 0.5 * (rho_u*rho_u + rho_v*rho_v) / rho
    p = (gamma - 1.0) * (E - kin_energy)
    return wp.max(p, p_floor)

@wp.func
def flux_x(q: wp.vec4) -> wp.vec4:
    """
    Computes the Euler Flux function F(q) in the x-direction.
    
    F(q) = [rho*u, rho*u^2 + p, rho*u*v, (E+p)*u]

    Args:
        q (wp.vec4): Conservative state vector.

    Returns:
        wp.vec4: The flux vector in x-direction.
    """
    rho = wp.max(q[0], rho_floor)
    p = pressure(q)
    u = q[1] / rho
    
    return wp.vec4(q[1], q[1]*u + p, q[2]*u, (q[3] + p)*u)

@wp.func
def flux_y(q: wp.vec4) -> wp.vec4:
    """
    Computes the Euler Flux function G(q) in the y-direction.
    
    G(q) = [rho*v, rho*v*u, rho*v^2 + p, (E+p)*v]

    Args:
        q (wp.vec4): Conservative state vector.

    Returns:
        wp.vec4: The flux vector in y-direction.
    """
    rho = wp.max(q[0], rho_floor)
    p = pressure(q)
    v = q[2] / rho
    
    return wp.vec4(q[2], q[1]*v, q[2]*v + p, (q[3] + p)*v)

@wp.func
def get_max_wave_speed(q: wp.vec4, nx: wp.float32, ny: wp.float32) -> wp.float32:
    """
    Calculates the acoustic wave speed |u_n| + c in the direction normal to a face.

    Args:
        q (wp.vec4): Conservative state vector.
        nx (float): X-component of the normal vector.
        ny (float): Y-component of the normal vector.

    Returns:
        float: The maximum wave speed.
    """
    rho = wp.max(q[0], rho_floor)
    p = pressure(q)
    c = wp.sqrt(gamma * p / rho)
    
    # Velocity normal to the face
    u_n = (q[1] * nx + q[2] * ny) / rho
    return wp.abs(u_n) + c

@wp.func
def rusanev_flux(q_l: wp.vec4, q_r: wp.vec4, nx: wp.float32, ny: wp.float32) -> wp.vec4:
    """
    Computes the Lax-Friedrichs / Rusanov numerical flux across an interface.
    
    F* = 0.5 * (F_l + F_r) - 0.5 * alpha * (Q_r - Q_l)
    where alpha is the maximum wave speed.

    Args:
        q_l (wp.vec4): State on the left/inner side of the face.
        q_r (wp.vec4): State on the right/outer side of the face.
        nx (float): Normal x-component.
        ny (float): Normal y-component.

    Returns:
        wp.vec4: The numerical flux vector.
    """
    # Fluxes projected onto normal
    F_l = flux_x(q_l) * nx + flux_y(q_l) * ny
    F_r = flux_x(q_r) * nx + flux_y(q_r) * ny
    
    # Wave speeds
    lambda_l = get_max_wave_speed(q_l, nx, ny)
    lambda_r = get_max_wave_speed(q_r, nx, ny)
    alpha = wp.max(lambda_l, lambda_r)
    
    # Numerical Flux: Avg(Flux) - Dissipation
    return 0.5 * (F_l + F_r) - 0.5 * alpha * (q_r - q_l)

# --- Solver Kernels ---

@wp.kernel
def compute_volume_term(
    q: wp.array(dtype=wp.vec4, ndim=2),      # (NumElems, Np)
    rhs: wp.array(dtype=wp.vec4, ndim=2),    # Output RHS
    Dr: wp.array(dtype=wp.float32, ndim=2),  # Differentiation Matrix r
    Ds: wp.array(dtype=wp.float32, ndim=2),  # Differentiation Matrix s
    rx: wp.array(dtype=wp.float32, ndim=1),  # Metric dr/dx
    ry: wp.array(dtype=wp.float32, ndim=1),  # Metric dr/dy
    sx: wp.array(dtype=wp.float32, ndim=1),  # Metric ds/dx
    sy: wp.array(dtype=wp.float32, ndim=1),  # Metric ds/dy
    Np: wp.int32                             # Number of points per element
):
    """
    Computes the divergence of the flux (volume integral) using the strong form for general curvilinear coordinates.
    
    div F = dF/dx + dG/dy
    Using the chain rule:
    dF/dx = (dF/dr * dr/dx) + (dF/ds * ds/dx)
    dG/dy = (dG/dr * dr/dy) + (dG/ds * ds/dy)
    
    The result is added to the RHS (Strong Form: rhs = -div F).

    Args:
        q (wp.array): State vector array.
        rhs (wp.array): Right-hand side accumulation array.
        Dr (wp.array): Differentiation matrix for reference coordinate r.
        Ds (wp.array): Differentiation matrix for reference coordinate s.
        rx (wp.array): Metric term dr/dx per element.
        ry (wp.array): Metric term dr/dy per element.
        sx (wp.array): Metric term ds/dx per element.
        sy (wp.array): Metric term ds/dy per element.
        Np (int): Number of nodes per element.
    """
    e, i = wp.tid() # Element e, Node i

    # Load Metrics for this element (constant for affine elements)
    dr_dx = rx[e]
    dr_dy = ry[e]
    ds_dx = sx[e]
    ds_dy = sy[e]

    dF_dr = wp.vec4(0.0)
    dF_ds = wp.vec4(0.0)
    dG_dr = wp.vec4(0.0)
    dG_ds = wp.vec4(0.0)

    # Matrix-Vector Multiplication: Sum over j to compute derivatives in reference space
    for j in range(Np):
        q_val = q[e, j]
        
        # Fluxes at node j
        F_val = flux_x(q_val)
        G_val = flux_y(q_val)
        
        # Accumulate derivatives
        # We need derivatives of both fluxes with respect to both ref coords
        dr = Dr[i, j]
        ds = Ds[i, j]
        
        dF_dr += F_val * dr
        dF_ds += F_val * ds
        
        dG_dr += G_val * dr
        dG_ds += G_val * ds

    # Apply Chain Rule to map derivatives to physical space
    dF_dx = dF_dr * dr_dx + dF_ds * ds_dx
    dG_dy = dG_dr * dr_dy + dG_ds * ds_dy
    
    # RHS update (Strong Form: rhs = -div F)
    rhs[e, i] = -(dF_dx + dG_dy) 


@wp.kernel
def compute_surface_term(
    q: wp.array(dtype=wp.vec4, ndim=2),        # State
    rhs: wp.array(dtype=wp.vec4, ndim=2),      # RHS to accumulate into
    connectivity: wp.array(dtype=wp.int32, ndim=2), # (NumElems, 4) -> NeighborID
    neighbor_face_indices: wp.array(dtype=wp.int32, ndim=2), # (NumElems, 4) -> NeighborFaceID
    face_map: wp.array(dtype=wp.int32, ndim=2),# (4, Nfp) -> Node Indices on faces
    LIFT: wp.array(dtype=wp.float32, ndim=2),  # (Np, 4*Nfp) Lift Matrix
    face_geo_factors: wp.array(dtype=wp.float32, ndim=3), # (NumElements, 4, 3) -> nx, ny, J_surf
    J: wp.array(dtype=wp.float32, ndim=1),     # Volume Jacobian
    bc_mask: wp.array(dtype=wp.int32, ndim=2), # (NumElements, 4) -> BC Type ID
    Nfp: wp.int32,                             # Number of face points
    t: wp.float32,                             # Current simulation time
    ramp_time: wp.float32                      # Time to ramp up inlet
):
    """
    Computes the surface integral (flux jump) and lifts it to the volume nodes.
    
    Iterates over all elements and their 4 faces. Calculates the numerical flux at the interface,
    computes the jump against the internal normal flux, and projects this jump to the volume nodes
    using the LIFT matrix. Supports unstructured meshes and boundary conditions.

    Args:
        q (wp.array): State vector array.
        rhs (wp.array): Right-hand side accumulation array.
        connectivity (wp.array): Element connectivity (neighbors).
        neighbor_face_indices (wp.array): Neighbor face indices.
        face_map (wp.array): Map from face index to local node indices.
        LIFT (wp.array): Lift operator matrix.
        face_geo_factors (wp.array): Geometric factors for faces (normals, surface Jacobian).
        J (wp.array): Volume Jacobian determinant per element.
        bc_mask (wp.array): Boundary condition type ID per face.
        Nfp (int): Number of face points.
    """
    e = wp.tid() # One thread per element

    # Pre-load Volume Jacobian
    vol_J = J[e]
    inv_J = 1.0 / vol_J

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
                
                # Assume standard conforming mesh orientation (reversed parameterization)
                # My node k corresponds to neighbor node (Nfp - 1 - k)
                neighbor_k = Nfp - 1 - k
                neighbor_node_idx = face_map[neighbor_face, neighbor_k]
                
                q_outer = q[neighbor_e, neighbor_node_idx]
            else:
                # --- Boundary Face ---
                bc_type = bc_mask[e, face_idx]
                
                if bc_type == BC_WALL:
                    # Slip Wall: Mirror velocity vector across the wall (remove normal component)
                    
                    rho = q_inner[0]
                    rhou = q_inner[1]
                    rhov = q_inner[2]
                    E = q_inner[3]
                    
                    # Momentum dot Normal
                    mom_dot_n = rhou * nx + rhov * ny
                    
                    rhou_ghost = rhou - 2.0 * mom_dot_n * nx
                    rhov_ghost = rhov - 2.0 * mom_dot_n * ny
                    
                    q_outer = wp.vec4(rho, rhou_ghost, rhov_ghost, E)
                    
                elif bc_type == BC_FARFIELD:
                    # Treat Farfield as Slip Wall (Symmetry) for stability on parallel boundaries.
                    # This enforces zero normal velocity, preventing inflow/outflow instabilities
                    # typical of grazing flow with simple extrapolation or Dirichlet BCs.
                    
                    rho = q_inner[0]
                    rhou = q_inner[1]
                    rhov = q_inner[2]
                    E = q_inner[3]
                    
                    # Momentum dot Normal
                    mom_dot_n = rhou * nx + rhov * ny
                    
                    # Reflected momentum: v_ghost = v - 2 * (v . n) * n
                    rhou_ghost = rhou - 2.0 * mom_dot_n * nx
                    rhov_ghost = rhov - 2.0 * mom_dot_n * ny
                    
                    q_outer = wp.vec4(rho, rhou_ghost, rhov_ghost, E)
                
                elif bc_type == BC_INLET:
                    # Dirichlet (Freestream)
                    q_outer = get_freestream_state(t, ramp_time)

                elif bc_type == BC_OUTLET:
                    # Subsonic Outlet: Fix Pressure, Extrapolate others
                    # p_back = 1.0
                    rho = q_inner[0]
                    rhou = q_inner[1]
                    rhov = q_inner[2]
                    
                    # Compute inner velocities to reconstruct Energy with new Pressure
                    # Kinetic Energy
                    # q_inner[3] is E_inner, we don't need it directly if we recompute E
                    
                    # E_outer = p_back / (gamma - 1) + 0.5 * (rho*u^2 + rho*v^2)
                    # 0.5 * rho * V^2 = 0.5 * (rhou^2 + rhov^2) / rho
                    
                    p_back = float(1.0)
                    kin_energy = 0.5 * (rhou*rhou + rhov*rhov) / rho
                    E_outer = p_back / (gamma - 1.0) + kin_energy
                    
                    q_outer = wp.vec4(rho, rhou, rhov, E_outer)

                else:
                    # Default/Fallback
                    q_outer = q_inner
            
            # 1. Numerical Flux (F*)
            f_star = rusanev_flux(q_inner, q_outer, nx, ny)
            
            # 2. Normal Flux from interior (F_n)
            f_n = flux_x(q_inner) * nx + flux_y(q_inner) * ny
            
            # 3. Flux Jump (F* - F_n) scaled by Surface Jacobian
            flux_jump = (f_n - f_star) * surf_J
            
            # 4. LIFTing: Add contribution to ALL volume nodes
            lift_col = face_idx * Nfp + k
            
            for i in range(q.shape[1]): # Iterate over all volume nodes (Np)
                lift_val = LIFT[i, lift_col]
                # Add to RHS: 1/J * LIFT * Jump
                val = lift_val * flux_jump * inv_J
                
                wp.atomic_add(rhs, e, i, val)

@wp.kernel
def compute_max_wave_speed(
    q: wp.array(dtype=wp.vec4, ndim=2),      # Shape: (num_elements, Np)
    max_speed: wp.array(dtype=wp.float32, ndim=1) # Shape: (1,)
):
    """
    Computes the maximum wave speed in the entire domain for CFL calculation.
    
    This kernel iterates over all nodes, calculates |u| + c, and updates a global maximum using atomic operations.

    Args:
        q (wp.array): State vector array.
        max_speed (wp.array): A scalar array (size 1) to store the result.
    """
    e, i = wp.tid()
    
    # Load state
    val = q[e, i]
    rho = wp.max(val[0], rho_floor)
    
    # Primitive variables
    u = val[1] / rho
    v = val[2] / rho
    p = pressure(val)
    
    # Sound speed c
    c = wp.sqrt(gamma * p / rho)
    
    # Velocity magnitude |u|
    vel_mag = wp.sqrt(u*u + v*v)
    
    # Wave speed lambda = |u| + c
    wave_speed = vel_mag + c
    
    # Update global maximum atomically
    wp.atomic_max(max_speed, 0, wave_speed)