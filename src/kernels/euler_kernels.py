import warp as wp

# --- Constants ---
gamma = 1.4
rho_floor = 1.0e-5
p_floor = 1.0e-5

BC_WALL = 1
BC_FARFIELD = 2

# --- Equation Kernels (Funktionen, die in anderen Kernels genutzt werden) ---

@wp.func
def pressure(q: wp.vec4) -> wp.float32:
    """Calculates the pressure from conservative variables with safety checks."""
    rho = wp.max(q[0], rho_floor)
    rho_u = q[1]
    rho_v = q[2]
    E = q[3]
    
    # Kinetische Energie: 0.5 * rho * (u^2 + v^2) = 0.5 * (rho_u^2 + rho_v^2) / rho
    kin_energy = 0.5 * (rho_u*rho_u + rho_v*rho_v) / rho
    p = (gamma - 1.0) * (E - kin_energy)
    return wp.max(p, p_floor)

@wp.func
def flux_x(q: wp.vec4) -> wp.vec4:
    """Flux function F(q) in x-direction."""
    rho = wp.max(q[0], rho_floor)
    p = pressure(q)
    u = q[1] / rho
    
    return wp.vec4(q[1], q[1]*u + p, q[2]*u, (q[3] + p)*u)

@wp.func
def flux_y(q: wp.vec4) -> wp.vec4:
    """Flux function G(q) in y-direction."""
    rho = wp.max(q[0], rho_floor)
    p = pressure(q)
    v = q[2] / rho
    
    return wp.vec4(q[2], q[1]*v, q[2]*v + p, (q[3] + p)*v)

@wp.func
def get_max_wave_speed(q: wp.vec4, nx: wp.float32, ny: wp.float32) -> wp.float32:
    """Calculates the acoustic wave speed |u_n| + c."""
    rho = wp.max(q[0], rho_floor)
    p = pressure(q)
    c = wp.sqrt(gamma * p / rho)
    
    # Velocity normal to the face
    u_n = (q[1] * nx + q[2] * ny) / rho
    return wp.abs(u_n) + c

@wp.func
def rusanev_flux(q_l: wp.vec4, q_r: wp.vec4, nx: wp.float32, ny: wp.float32) -> wp.vec4:
    """Computes the Lax-Friedrichs / Rusanov numerical flux."""
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
    Computes -div(F) using the strong form for general curvilinear coordinates:
    div F = dF/dx + dG/dy
    
    Chain rule:
    dF/dx = (dF/dr * dr/dx) + (dF/ds * ds/dx)
    dG/dy = (dG/dr * dr/dy) + (dG/ds * ds/dy)
    
    rhs = - (dF_dx + dG_dy)
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

    # Matrix-Vector Multiplication: Sum over j
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

    # Apply Chain Rule
    dF_dx = dF_dr * dr_dx + dF_ds * ds_dx
    dG_dy = dG_dr * dr_dy + dG_ds * ds_dy
    
    # RHS update (Strong Form: rhs = -div F)
    rhs[e, i] = -(dF_dx + dG_dy) 


@wp.kernel
def compute_surface_term(
    q: wp.array(dtype=wp.vec4, ndim=2),        # State
    rhs: wp.array(dtype=wp.vec4, ndim=2),      # RHS to accumulate into
    connectivity: wp.array(dtype=wp.int32, ndim=2), # (NumElems, 4) -> NeighborID
    face_map: wp.array(dtype=wp.int32, ndim=2),# (4, Nfp) -> Node Indices on faces
    LIFT: wp.array(dtype=wp.float32, ndim=2),  # (Np, 4*Nfp) Lift Matrix
    face_geo_factors: wp.array(dtype=wp.float32, ndim=3), # (NumElements, 4, 3) -> nx, ny, J_surf
    J: wp.array(dtype=wp.float32, ndim=1),     # Volume Jacobian
    bc_mask: wp.array(dtype=wp.int32, ndim=2), # (NumElements, 4) -> BC Type ID
    Nfp: wp.int32                              # Number of face points
):
    """
    Computes the Flux Jump and Lifts it to the volume nodes.
    Iterates over elements, calculates fluxes on all 4 faces, adds to RHS.
    Supports unstructured meshes by reading normals and scaling from face_geo_factors.
    Handles Boundary Conditions via bc_mask.
    """
    e = wp.tid() # One thread per element

    # Pre-load Volume Jacobian
    vol_J = J[e]
    inv_J = 1.0 / vol_J

    # Loop over all 4 faces
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
                # Find matching node on neighbor face.
                # If I am face 1 (Right), neighbor sees me as face 3 (Left).
                # Matching node ordering reverses (k -> Nfp-1-k).
                neighbor_face = (face_idx + 2) % 4
                
                neighbor_node_idx = face_map[neighbor_face, k]
                q_outer = q[neighbor_e, neighbor_node_idx]
            else:
                # --- Boundary Face ---
                bc_type = bc_mask[e, face_idx]
                
                if bc_type == BC_WALL:
                    # Slip Wall: Mirror velocity vector across the wall (remove normal component)
                    # v_ghost = v - 2(v . n)n
                    
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
                    # Extrapolation (Zero-Gradient / Transmissive)
                    q_outer = q_inner
                
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
    Berechnet die maximale Wellengeschwindigkeit im gesamten Gebiet
    und speichert das Maximum in max_speed[0].
    """
    e, i = wp.tid()
    
    # Zustand laden
    val = q[e, i]
    rho = wp.max(val[0], rho_floor)
    
    # Primitive Variablen
    u = val[1] / rho
    v = val[2] / rho
    p = pressure(val)
    
    # Schallgeschwindigkeit c
    c = wp.sqrt(gamma * p / rho)
    
    # Betrag der Geschwindigkeit |u|
    vel_mag = wp.sqrt(u*u + v*v)
    
    # Wellengeschwindigkeit lambda = |u| + c
    wave_speed = vel_mag + c
    
    # Globales Maximum schreiben
    # Hinweis: atomic_max funktioniert in neueren Warp-Versionen auch für floats.
    wp.atomic_max(max_speed, 0, wave_speed)