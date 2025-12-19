import warp as wp

# --- Constants ---
gamma = 1.4
rho_floor = 1.0e-5
p_floor = 1.0e-5

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
    rx: wp.array(dtype=wp.float32, ndim=1),  # Metric dr/dx (per Element)
    sy: wp.array(dtype=wp.float32, ndim=1),  # Metric ds/dy (per Element)
    Np: wp.int32                             # Number of points per element
):
    """
    Computes -div(F) using the strong form:
    rhs = - (rx * Dr * F + sy * Ds * G)
    """
    e, i = wp.tid() # Element e, Node i

    # Load Metrics for this element
    # Für kartesische Gitter sind rx, sy konstant im Element
    dr_dx = rx[e]
    ds_dy = sy[e]

    # Wir berechnen die Ableitung an Node i als Vektorprodukt
    dF_dr = wp.vec4(0.0)
    dG_ds = wp.vec4(0.0)

    # Matrix-Vector Multiplication: Sum over j
    # (Kann optimiert werden mit Shared Memory, aber Warp macht das oft automatisch gut)
    for j in range(Np):
        q_val = q[e, j]
        
        # Fluxes at node j
        F_val = flux_x(q_val)
        G_val = flux_y(q_val)
        
        # Accumulate derivatives
        dF_dr += F_val * Dr[i, j]
        dG_ds += G_val * Ds[i, j]

    # Chain rule: dF/dx = dF/dr * dr/dx  (da dr/dy = 0 etc. bei Rechtecken)
    div_F = dF_dr * dr_dx + dG_ds * ds_dy
    
    # RHS update (Strong Form: rhs = -div F)
    # Wir benutzen atomic_sub, falls wir diesen Kernel mehrmals aufrufen würden, 
    # aber hier reicht direktes Zuweisen, wenn wir RHS vorher nullen.
    rhs[e, i] = -div_F 


@wp.kernel
def compute_surface_term(
    q: wp.array(dtype=wp.vec4, ndim=2),        # State
    rhs: wp.array(dtype=wp.vec4, ndim=2),      # RHS to accumulate into
    connectivity: wp.array(dtype=wp.int32, ndim=2), # (NumElems, 4) -> NeighborID
    face_map: wp.array(dtype=wp.int32, ndim=2),# (4, Nfp) -> Node Indices on faces
    LIFT: wp.array(dtype=wp.float32, ndim=2),  # (Np, 4*Nfp) Lift Matrix
    Js_x: wp.array(dtype=wp.float32, ndim=1),  # Surface Jacobian x-face (dy/2)
    Js_y: wp.array(dtype=wp.float32, ndim=1),  # Surface Jacobian y-face (dx/2)
    J: wp.array(dtype=wp.float32, ndim=1),     # Volume Jacobian
    Nfp: wp.int32                              # Number of face points
):
    """
    Computes the Flux Jump and Lifts it to the volume nodes.
    Iterates over elements, calculates fluxes on all 4 faces, adds to RHS.
    """
    e = wp.tid() # One thread per element

    # Pre-load Geometric Factors
    vol_J = J[e]
    inv_J = 1.0 / vol_J
    surf_J_x = Js_x[e] # For vertical faces (Left/Right) -> Length involves dy
    surf_J_y = Js_y[e] # For horizontal faces (Bottom/Top) -> Length involves dx

    # Loop over all 4 faces
    for face_idx in range(4):
        
        # Determine Face Normal and Surface Jacobian
        nx = 0.0
        ny = 0.0
        surf_J = 0.0
        
        if face_idx == 0:   # Bottom (y=-1)
            nx = 0.0; ny = -1.0; surf_J = surf_J_y
        elif face_idx == 1: # Right (x=+1)
            nx = 1.0; ny = 0.0; surf_J = surf_J_x
        elif face_idx == 2: # Top (y=+1)
            nx = 0.0; ny = 1.0; surf_J = surf_J_y
        elif face_idx == 3: # Left (x=-1)
            nx = -1.0; ny = 0.0; surf_J = surf_J_x

        # Get Neighbor Info
        # Connectivity speichert hier [neighbor_id, neighbor_face_index]
        # Vereinfachung: Wir nehmen an, connectivity ist (NumElems, 4) mit NeighborIDs
        # und das Gitter ist konform -> wir finden die Indices rechnerisch.
        
        neighbor_e = connectivity[e, face_idx] 
        
        # Loop over nodes on this face
        for k in range(Nfp):
            # Node Index on current element
            node_idx_local = face_map[face_idx, k]
            
            q_inner = q[e, node_idx_local]
            q_outer = q_inner # Default for boundary (e.g. wall)
            
            if neighbor_e >= 0:
                # Find matching node on neighbor face.
                # Standard Structured Mesh logic:
                # If I am face 1 (Right), neighbor sees me as face 3 (Left).
                # Matching node ordering usually reverses (k -> Nfp-1-k) in 2D to match coordinates.
                neighbor_face = (face_idx + 2) % 4
                #neighbor_node_idx = face_map[neighbor_face, Nfp - 1 - k]
                neighbor_node_idx = face_map[neighbor_face, k]
                q_outer = q[neighbor_e, neighbor_node_idx]
            
            # 1. Numerical Flux (F*)
            f_star = rusanev_flux(q_inner, q_outer, nx, ny)
            
            # 2. Normal Flux from interior (F_n)
            f_n = flux_x(q_inner) * nx + flux_y(q_inner) * ny
            
            # 3. Flux Jump (F* - F_n) scaled by Surface Jacobian
            flux_jump = (f_n - f_star) * surf_J
            
            # 4. LIFTing: Add contribution to ALL volume nodes
            # LIFT matrix maps surface node k (on face f) to volume node i
            # Row in LIFT: volume node i, Col in LIFT: face_idx * Nfp + k
            lift_col = face_idx * Nfp + k
            
            for i in range(q.shape[1]): # Iterate over all volume nodes (Np)
                lift_val = LIFT[i, lift_col]
                # Add to RHS: 1/J * LIFT * Jump
                val = lift_val * flux_jump * inv_J
                
                # Atomic add needed because multiple faces contribute to same node
                wp.atomic_add(rhs, e, i, val)

@wp.kernel
def rk_step(
    q_old: wp.array(dtype=wp.vec4, ndim=2),
    q_new: wp.array(dtype=wp.vec4, ndim=2),
    rhs: wp.array(dtype=wp.vec4, ndim=2),
    dt: wp.float32,
    a: wp.float32, # RK parameter for q_old
    b: wp.float32  # RK parameter for (q_new + dt*rhs)
):
    e, i = wp.tid()
    q_new[e, i] = a * q_old[e, i] + b * (q_new[e, i] + dt * rhs[e, i])

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

@wp.kernel
def rk_stage_1(
    q: wp.array(dtype=wp.vec4, ndim=2),
    rhs: wp.array(dtype=wp.vec4, ndim=2),
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)
):
    """
    Stage 1: Q(1) = Q_n + dt * RHS(Q_n)
    """
    e, i = wp.tid()
    q_out[e, i] = q[e, i] + dt * rhs[e, i]

@wp.kernel
def rk_stage_2(
    q: wp.array(dtype=wp.vec4, ndim=2),      # Q_n (Startzustand)
    q_1: wp.array(dtype=wp.vec4, ndim=2),    # Q(1) aus Stage 1
    rhs: wp.array(dtype=wp.vec4, ndim=2),    # RHS(Q(1))
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)   # Ziel für Q(2)
):
    """
    Stage 2: Q(2) = 3/4 * Q_n + 1/4 * (Q(1) + dt * RHS(Q(1)))
    """
    e, i = wp.tid()
    # 0.75 * Q_n + 0.25 * (Q_1 + dt * RHS)
    q_out[e, i] = 0.75 * q[e, i] + 0.25 * (q_1[e, i] + dt * rhs[e, i])

@wp.kernel
def rk_stage_3(
    q: wp.array(dtype=wp.vec4, ndim=2),      # Q_n
    q_2: wp.array(dtype=wp.vec4, ndim=2),    # Q(2) aus Stage 2
    rhs: wp.array(dtype=wp.vec4, ndim=2),    # RHS(Q(2))
    dt: wp.float32,
    q_out: wp.array(dtype=wp.vec4, ndim=2)   # Ziel für Q(n+1) (überschreibt oft self.Q)
    ):
    """
    Stage 3: Q(n+1) = 1/3 * Q_n + 2/3 * (Q(2) + dt * RHS(Q(2)))
    """
    e, i = wp.tid()
    # 1/3 * Q_n + 2/3 * (Q_2 + dt * RHS)
    q_out[e, i] = (1.0/3.0) * q[e, i] + (2.0/3.0) * (q_2[e, i] + dt* rhs[e, i])