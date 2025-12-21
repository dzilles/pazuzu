import numpy as np

# Gas constant
gamma = 1.4

def pressure(Q):
    """
    Computes the pressure from the conservative state vector using the ideal gas law.
    
    Equation of state: p = (gamma - 1) * (E - 0.5 * rho * (u^2 + v^2))

    Args:
        Q (np.array): Conservative variables [rho, rho*u, rho*v, E].

    Returns:
        np.array: The pressure.
    """
    rho, rho_u, rho_v, E = Q
    # Soft clip to avoid division by zero (better than hard clip for optimization, but fine here)
    rho = np.maximum(rho, 1e-12)
    
    # Optimization: 0.5 * (rho_u^2 + rho_v^2) / rho is numerically more stable 
    # than 0.5 * rho * (u^2 + v^2), as u=rho_u/rho -> u^2 = rho_u^2/rho^2
    # This saves one division and preserves accuracy for small rho.
    kinetic_energy = 0.5 * (rho_u**2 + rho_v**2) / rho
    p_val = (gamma - 1) * (E - kinetic_energy)
    
    return np.maximum(p_val, 1e-12)

def euler_fluxes(Q):
    """
    Computes both Euler Flux vectors F(Q) and G(Q) simultaneously.
    
    F(Q) = [rho*u, rho*u^2 + p, rho*u*v, (E+p)*u]
    G(Q) = [rho*v, rho*v*u, rho*v^2 + p, (E+p)*v]

    Args:
        Q (np.array): Conservative variables.

    Returns:
        tuple: (F, G) numpy arrays containing the fluxes in x and y directions.
    """
    rho, rho_u, rho_v, E = Q
    rho = np.maximum(rho, 1e-12)
    p = pressure(Q)
    
    # Inverse density
    inv_rho = 1.0 / rho
    u = rho_u * inv_rho
    v = rho_v * inv_rho
    
    # Flux F (x-dir)
    f1 = rho_u
    f2 = rho_u * u + p
    f3 = rho_u * v     # rho_v * u is equivalent to rho_u * v
    f4 = (E + p) * u
    
    # Flux G (y-dir)
    g1 = rho_v
    g2 = rho_v * u     # rho_u * v
    g3 = rho_v * v + p
    g4 = (E + p) * v
    
    F = np.array([f1, f2, f3, f4])
    G = np.array([g1, g2, g3, g4])
    return F, G

def get_max_eigenvalue(Q, normal_vector):
    """
    Computes the maximum eigenvalue (wave speed) of the system in a given normal direction.
    
    lambda_max = |u_n| + c

    Args:
        Q (np.array): Conservative variables.
        normal_vector (tuple): (nx, ny) components of the normal vector.

    Returns:
        np.array: The maximum wave speed.
    """
    rho = Q[0]
    # Avoid division by near-zero
    rho = np.maximum(rho, 1e-12)
    p = pressure(Q)
    
    # Sound speed
    c = np.sqrt(gamma * p / rho)
    
    nx, ny = normal_vector
    # Velocity projected on normal
    # u_n = (rho_u * nx + rho_v * ny) / rho
    u_n = (Q[1] * nx + Q[2] * ny) / rho
    
    return np.abs(u_n) + c

def lax_friedrichs_flux(q_l, q_r, normal_vector):
    """
    Computes the Lax-Friedrichs / Rusanov numerical flux projected on the normal.
    
    F* = 0.5 * (F_L + F_R) - 0.5 * alpha * (Q_R - Q_L)

    Args:
        q_l (np.array): State on the left side of the interface.
        q_r (np.array): State on the right side of the interface.
        normal_vector (tuple): Normal vector (nx, ny).

    Returns:
        np.array: The numerical flux vector.
    """
    nx, ny = normal_vector
    
    # Compute fluxes
    # We use the combined function to save computation time
    F_l, G_l = euler_fluxes(q_l)
    F_r, G_r = euler_fluxes(q_r)
    
    # Projected fluxes (Flux dot Normal)
    flux_n_l = F_l * nx + G_l * ny
    flux_n_r = F_r * nx + G_r * ny
    
    # Max eigenvalue (wave speed)
    lambda_l = get_max_eigenvalue(q_l, normal_vector)
    lambda_r = get_max_eigenvalue(q_r, normal_vector)
    
    # Global max wave speed at the interface
    alpha = np.maximum(lambda_l, lambda_r)
    
    # Broadcasting Safety: ensure alpha has shape (1, N) if q is (4, N)
    # This prevents errors if dimensions are flipped
    if alpha.ndim == 1:
        alpha_b = alpha[None, :] 
    else:
        alpha_b = alpha

    # Numerical Flux Formula
    # F_num = 0.5 * (F_L + F_R) - 0.5 * alpha * (Q_R - Q_L)
    flux = 0.5 * (flux_n_l + flux_n_r) - 0.5 * alpha_b * (q_r - q_l)
    
    return flux

def primitive_to_conservative(prim):
    """
    Converts primitive variables [rho, u, v, p] to conservative variables [rho, rho*u, rho*v, E].

    Args:
        prim (np.array): Primitive variables.

    Returns:
        np.array: Conservative variables.
    """
    rho, u, v, p = prim
    rho_u = rho * u
    rho_v = rho * v
    kinetic_energy = 0.5 * rho * (u**2 + v**2)
    E = p / (gamma - 1) + kinetic_energy
    return np.array([rho, rho_u, rho_v, E])

def conservative_to_primitive(Q):
    """
    Converts conservative variables [rho, rho*u, rho*v, E] back to primitive [rho, u, v, p].
    Used for visualization and post-processing.

    Args:
        Q (np.array): Conservative variables.

    Returns:
        np.array: Primitive variables.
    """
    rho, rho_u, rho_v, E = Q
    
    # Avoid division by zero
    rho_safe = np.maximum(rho, 1e-12)
    
    u = rho_u / rho_safe
    v = rho_v / rho_safe
    p_val = pressure(Q) # Uses the existing pressure function
    
    return np.array([rho, u, v, p_val])