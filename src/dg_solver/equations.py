import numpy as np

# Gas constant
gamma = 1.4

def pressure(Q):
    rho, rho_u, rho_v, E = Q
    # Soft clip ist besser als hard clip, aber für den Anfang okay
    rho = np.maximum(rho, 1e-12)
    
    # Optimierung: 0.5 * (rho_u^2 + rho_v^2) / rho ist numerisch stabiler 
    # als 0.5 * rho * (u^2 + v^2), da u=rho_u/rho -> u^2 = rho_u^2/rho^2
    # Das spart eine Division und ist bei kleinem rho genauer.
    kinetic_energy = 0.5 * (rho_u**2 + rho_v**2) / rho
    p_val = (gamma - 1) * (E - kinetic_energy)
    
    return np.maximum(p_val, 1e-12)

def euler_fluxes(Q):
    """
    Computes both F and G at once to save unpack operations/pressure calcs.
    Returns F, G
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
    f3 = rho_u * v     # rho_v * u ist das gleiche wie rho_u * v
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
    rho = Q[0]
    # Vermeide Division durch fast-Null
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
    Computes the Lax-Friedrichs numerical flux projected on the normal.
    """
    nx, ny = normal_vector
    
    # Compute fluxes
    # Wir nutzen hier die kombinierte Funktion, um Rechenzeit zu sparen
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
    
    # Broadcasting Safety: alpha auf (1, N) bringen falls q (4, N) ist
    # Das verhindert Fehler, falls die Dimensionen mal drehen
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
    Used for visualization.
    """
    rho, rho_u, rho_v, E = Q
    
    # Avoid division by zero
    rho_safe = np.maximum(rho, 1e-12)
    
    u = rho_u / rho_safe
    v = rho_v / rho_safe
    p_val = pressure(Q) # Nutzt die bereits vorhandene pressure-Funktion
    
    return np.array([rho, u, v, p_val])