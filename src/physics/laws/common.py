import numpy as np

# Gas constant
gamma = 1.4

def pressure(q):
    """
    Computes the pressure from the conservative state vector using the ideal gas law.
    
    Equation of state: p = (gamma - 1) * (E - 0.5 * rho * (u^2 + v^2))

    Args:
        q (np.array): Conservative variables [rho, rho*u, rho*v, E].

    Returns:
        np.array: The pressure.
    """
    rho, rho_u, rho_v, E = q
    # Soft clip to avoid division by zero (better than hard clip for optimization, but fine here)
    rho = np.maximum(rho, 1e-12)
    
    # Optimization: 0.5 * (rho_u^2 + rho_v^2) / rho is numerically more stable 
    # than 0.5 * rho * (u^2 + v^2), as u=rho_u/rho -> u^2 = rho_u^2/rho^2
    # This saves one division and preserves accuracy for small rho.
    kinetic_energy = 0.5 * (rho_u**2 + rho_v**2) / rho
    p_val = (gamma - 1) * (E - kinetic_energy)
    
    return np.maximum(p_val, 1e-12)

def euler_fluxes(q):
    """
    Computes both Euler Flux vectors F(q) and G(q) simultaneously.
    
    F(q) = [rho*u, rho*u^2 + p, rho*u*v, (E+p)*u]
    G(q) = [rho*v, rho*v*u, rho*v^2 + p, (E+p)*v]

    Args:
        q (np.array): Conservative variables.

    Returns:
        tuple: (F, G) numpy arrays containing the fluxes in x and y directions.
    """
    rho, rho_u, rho_v, E = q
    rho = np.maximum(rho, 1e-12)
    p = pressure(q)
    
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

def get_max_eigenvalue(q, normal_vector):
    """
    Computes the maximum eigenvalue (wave speed) of the system in a given normal direction.
    
    lambda_max = |u_n| + c

    Args:
        q (np.array): Conservative variables.
        normal_vector (tuple): (nx, ny) components of the normal vector.

    Returns:
        np.array: The maximum wave speed.
    """
    rho = q[0]
    # Avoid division by near-zero
    rho = np.maximum(rho, 1e-12)
    p = pressure(q)
    
    # Sound speed
    c = np.sqrt(gamma * p / rho)
    
    nx, ny = normal_vector
    # Velocity projected on normal
    # u_n = (rho_u * nx + rho_v * ny) / rho
    u_n = (q[1] * nx + q[2] * ny) / rho
    
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

def hllc_flux(q_l, q_r, normal_vector):
    """
    Computes the HLLC (Harten-Lax-van Leer-Contact) numerical flux.
    Resolves contact discontinuities better than Rusanov.

    Args:
        q_l (np.array): State on the left (4, ...).
        q_r (np.array): State on the right (4, ...).
        normal_vector (tuple): Normal vector (nx, ny).

    Returns:
        np.array: The HLLC flux.
    """
    nx, ny = normal_vector
    
    # 1. Primitives
    rho_l = np.maximum(q_l[0], 1e-12)
    u_l = q_l[1] / rho_l
    v_l = q_l[2] / rho_l
    p_l = pressure(q_l)
    
    rho_r = np.maximum(q_r[0], 1e-12)
    u_r = q_r[1] / rho_r
    v_r = q_r[2] / rho_r
    p_r = pressure(q_r)

    # Normal velocities
    un_l = u_l * nx + v_l * ny
    un_r = u_r * nx + v_r * ny

    # Sound speeds
    c_l = np.sqrt(gamma * p_l / rho_l)
    c_r = np.sqrt(gamma * p_r / rho_r)

    # 2. Wave Speed Estimates
    s_l = np.minimum(un_l - c_l, un_r - c_r)
    s_r = np.maximum(un_l + c_l, un_r + c_r)

    # 3. Compute Flux
    # We use numpy masking to handle the different regions
    
    # Contact Wave Speed S_star
    denom = rho_l * (s_l - un_l) - rho_r * (s_r - un_r)
    denom = np.where(np.abs(denom) < 1e-10, 1e-10, denom) # Avoid div/0
    s_star = (p_r - p_l + rho_l * un_l * (s_l - un_l) - rho_r * un_r * (s_r - un_r)) / denom

    # Calculate standard fluxes
    F_l, G_l = euler_fluxes(q_l)
    F_r, G_r = euler_fluxes(q_r)
    flux_l = F_l * nx + G_l * ny
    flux_r = F_r * nx + G_r * ny

    # Create output array
    flux = np.zeros_like(q_l)

    # Masks for regions
    # Region 1: S_L >= 0 (Supersonic Right)
    mask_l = s_l >= 0
    # Region 4: S_R <= 0 (Supersonic Left)
    mask_r = s_r <= 0
    # Star Region
    mask_star = ~(mask_l | mask_r)
    
    # Apply simple fluxes
    # Note: If q_l is 1D (shape (4,)), masks are scalars. If (4, N), masks are (N,)
    if mask_l.ndim == 0:
        if mask_l: return flux_l
        if mask_r: return flux_r
        # Star logic for scalar input
        if s_star >= 0:
            # Star Left
            D_l = (s_l - un_l) / (s_l - s_star)
            factor_l = rho_l * D_l
            u_star_l = u_l + (s_star - un_l) * nx
            v_star_l = v_l + (s_star - un_l) * ny
            E_star_l = (q_l[3]/rho_l) + (s_star - un_l) * (s_star + p_l / (rho_l * (s_l - un_l)))
            
            U_star_l = np.array([factor_l, factor_l*u_star_l, factor_l*v_star_l, factor_l*E_star_l])
            return flux_l + s_l * (U_star_l - q_l)
        else:
            # Star Right
            D_r = (s_r - un_r) / (s_r - s_star)
            factor_r = rho_r * D_r
            u_star_r = u_r + (s_star - un_r) * nx
            v_star_r = v_r + (s_star - un_r) * ny
            E_star_r = (q_r[3]/rho_r) + (s_star - un_r) * (s_star + p_r / (rho_r * (s_r - un_r)))
            
            U_star_r = np.array([factor_r, factor_r*u_star_r, factor_r*v_star_r, factor_r*E_star_r])
            return flux_r + s_r * (U_star_r - q_r)
    else:
        # Vectorized logic
        flux[:, mask_l] = flux_l[:, mask_l]
        flux[:, mask_r] = flux_r[:, mask_r]
        
        # Star Regions
        if np.any(mask_star):
            mask_sl = mask_star & (s_star >= 0)
            mask_sr = mask_star & (s_star < 0)
            
            # Star Left
            if np.any(mask_sl):
                sl_valid = s_l[mask_sl]
                unl_valid = un_l[mask_sl]
                ss_valid = s_star[mask_sl]
                
                D_l = (sl_valid - unl_valid) / (sl_valid - ss_valid)
                rho_l_v = rho_l[mask_sl]
                factor_l = rho_l_v * D_l
                
                u_star_l = u_l[mask_sl] + (ss_valid - unl_valid) * nx
                v_star_l = v_l[mask_sl] + (ss_valid - unl_valid) * ny
                
                E_term = (q_l[3, mask_sl]/rho_l_v) + (ss_valid - unl_valid) * (ss_valid + p_l[mask_sl] / (rho_l_v * (sl_valid - unl_valid)))
                
                U_star_l = np.stack([
                    factor_l,
                    factor_l * u_star_l,
                    factor_l * v_star_l,
                    factor_l * E_term
                ])
                
                flux[:, mask_sl] = flux_l[:, mask_sl] + sl_valid * (U_star_l - q_l[:, mask_sl])

            # Star Right
            if np.any(mask_sr):
                sr_valid = s_r[mask_sr]
                unr_valid = un_r[mask_sr]
                ss_valid = s_star[mask_sr]
                
                D_r = (sr_valid - unr_valid) / (sr_valid - ss_valid)
                rho_r_v = rho_r[mask_sr]
                factor_r = rho_r_v * D_r
                
                u_star_r = u_r[mask_sr] + (ss_valid - unr_valid) * nx
                v_star_r = v_r[mask_sr] + (ss_valid - unr_valid) * ny
                
                E_term = (q_r[3, mask_sr]/rho_r_v) + (ss_valid - unr_valid) * (ss_valid + p_r[mask_sr] / (rho_r_v * (sr_valid - unr_valid)))
                
                U_star_r = np.stack([
                    factor_r,
                    factor_r * u_star_r,
                    factor_r * v_star_r,
                    factor_r * E_term
                ])
                
                flux[:, mask_sr] = flux_r[:, mask_sr] + sr_valid * (U_star_r - q_r[:, mask_sr])

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

def conservative_to_primitive(q):
    """
    Converts conservative variables [rho, rho*u, rho*v, E] back to primitive [rho, u, v, p].
    Used for visualization and post-processing.

    Args:
        q (np.array): Conservative variables.

    Returns:
        np.array: Primitive variables.
    """
    rho, rho_u, rho_v, E = q
    
    # Avoid division by zero
    rho_safe = np.maximum(rho, 1e-12)
    
    u = rho_u / rho_safe
    v = rho_v / rho_safe
    p_val = pressure(q) # Uses the existing pressure function
    
    return np.array([rho, u, v, p_val])