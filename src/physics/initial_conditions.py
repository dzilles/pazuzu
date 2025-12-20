import numpy as np

def vortex(x, y, t=0.0):
    """
    2D Isentropic Vortex problem.
    The solution is a vortex that advects diagonally across the domain.
    Exact solution: q(x, y, t) = q_0(x - u_inf*t, y - v_inf*t)
    """
    gamma = 1.4
    beta = 5.0  # Vortex strength
    
    u_inf = 1.0
    v_inf = 1.0
    
    # Advect coordinates back to t=0
    x0 = x - u_inf * t
    y0 = y - v_inf * t
    
    # Center of vortex at t=0
    xc, yc = 5.0, 5.0
    
    r_sq = (x0 - xc)**2 + (y0 - yc)**2
    
    du = -(beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * (y0 - yc)
    dv = (beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * (x0 - xc)
    
    u = u_inf + du
    v = v_inf + dv
    
    T_inf = 1.0
    T = T_inf - ((gamma - 1) * beta**2 / (8 * gamma * np.pi**2)) * np.exp(1 - r_sq)
    
    rho = T**(1.0 / (gamma - 1))
    p = rho**gamma
    
    return rho, u, v, p

def uniform(x, y, t=0.0):
    """
    Uniform flow.
    """
    rho = 1.0
    u = 1.0
    v = 0.0
    p = 1.0
    return rho, u, v, p
