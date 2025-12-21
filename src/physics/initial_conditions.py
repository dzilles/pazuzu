import numpy as np

def vortex(x, y, t=0.0):
    """
    Computes the exact solution for the 2D Isentropic Vortex problem at time t.
    
    The vortex is a common test case for Euler solvers. It consists of a mean flow
    (u_inf, v_inf) with a superimposed perturbation that satisfies the Euler equations.
    The vortex advects with the mean flow without changing shape.

    Args:
        x (float): Physical x-coordinate.
        y (float): Physical y-coordinate.
        t (float, optional): Simulation time. Defaults to 0.0.

    Returns:
        tuple: (rho, u, v, p) - Density, x-velocity, y-velocity, pressure.
    """
    gamma = 1.4
    beta = 5.0  # Vortex strength
    
    u_inf = 1.0
    v_inf = 1.0
    
    # Advect coordinates back to t=0 to find the vortex position in the moving frame
    # Exact solution property: q(x, y, t) = q_0(x - u_inf*t, y - v_inf*t)
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
    Computes a uniform flow field state.

    Args:
        x (float): Physical x-coordinate.
        y (float): Physical y-coordinate.
        t (float, optional): Simulation time. Defaults to 0.0.

    Returns:
        tuple: (rho, u, v, p) with rho=1.0, u=1.0, v=0.0, p=1.0.
    """
    rho = 1.0
    u = 1.0
    v = 0.0
    p = 1.0
    return rho, u, v, p

def rest(x, y, t=0.0):
    """
    Computes a state at rest.

    Args:
        x (float): Physical x-coordinate.
        y (float): Physical y-coordinate.
        t (float, optional): Simulation time. Defaults to 0.0.

    Returns:
        tuple: (rho, u, v, p) with rho=1.0, u=0.0, v=0.0, p=1.0.
    """
    return 1.0, 0.0, 0.0, 1.0
