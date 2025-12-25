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

def acoustic_pulse(x, y, t=0.0):
    """
    Gaussian pressure pulse in a fluid at rest.
    Used to verify non-reflecting (characteristic) boundary conditions.
    """
    gamma = 1.4
    rho_inf = 1.0
    p_inf = 1.0
    epsilon = 0.2
    sigma = 0.1
    
    r_sq = x*x + y*y
    
    # Pressure perturbation
    p = p_inf + epsilon * np.exp(-r_sq / (2.0 * sigma**2))
    
    # Isentropic density: rho = rho_inf * (p/p_inf)^(1/gamma)
    rho = rho_inf * np.power(p / p_inf, 1.0 / gamma)
    
    u = 0.0
    v = 0.0
    
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

def sod_shock_tube(x, y, t=0.0):
    """
    Computes the initial state for the Sod Shock Tube problem.

    Args:
        x (float): Physical x-coordinate.
        y (float): Physical y-coordinate.
        t (float, optional): Simulation time. Defaults to 0.0.

    Returns:
        tuple: (rho, u, v, p)
    """
    if x < 0.5:
        # Left State
        return 1.0, 0.0, 0.0, 1.0
    else:
        # Right State
        return 0.125, 0.0, 0.0, 0.125

def double_mach_reflection(x, y, t=0.0):
    """
    Woodward & Colella (1984) Double Mach Reflection setup.
    A Mach 10 shock hits a 30-degree wedge (rotated so wedge is on x-axis).
    """
    # Angle of the shock with the wall (x-axis) is 60 degrees.
    # tan(60) = sqrt(3)
    # The shock line at t=0: x = 1/6 + y / sqrt(3)
    # Shock speed along x-axis: V_x = 10 / sin(60) = 20 / sqrt(3)
    
    sin_60 = np.sqrt(3.0) / 2.0
    tan_60 = np.sqrt(3.0)
    
    shock_x = 1.0/6.0 + y / tan_60 + (10.0 / sin_60) * t
    
    if x < shock_x:
        # Post-shock (Left) State
        rho = 8.0
        u = 8.25 * np.cos(np.radians(30))
        v = -8.25 * np.sin(np.radians(30))
        p = 116.5
    else:
        # Pre-shock (Right) State
        rho = 1.4
        u = 0.0
        v = 0.0
        p = 1.0
        
    return rho, u, v, p