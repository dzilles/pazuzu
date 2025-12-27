import numpy as np

def vortex(x, y, t=0.0, **kwargs):
    """
    Computes the exact solution for the 2D Isentropic Vortex problem at time t.
    """
    gamma = kwargs.get("gamma", 1.4)
    beta = kwargs.get("beta", 5.0)  # Vortex strength
    
    u_inf = kwargs.get("u_inf", 1.0)
    v_inf = kwargs.get("v_inf", 1.0)
    
    # Center of vortex at t=0
    xc = kwargs.get("xc", 5.0)
    yc = kwargs.get("yc", 5.0)
    
    # Advect coordinates back to t=0 to find the vortex position in the moving frame
    x0 = x - u_inf * t
    y0 = y - v_inf * t
    
    r_sq = (x0 - xc)**2 + (y0 - yc)**2
    
    du = -(beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * (y0 - yc)
    dv = (beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * (x0 - xc)
    
    u = u_inf + du
    v = v_inf + dv
    
    T_inf = kwargs.get("T_inf", 1.0)
    T = T_inf - ((gamma - 1) * beta**2 / (8 * gamma * np.pi**2)) * np.exp(1 - r_sq)
    
    rho = T**(1.0 / (gamma - 1))
    p = rho**gamma
    
    return rho, u, v, p

def uniform(x, y, t=0.0, **kwargs):
    """
    Computes a uniform flow field state.
    """
    rho = kwargs.get("rho", 1.0)
    u = kwargs.get("u", 1.0)
    v = kwargs.get("v", 0.0)
    p = kwargs.get("p", 1.0)
    return rho, u, v, p

def acoustic_pulse(x, y, t=0.0, **kwargs):
    """
    Gaussian pressure pulse in a fluid at rest.
    """
    gamma = kwargs.get("gamma", 1.4)
    rho_inf = kwargs.get("rho_inf", 1.0)
    p_inf = kwargs.get("p_inf", 1.0)
    epsilon = kwargs.get("epsilon", 0.2)
    sigma = kwargs.get("sigma", 0.1)
    
    # Center of pulse
    xc = kwargs.get("xc", 0.0)
    yc = kwargs.get("yc", 0.0)
    
    r_sq = (x - xc)**2 + (y - yc)**2
    
    # Pressure perturbation
    p = p_inf + epsilon * np.exp(-r_sq / (2.0 * sigma**2))
    
    # Isentropic density: rho = rho_inf * (p/p_inf)^(1/gamma)
    rho = rho_inf * np.power(p / p_inf, 1.0 / gamma)
    
    u = kwargs.get("u_inf", 0.0)
    v = kwargs.get("v_inf", 0.0)
    
    return rho, u, v, p

def rest(x, y, t=0.0, **kwargs):
    """
    Computes a state at rest.
    """
    rho = kwargs.get("rho", 1.0)
    p = kwargs.get("p", 1.0)
    return rho, 0.0, 0.0, p

def sod_shock_tube(x, y, t=0.0, **kwargs):
    """
    Computes the initial state for the Sod Shock Tube problem.
    """
    # X-discontinuity location
    x0 = kwargs.get("x0", 0.5)
    
    if x < x0:
        # Left State
        rho = kwargs.get("rho_l", 1.0)
        u = kwargs.get("u_l", 0.0)
        v = kwargs.get("v_l", 0.0)
        p = kwargs.get("p_l", 1.0)
    else:
        # Right State
        rho = kwargs.get("rho_r", 0.125)
        u = kwargs.get("u_r", 0.0)
        v = kwargs.get("v_r", 0.0)
        p = kwargs.get("p_r", 0.125)
        
    return rho, u, v, p

def double_mach_reflection(x, y, t=0.0, **kwargs):

    """

    Woodward & Colella (1984) Double Mach Reflection setup.

    """

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



# Initial Condition Registry

_IC_REGISTRY = {

    "vortex": vortex,

    "uniform": uniform,

    "rest": rest,

    "sod_shock_tube": sod_shock_tube,

    "double_mach_reflection": double_mach_reflection,

    "acoustic_pulse": acoustic_pulse

}



def get_ic_function(name: str):

    """

    Returns the initial condition function corresponding to the given name.

    

    Args:

        name (str): The name of the initial condition.

        

    Returns:

        callable: The initial condition function.

        

    Raises:

        ValueError: If the initial condition name is not recognized.

    """

    if name not in _IC_REGISTRY:

        raise ValueError(

            f"Unknown initial condition: '{name}'. "

            f"Available ICs: {list(_IC_REGISTRY.keys())}"

        )

    return _IC_REGISTRY[name]



def register_ic(name: str, func: callable):

    """Registers a new initial condition function."""

    _IC_REGISTRY[name] = func
