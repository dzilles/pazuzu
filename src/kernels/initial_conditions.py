import warp as wp
from src.kernels.structs import EquationParams32

@wp.kernel
def init_isentropic_vortex(
    x: wp.array(dtype=float, ndim=2),
    y: wp.array(dtype=float, ndim=2),
    q: wp.array(dtype=wp.vec4, ndim=2),
    active_indices: wp.array(dtype=int),
    num_active: int,
    params: EquationParams32,
    t: float,
    beta: float,
    radius: float
):
    # Launch dimensions: (num_active, Np)
    block_idx, node_idx = wp.tid()
    
    pool_idx = active_indices[block_idx]
    
    xx = x[pool_idx, node_idx]
    yy = y[pool_idx, node_idx]
    
    # Vortex Parameters
    # Advecting with u_inf, v_inf
    x0 = params.u_inf * t 
    y0 = params.v_inf * t
    gamma = params.gamma
    
    dx = xx - x0
    dy = yy - y0
    r2 = dx*dx + dy*dy
    r2_scaled = r2 / (radius * radius)
    
    # f(x,y) = (beta / 2pi) * exp(0.5 * (1 - r^2))
    
    # Standard Isentropic Vortex:
    # u_inf = 1, v_inf = 0.
    # du = - (S / 2pi) * dy * exp(0.5*(1-r2))
    # dv =   (S / 2pi) * dx * exp(0.5*(1-r2))
    # T = 1 - (gamma-1) * (S^2 / (8*pi*pi)) * exp(1-r2)
    # rho = T^(1/(gamma-1))
    # p = rho^gamma
    
    S_2pi = beta / (2.0 * 3.14159265359)
    exp_term = wp.exp(0.5 * (1.0 - r2_scaled))
    
    du = -S_2pi * dy * exp_term
    dv =  S_2pi * dx * exp_term
    
    u = params.u_inf + du
    v = params.v_inf + dv
    
    T_sub = (gamma - 1.0) * 0.5 * (S_2pi * S_2pi) * wp.exp(1.0 - r2_scaled)
    T = 1.0 - T_sub
    
    rho = wp.pow(T, 1.0 / (gamma - 1.0))
    p = wp.pow(rho, gamma)
    
    E = p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v)
    
    q[pool_idx, node_idx] = wp.vec4(rho, rho*u, rho*v, E)
