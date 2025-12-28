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
    t: float
):
    # Launch dimensions: (num_active, Np)
    block_idx, node_idx = wp.tid()
    
    pool_idx = active_indices[block_idx]
    
    xx = x[pool_idx, node_idx]
    yy = y[pool_idx, node_idx]
    
    # Vortex Parameters
    x0 = 0.0 + t # Advecting with u=1
    y0 = 0.0
    beta = 5.0
    gamma = params.gamma
    
    dx = xx - x0
    dy = yy - y0
    r2 = dx*dx + dy*dy
    
    # f(x,y) = (beta / 2pi) * exp(0.5 * (1 - r^2))
    # But standard is: delta T = - (gamma-1)/(2gamma) * beta^2 * exp(1-r^2)
    # let's use standard:
    
    f = (1.0 - r2)
    
    # Periodic BC handling (heuristic)? 
    # For now assume domain is large enough or conformant.
    
    # Temperature and Density
    # T = 1 - (gamma-1)/2 * M_vortex^2 * exp(1-r^2)
    # Using beta as strength.
    
    S = 13.5 # Strength
    # Common test case:
    # u = 1 - S * y * exp(0.5*(1-r^2))
    # v = 0 + S * x * exp(0.5*(1-r^2)) (Assuming center 0,0)
    # T = 1 - (gamma-1)*S^2/(8*pi^2) * exp(1-r^2) ? 
    
    # Let's use the Gatton/Project standard if known.
    # Standard Isentropic Vortex:
    # u_inf = 1, v_inf = 0.
    # du = - (S / 2pi) * dy * exp(0.5*(1-r2))
    # dv =   (S / 2pi) * dx * exp(0.5*(1-r2))
    # T = 1 - (gamma-1) * (S^2 / (8*pi*pi)) * exp(1-r2)
    # rho = T^(1/(gamma-1))
    # p = rho^gamma
    
    S_2pi = 5.0 / (2.0 * 3.14159265359)
    exp_term = wp.exp(0.5 * (1.0 - r2))
    
    du = -S_2pi * dy * exp_term
    dv =  S_2pi * dx * exp_term
    
    u = 1.0 + du
    v = 0.0 + dv
    
    T_sub = (gamma - 1.0) * 0.5 * (S_2pi * S_2pi) * wp.exp(1.0 - r2)
    T = 1.0 - T_sub
    
    rho = wp.pow(T, 1.0 / (gamma - 1.0))
    p = wp.pow(rho, gamma)
    
    E = p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v)
    
    q[pool_idx, node_idx] = wp.vec4(rho, rho*u, rho*v, E)
