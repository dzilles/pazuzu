import warp as wp
from src.kernels import boundary_conditions as bc
from typing import Any

@wp.kernel
def init_isentropic_vortex(
    x: wp.array(dtype=Any, ndim=2),
    y: wp.array(dtype=Any, ndim=2),
    q: wp.array(dtype=Any, ndim=2),
    active_indices: wp.array(dtype=int),
    num_active: int,
    params: Any,
    t: Any,
    beta: Any,
    radius: Any,
    center_x: Any,
    center_y: Any
):
    # Launch dimensions: (num_active, Np)
    block_idx, node_idx = wp.tid()
    
    pool_idx = active_indices[block_idx]
    
    xx = x[pool_idx, node_idx]
    yy = y[pool_idx, node_idx]
    
    # Vortex Parameters
    # Advecting with u_inf, v_inf
    x0 = center_x + params.u_inf * t 
    y0 = center_y + params.v_inf * t
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
    
    template = beta
    one = bc.get_one_generic(template)
    half = bc.get_half_generic(template)
    two = one + one
    pi = bc.get_any_generic(template, 3.141592653589793)

    S_2pi = beta / (two * pi)
    exp_term = wp.exp(half * (one - r2_scaled))
    
    # Scale perturbation by radius to keep beta as peak velocity
    du = -S_2pi * (dy / radius) * exp_term
    dv =  S_2pi * (dx / radius) * exp_term
    
    u = params.u_inf + du
    v = params.v_inf + dv
    
    # Correct Isentropic Relation: T = 1 - ((gamma-1)/gamma) * (S^2/8pi^2) * exp(...)
    T_sub = (gamma - one) / gamma * half * (S_2pi * S_2pi) * wp.exp(one - r2_scaled)
    T = one - T_sub
    
    rho = wp.pow(T, one / (gamma - one))
    p = wp.pow(rho, gamma)
    
    E = p / (gamma - one) + half * rho * (u*u + v*v)
    
    q[pool_idx, node_idx] = bc.make_vec4_generic(rho, rho*u, rho*v, E)
