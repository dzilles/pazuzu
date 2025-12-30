import warp as wp
from typing import Any
import src.kernels.boundary_conditions as bc
from src.kernels.structs import EquationParams32, EquationParams64

@wp.kernel
def compute_max_wave_speed(
    q: wp.array(dtype=Any, ndim=2),
    active_indices: wp.array(dtype=int),
    num_active: int,
    block_levels: wp.array(dtype=int),
    root_bounds: Any, # Also Any for vec4f/vec4d
    params: Any,
    max_inv_dt: wp.array(dtype=Any) # Output: size 1
):
    tid = wp.tid()
    Np = q.shape[1]
    
    # Grid Stride Loop
    count = num_active * Np
    if tid >= count:
        return

    # Decompose Index
    block_offset = tid // Np
    node_local = tid % Np
    
    pool_idx = active_indices[block_offset]
    level = block_levels[pool_idx]
    
    # Calculate dx, dy
    grid_dim_int = 1 << level
    
    domain_w = root_bounds[2] - root_bounds[0]
    domain_h = root_bounds[3] - root_bounds[1]
    
    template = domain_w
    f_grid_dim = bc.get_any_generic(template, grid_dim_int)

    dx = domain_w / f_grid_dim
    dy = domain_h / f_grid_dim
    min_spacing = wp.min(dx, dy)
    
    # Load State
    q_val = q[pool_idx, node_local]
    rho = q_val[0]
    rhou = q_val[1]
    rhov = q_val[2]
    E = q_val[3]
    
    # Generic Constants
    one = bc.get_one_generic(rho)
    half = bc.get_half_generic(rho)
    
    # Primitives
    # Guard against vacuum
    rho = wp.max(rho, params.rho_floor)
    inv_rho = one / rho
    u = rhou * inv_rho
    v = rhov * inv_rho
    
    p = (params.gamma - one) * (E - half * rho * (u*u + v*v))
    p = wp.max(p, params.p_floor)
    
    # Sound Speed
    c = wp.sqrt(params.gamma * p * inv_rho)
    
    # Velocity Magnitude
    vel = wp.sqrt(u*u + v*v)
    
    # Max Wave Speed
    lambda_val = vel + c
    
    # Inverse DT contribution
    inv_dt = lambda_val / min_spacing
    
    # Atomic Max
    wp.atomic_max(max_inv_dt, 0, inv_dt)
