import warp as wp
from src.kernels.utils import make_vec4_generic
from typing import Any

@wp.func
def pressure(q: Any, params: Any):
    """
    Calculates the pressure from conservative variables with safety checks.
    """
    # Type promotion helper
    zero = q[0] - q[0]
    
    rho = wp.max(q[0], zero + params.rho_floor)
    rho_u = q[1]
    rho_v = q[2]
    E = q[3]
    
    kin_energy = params.half * (rho_u*rho_u + rho_v*rho_v) / rho
    p = (params.gamma - params.one) * (E - kin_energy)
    return wp.max(p, zero + params.p_floor)

@wp.func
def flux_x(q: Any, params: Any):
    """
    Computes the Euler Flux function F(q) in the x-direction.
    """
    zero = q[0] - q[0]
    rho = wp.max(q[0], zero + params.rho_floor)
    p = pressure(q, params)
    u = q[1] / rho
    v = q[2] / rho
    
    # Ensure energy is consistent with floored pressure
    E_consistent = p / (params.gamma - params.one) + params.half * rho * (u*u + v*v)
    
    return make_vec4_generic(q[1], q[1]*u + p, q[2]*u, (E_consistent + p)*u)

@wp.func
def flux_y(q: Any, params: Any):
    """
    Computes the Euler Flux function G(q) in the y-direction.
    """
    zero = q[0] - q[0]
    rho = wp.max(q[0], zero + params.rho_floor)
    p = pressure(q, params)
    u = q[1] / rho
    v = q[2] / rho
    
    # Ensure energy is consistent with floored pressure
    E_consistent = p / (params.gamma - params.one) + params.half * rho * (u*u + v*v)
    
    return make_vec4_generic(q[2], q[1]*v, q[2]*v + p, (E_consistent + p)*v)