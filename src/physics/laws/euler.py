import warp as wp
from src.kernels import boundary_conditions as bc
from src.kernels.structs import EquationParams32
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
    
    return bc.make_vec4_generic(q[1], q[1]*u + p, q[2]*u, (q[3] + p)*u)

@wp.func
def flux_y(q: Any, params: Any):
    """
    Computes the Euler Flux function G(q) in the y-direction.
    """
    zero = q[0] - q[0]
    rho = wp.max(q[0], zero + params.rho_floor)
    p = pressure(q, params)
    v = q[2] / rho
    
    return bc.make_vec4_generic(q[2], q[1]*v, q[2]*v + p, (q[3] + p)*v)

@wp.func
def rusanov_flux(q_L: Any, q_R: Any, nx: float, ny: float, params: Any):
    """
    Computes the Rusanov (LLF) numerical flux at an interface with normal (nx, ny).
    """
    # Type promotion helper
    zero = q_L[0] - q_L[0]

    # Fluxes
    F_L = flux_x(q_L, params)
    G_L = flux_y(q_L, params)
    Fn_L = F_L * nx + G_L * ny
    
    F_R = flux_x(q_R, params)
    G_R = flux_y(q_R, params)
    Fn_R = F_R * nx + G_R * ny
    
    # Wave Speeds
    # L
    rho_L = wp.max(q_L[0], zero + params.rho_floor)
    p_L = pressure(q_L, params)
    c_L = wp.sqrt(params.gamma * p_L / rho_L)
    vn_L = (q_L[1] * nx + q_L[2] * ny) / rho_L
    max_ev_L = wp.abs(vn_L) + c_L
    
    # R
    rho_R = wp.max(q_R[0], zero + params.rho_floor)
    p_R = pressure(q_R, params)
    c_R = wp.sqrt(params.gamma * p_R / rho_R)
    vn_R = (q_R[1] * nx + q_R[2] * ny) / rho_R
    max_ev_R = wp.abs(vn_R) + c_R
    
    lambda_max = wp.max(max_ev_L, max_ev_R)
    
    # Rusanov Formula
    # F* = 0.5 * (Fn_L + Fn_R) - 0.5 * lambda * (q_R - q_L)
    return 0.5 * (Fn_L + Fn_R) - 0.5 * lambda_max * (q_R - q_L)
