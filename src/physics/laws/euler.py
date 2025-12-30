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
    v = q[2] / rho
    
    # Ensure energy is consistent with floored pressure
    E_consistent = p / (params.gamma - params.one) + params.half * rho * (u*u + v*v)
    
    return bc.make_vec4_generic(q[1], q[1]*u + p, q[2]*u, (E_consistent + p)*u)

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
    
    return bc.make_vec4_generic(q[2], q[1]*v, q[2]*v + p, (E_consistent + p)*v)

@wp.func
def rusanov_flux(q_L: Any, q_R: Any, nx: Any, ny: Any, params: Any):
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
    return params.half * (Fn_L + Fn_R) - params.half * lambda_max * (q_R - q_L)

@wp.func
def hllc_flux(q_L: Any, q_R: Any, nx: Any, ny: Any, params: Any):
    zero = q_L[0] - q_L[0]
    one = params.one
    
    # Primitives L
    rho_L = wp.max(q_L[0], zero + params.rho_floor)
    u_L = q_L[1] / rho_L
    v_L = q_L[2] / rho_L
    p_L = pressure(q_L, params)
    c_L = wp.sqrt(params.gamma * p_L / rho_L)
    vn_L = u_L * nx + v_L * ny
    E_L = q_L[3]
    
    # Primitives R
    rho_R = wp.max(q_R[0], zero + params.rho_floor)
    u_R = q_R[1] / rho_R
    v_R = q_R[2] / rho_R
    p_R = pressure(q_R, params)
    c_R = wp.sqrt(params.gamma * p_R / rho_R)
    vn_R = u_R * nx + v_R * ny
    E_R = q_R[3]
    
    # Wave Speeds (Davis estimate)
    S_L = wp.min(vn_L - c_L, vn_R - c_R)
    S_R = wp.max(vn_L + c_L, vn_R + c_R)
    
    # If supersonic, return F_L or F_R
    if S_L >= zero:
        F_L = flux_x(q_L, params)
        G_L = flux_y(q_L, params)
        return F_L * nx + G_L * ny
    elif S_R <= zero:
        F_R = flux_x(q_R, params)
        G_R = flux_y(q_R, params)
        return F_R * nx + G_R * ny
        
    # Subsonic - Compute S_star
    # Denominator check
    denom = rho_L * (S_L - vn_L) - rho_R * (S_R - vn_R)
    
    # Check for singularity if fallback is enabled
    if params.hllc_fallback == 1 and wp.abs(denom) < 1e-10:
        # Print warning (Warp kernel print)
        print("Warning: HLLC singularity detected. Switching to Rusanov.")
        return rusanov_flux(q_L, q_R, nx, ny, params)

    rho_term = (p_R - p_L + rho_L * vn_L * (S_L - vn_L) - rho_R * vn_R * (S_R - vn_R))
    S_star = rho_term / denom 
    
    # HLLC Flux
    if S_star >= zero:
        # Left Star
        factor_L = rho_L * (S_L - vn_L) / (S_L - S_star)
        
        # q*_L components constructed to satisfy Rankine-Hugoniot
        qs_rho_L = one
        qs_u_L = u_L + (S_star - vn_L) * nx
        qs_v_L = v_L + (S_star - vn_L) * ny
        qs_E_L = E_L/rho_L + (S_star - vn_L) * (S_star + p_L/(rho_L*(S_L - vn_L)))
        
        q_star_L = bc.make_vec4_generic(qs_rho_L, qs_u_L, qs_v_L, qs_E_L) * factor_L
        
        F_L = flux_x(q_L, params)
        G_L = flux_y(q_L, params)
        Fn_L = F_L * nx + G_L * ny
        
        return Fn_L + (q_star_L - q_L) * S_L

    else:
        # Right Star
        factor_R = rho_R * (S_R - vn_R) / (S_R - S_star)
        
        qs_rho_R = one
        qs_u_R = u_R + (S_star - vn_R) * nx
        qs_v_R = v_R + (S_star - vn_R) * ny
        qs_E_R = E_R/rho_R + (S_star - vn_R) * (S_star + p_R/(rho_R*(S_R - vn_R)))
        
        q_star_R = bc.make_vec4_generic(qs_rho_R, qs_u_R, qs_v_R, qs_E_R) * factor_R
        
        F_R = flux_x(q_R, params)
        G_R = flux_y(q_R, params)
        Fn_R = F_R * nx + G_R * ny
        
        return Fn_R + (q_star_R - q_R) * S_R