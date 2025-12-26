import warp as wp
from src.kernels import boundary_conditions as bc
from src.physics.models.euler_laws import pressure
from typing import Any

@wp.func
def temperature(q: Any, params: Any):
    """
    Calculates temperature T = p / (rho * R)
    Assuming R = Cp - Cv and gamma = Cp/Cv -> R = Cp * (1 - 1/gamma)
    Actually, we can use p / (rho * R) if we have R.
    Using the relation: p = rho * R * T  => T = p / (rho * R)
    R = Cp * (gamma - 1) / gamma
    """
    zero = q[0] - q[0]
    rho = wp.max(q[0], zero + params.rho_floor)
    p = pressure(q, params)
    
    return p / (rho * params.gas_constant)

@wp.func
def viscous_flux_x(q: Any, grad_u: Any, grad_v: Any, grad_T: Any, params: Any):
    """
    Computes the viscous flux in the x-direction.
    Fv = [0, tau_xx, tau_xy, u*tau_xx + v*tau_xy - qx]
    """
    # Dynamic viscosity mu
    mu = params.mu
    
    # Divergence of velocity (2D)
    div_u = grad_u[0] + grad_v[1]
    
    # Stress tensor components (Stokes' hypothesis: lambda = -2/3 * mu)
    # tau_xx = 2 * mu * (du/dx - 1/3 * div(u))
    # tau_xy = mu * (du/dy + dv/dx)
    
    two = params.one + params.one
    three = two + params.one
    one_third = params.one / three
    
    tau_xx = two * mu * (grad_u[0] - one_third * div_u)
    tau_xy = mu * (grad_u[1] + grad_v[0])
    
    # Heat flux qx = -kappa * dT/dx
    # kappa = mu * Cp / Pr
    kappa = mu * params.cp / params.prandtl
    qx = -kappa * grad_T[0]
    
    # Velocity components
    zero = q[0] - q[0]
    rho = wp.max(q[0], zero + params.rho_floor)
    u = q[1] / rho
    v = q[2] / rho
    
    return bc.make_vec4_generic(
        zero,
        tau_xx,
        tau_xy,
        u * tau_xx + v * tau_xy - qx
    )

@wp.func
def viscous_flux_y(q: Any, grad_u: Any, grad_v: Any, grad_T: Any, params: Any):
    """
    Computes the viscous flux in the y-direction.
    Gv = [0, tau_yx, tau_yy, u*tau_yx + v*tau_yy - qy]
    """
    mu = params.mu
    div_u = grad_u[0] + grad_v[1]
    
    two = params.one + params.one
    three = two + params.one
    one_third = params.one / three
    
    tau_yy = two * mu * (grad_v[1] - one_third * div_u)
    tau_yx = mu * (grad_u[1] + grad_v[0]) # tau_yx = tau_xy
    
    kappa = mu * params.cp / params.prandtl
    qy = -kappa * grad_T[1]
    
    zero = q[0] - q[0]
    rho = wp.max(q[0], zero + params.rho_floor)
    u = q[1] / rho
    v = q[2] / rho
    
    return bc.make_vec4_generic(
        zero,
        tau_yx,
        tau_yy,
        u * tau_yx + v * tau_yy - qy
    )
