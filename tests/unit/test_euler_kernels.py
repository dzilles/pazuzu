import numpy as np
import pytest
import warp as wp
import os
import sys
from typing import Any

# Add project root (parent of src) to path to allow 'src.' imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.kernels import euler_kernels
from src.kernels.structs import EquationParams32

@pytest.fixture(scope="module")
def device():
    wp.init()
    return "cpu"

@pytest.fixture(scope="module")
def params():
    p = EquationParams32()
    p.gamma = 1.4
    p.rho_floor = 1e-6
    p.p_floor = 1e-6
    p.half = 0.5
    p.one = 1.0
    return p

# Helper kernels to test wp.func
@wp.kernel
def kernel_pressure(
    q: wp.array(dtype=wp.vec4),
    p_out: wp.array(dtype=float),
    params: Any
):
    tid = wp.tid()
    p_out[tid] = euler_kernels.pressure(q[tid], params)

@wp.kernel
def kernel_flux_x(
    q: wp.array(dtype=wp.vec4),
    flux_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    flux_out[tid] = euler_kernels.flux_x(q[tid], params)

@wp.kernel
def kernel_flux_y(
    q: wp.array(dtype=wp.vec4),
    flux_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    flux_out[tid] = euler_kernels.flux_y(q[tid], params)

@wp.kernel
def kernel_rusanov_flux(
    q_l: wp.array(dtype=wp.vec4),
    q_r: wp.array(dtype=wp.vec4),
    nx: float,
    ny: float,
    flux_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    ql = q_l[tid]
    qr = q_r[tid]
    
    # Compute normal fluxes needed for the rusanov_flux signature
    fln = euler_kernels.flux_x(ql, params) * nx + euler_kernels.flux_y(ql, params) * ny
    frn = euler_kernels.flux_x(qr, params) * nx + euler_kernels.flux_y(qr, params) * ny
    
    flux_out[tid] = euler_kernels.rusanov_flux(ql, qr, fln, frn, nx, ny, params)

def test_pressure_calculation(device, params):
    # rho=1, u=1, v=0, p=1 => E = p/(gamma-1) + 0.5*rho*u^2 = 1/0.4 + 0.5 = 2.5 + 0.5 = 3.0
    q_val = [1.0, 1.0, 0.0, 3.0]
    q_wp = wp.array([q_val], dtype=wp.vec4, device=device)
    p_out_wp = wp.zeros(1, dtype=float, device=device)
    
    wp.launch(
        kernel=kernel_pressure,
        dim=1,
        inputs=[q_wp, p_out_wp, params],
        device=device
    )
    
    assert pytest.approx(p_out_wp.numpy()[0]) == 1.0

def test_flux_x_calculation(device, params):
    # rho=1, u=2, v=0, p=1
    # q = [1, 2, 0, E]
    # E = 1/0.4 + 0.5*1*2^2 = 2.5 + 2 = 4.5
    q_val = [1.0, 2.0, 0.0, 4.5]
    q_wp = wp.array([q_val], dtype=wp.vec4, device=device)
    flux_out_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    
    wp.launch(
        kernel=kernel_flux_x,
        dim=1,
        inputs=[q_wp, flux_out_wp, params],
        device=device
    )
    
    # Flux F = [rho*u, rho*u^2 + p, rho*u*v, (E+p)*u]
    # F = [2, 1*4 + 1, 0, (4.5+1)*2] = [2, 5, 0, 11]
    expected = np.array([2.0, 5.0, 0.0, 11.0])
    np.testing.assert_allclose(flux_out_wp.numpy()[0], expected, atol=1e-6)

def get_random_state(subsonic=True, gamma=1.4):
    """Generate a random physically valid state."""
    rho = np.random.uniform(0.5, 1.5)
    p = np.random.uniform(0.5, 1.5)
    c = np.sqrt(gamma * p / rho)
    
    if subsonic:
        u = np.random.uniform(-0.5 * c, 0.5 * c)
        v = np.random.uniform(-0.5 * c, 0.5 * c)
    else:
        u = np.random.uniform(1.5 * c, 2.5 * c)
        v = np.random.uniform(1.5 * c, 2.5 * c)
        
    E = p / (gamma - 1.0) + 0.5 * rho * (u**2 + v**2)
    return [rho, rho*u, rho*v, E]

@pytest.mark.parametrize("state_type", ["subsonic", "supersonic"])
@pytest.mark.parametrize("seed", [42, 123])
def test_flux_consistency(state_type, seed, device, params):
    """Verify that Rusanov flux with identical states equals physical flux: Frusanov(q,q,n) = Fx(q)*nx + Fy(q)*ny"""
    np.random.seed(seed)
    q_val = get_random_state(subsonic=(state_type == "subsonic"), gamma=params.gamma)
    
    # Random normal vector
    theta = np.random.uniform(0, 2 * np.pi)
    nx, ny = np.cos(theta), np.sin(theta)
    
    q_wp = wp.array([q_val], dtype=wp.vec4, device=device)
    flux_rusanov_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    flux_x_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    flux_y_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    
    # Compute Rusanov flux
    wp.launch(kernel=kernel_rusanov_flux, dim=1, inputs=[q_wp, q_wp, float(nx), float(ny), flux_rusanov_wp, params], device=device)
    
    # Compute physical fluxes
    wp.launch(kernel=kernel_flux_x, dim=1, inputs=[q_wp, flux_x_wp, params], device=device)
    wp.launch(kernel=kernel_flux_y, dim=1, inputs=[q_wp, flux_y_wp, params], device=device)
    
    expected_flux = flux_x_wp.numpy()[0] * nx + flux_y_wp.numpy()[0] * ny
    np.testing.assert_allclose(flux_rusanov_wp.numpy()[0], expected_flux, atol=1e-6)

def test_flux_symmetry(device, params):
    """Verify rusanov_flux(q_L, q_R, n) == -rusanov_flux(q_R, q_L, -n)"""
    np.random.seed(42)
    q_l = get_random_state(subsonic=True)
    q_r = get_random_state(subsonic=True)
    
    nx, ny = 0.6, 0.8 # Some normal
    
    q_l_wp = wp.array([q_l], dtype=wp.vec4, device=device)
    q_r_wp = wp.array([q_r], dtype=wp.vec4, device=device)
    
    flux_lr_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    flux_rl_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    
    wp.launch(kernel=kernel_rusanov_flux, dim=1, inputs=[q_l_wp, q_r_wp, nx, ny, flux_lr_wp, params], device=device)
    wp.launch(kernel=kernel_rusanov_flux, dim=1, inputs=[q_r_wp, q_l_wp, -nx, -ny, flux_rl_wp, params], device=device)
    
    np.testing.assert_allclose(flux_lr_wp.numpy()[0], -flux_rl_wp.numpy()[0], atol=1e-6)

@wp.kernel
def kernel_hllc_flux(
    q_l: wp.array(dtype=wp.vec4),
    q_r: wp.array(dtype=wp.vec4),
    nx: float,
    ny: float,
    flux_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    ql = q_l[tid]
    qr = q_r[tid]
    
    # Compute normal fluxes
    fln = euler_kernels.flux_x(ql, params) * nx + euler_kernels.flux_y(ql, params) * ny
    frn = euler_kernels.flux_x(qr, params) * nx + euler_kernels.flux_y(qr, params) * ny
    
    flux_out[tid] = euler_kernels.hllc_flux(ql, qr, fln, frn, nx, ny, params)

class TestRiemannSolvers:
    def test_hllc_consistency(self, device, params):
        """HLLC flux with identical states should match physical flux."""
        q_val = [1.0, 1.2, 0.0, 3.22] # rho=1, u=1.2, p=1
        q_wp = wp.array([q_val], dtype=wp.vec4, device=device)
        flux_out_wp = wp.zeros(1, dtype=wp.vec4, device=device)
        
        nx, ny = 1.0, 0.0
        wp.launch(kernel=kernel_hllc_flux, dim=1, inputs=[q_wp, q_wp, nx, ny, flux_out_wp, params], device=device)
        
        # Expected physical flux Fx
        flux_x_wp = wp.zeros(1, dtype=wp.vec4, device=device)
        wp.launch(kernel=kernel_flux_x, dim=1, inputs=[q_wp, flux_x_wp, params], device=device)
        
        np.testing.assert_allclose(flux_out_wp.numpy()[0], flux_x_wp.numpy()[0], atol=1e-6)

    def test_hllc_supersonic_upwinding(self, device, params):
        """For supersonic flow to the right, HLLC should return flux(Q_left)."""
        # Mach number >> 1
        q_l = [1.0, 5.0, 0.0, 1.0/0.4 + 0.5*25.0] # rho=1, u=5, p=1
        q_r = [1.0, 4.0, 0.0, 1.0/0.4 + 0.5*16.0] # rho=1, u=4, p=1
        
        q_l_wp = wp.array([q_l], dtype=wp.vec4, device=device)
        q_r_wp = wp.array([q_r], dtype=wp.vec4, device=device)
        flux_out_wp = wp.zeros(1, dtype=wp.vec4, device=device)
        
        nx, ny = 1.0, 0.0
        wp.launch(kernel=kernel_hllc_flux, dim=1, inputs=[q_l_wp, q_r_wp, nx, ny, flux_out_wp, params], device=device)
        
        # Expected: Fx(q_l)
        flux_x_l_wp = wp.zeros(1, dtype=wp.vec4, device=device)
        wp.launch(kernel=kernel_flux_x, dim=1, inputs=[q_l_wp, flux_x_l_wp, params], device=device)
        
        np.testing.assert_allclose(flux_out_wp.numpy()[0], flux_x_l_wp.numpy()[0], atol=1e-6)

    def test_hllc_contact_discontinuity(self, device, params):
        """Preserve contact discontinuity (same p and u, different rho)."""
        # rho_l=1.0, rho_r=2.0, u=1.0, p=1.0
        q_l = [1.0, 1.0, 0.0, 2.5 + 0.5]
        q_r = [2.0, 2.0, 0.0, 2.5 + 1.0]
        
        q_l_wp = wp.array([q_l], dtype=wp.vec4, device=device)
        q_r_wp = wp.array([q_r], dtype=wp.vec4, device=device)
        flux_out_wp = wp.zeros(1, dtype=wp.vec4, device=device)
        
        nx, ny = 1.0, 0.0
        wp.launch(kernel=kernel_hllc_flux, dim=1, inputs=[q_l_wp, q_r_wp, nx, ny, flux_out_wp, params], device=device)
        
        # Across a contact moving at S*=u_l=u_r, the pressure should be preserved.
        # Numerical flux for rho should be rho_l * S* (if S* > 0)
        # F_rho = 1.0 * 1.0 = 1.0
        # F_rhou = 1.0 * 1.0^2 + 1.0 = 2.0
        # F_rhov = 0.0
        # F_E = (E_l + p) * u = (3.0 + 1.0) * 1.0 = 4.0
        
        expected = np.array([1.0, 2.0, 0.0, 4.0])
        np.testing.assert_allclose(flux_out_wp.numpy()[0], expected, atol=1e-6)
