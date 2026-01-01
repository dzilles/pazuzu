import pytest
import numpy as np
import warp as wp
from src.kernels import navier_stokes_kernels as nsk
from src.kernels.structs import EquationParams32
from src.physics.laws import navier_stokes as nsl
from typing import Any

# Helper kernels to test wp.func
@wp.kernel
def kernel_viscous_flux(
    q: wp.array(dtype=wp.vec4),
    grad_u: wp.array(dtype=wp.vec2),
    grad_v: wp.array(dtype=wp.vec2),
    grad_T: wp.array(dtype=wp.vec2),
    flux_x_out: wp.array(dtype=wp.vec4),
    flux_y_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    flux_x_out[tid] = nsl.viscous_flux_x(q[tid], grad_u[tid], grad_v[tid], grad_T[tid], params)
    flux_y_out[tid] = nsl.viscous_flux_y(q[tid], grad_u[tid], grad_v[tid], grad_T[tid], params)

class TestNavierStokesKernels:
        
    def test_primitive_gradients_constant(self, device):
        """
        Verifies that gradients of a constant field are zero.
        """
        wp.init()
        
        # Setup 1 Element, N=1 (4 nodes)
        num_elements = 1
        N = 1
        Np = (N+1)**2
        
        # Constant field
        rho = 1.0
        u = 2.0
        v = 0.5
        p = 1.0
        # E = p/(gamma-1) + 0.5*rho*(u^2+v^2)
        # We don't need E for gradient of primitives u, v, T.
        # But we pass Q.
        
        q_host = np.zeros((num_elements, Np, 4), dtype=np.float32)
        q_host[:, :, 0] = rho
        q_host[:, :, 1] = rho * u
        q_host[:, :, 2] = rho * v
        q_host[:, :, 3] = 10.0 # Arbitrary E
        
        q = wp.array(q_host, dtype=wp.vec4, device=device)
        
        grad_u = wp.zeros((num_elements, Np), dtype=wp.vec2, device=device)
        grad_v = wp.zeros((num_elements, Np), dtype=wp.vec2, device=device)
        grad_T = wp.zeros((num_elements, Np), dtype=wp.vec2, device=device)
        
        # Identity Metrics (dx/dr = 1) -> Dr = Dx
        # Use simple D matrices (zero for constant check? No, actual D)
        # But if field is constant, D * q = 0 regardless of metric.
        Dr = wp.array(np.zeros((Np, Np), dtype=np.float32), dtype=wp.float32, device=device)
        Ds = wp.array(np.zeros((Np, Np), dtype=np.float32), dtype=wp.float32, device=device)
        
        # Metrics don't matter if derivative is zero
        rx = wp.array(np.ones((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=device)
        ry = wp.array(np.zeros((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=device)
        sx = wp.array(np.zeros((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=device)
        sy = wp.array(np.ones((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=device)
        
        params = EquationParams32()
        params.gamma = 1.4
        params.rho_floor = 1e-5
        params.p_floor = 1e-5
        params.one = 1.0
        params.gas_constant = 1.0
        
        wp.launch(
            kernel=nsk.compute_primitive_gradients_volume,
            dim=num_elements * Np,
            inputs=[q, grad_u, grad_v, grad_T, Dr, Ds, rx, ry, sx, sy, Np, params],
            device=device
        )
        # Check
        assert np.allclose(grad_u.numpy(), 0.0)
        assert np.allclose(grad_v.numpy(), 0.0)
        # T depends on P and Rho. Both constant.
        assert np.allclose(grad_T.numpy(), 0.0)

    def test_viscous_flux_calculation(self, device):
        """Test viscous flux calculation for zero gradients and shear flow."""
        params = EquationParams32()
        params.gamma = 1.4
        params.rho_floor = 1e-6
        params.p_floor = 1e-6
        params.one = 1.0
        params.mu = 0.1
        params.cp = 1000.0
        params.prandtl = 0.72
        params.gas_constant = 287.0
        
        # Case 1: Zero Gradients
        q_val = [1.0, 1.0, 0.0, 2.5] # rho=1, u=1, v=0, p=1
        q_wp = wp.array([q_val], dtype=wp.vec4, device=device)
        grad_u_wp = wp.zeros(1, dtype=wp.vec2, device=device)
        grad_v_wp = wp.zeros(1, dtype=wp.vec2, device=device)
        grad_T_wp = wp.zeros(1, dtype=wp.vec2, device=device)
        
        flux_x_out = wp.zeros(1, dtype=wp.vec4, device=device)
        flux_y_out = wp.zeros(1, dtype=wp.vec4, device=device)
        
        wp.launch(
            kernel=kernel_viscous_flux,
            dim=1,
            inputs=[q_wp, grad_u_wp, grad_v_wp, grad_T_wp, flux_x_out, flux_y_out, params],
            device=device
        )
        
        # All viscous fluxes should be zero
        np.testing.assert_allclose(flux_x_out.numpy()[0], 0.0, atol=1e-7)
        np.testing.assert_allclose(flux_y_out.numpy()[0], 0.0, atol=1e-7)
        
        # Case 2: Shear Flow
        # u = y => grad_u = (du/dx, du/dy) = (0, 1)
        # v = 0 => grad_v = (0, 0)
        # T = const => grad_T = (0, 0)
        # mu = 0.1
        # tau_xy = mu * (du/dy + dv/dx) = 0.1 * (1 + 0) = 0.1
        # tau_yx = tau_xy = 0.1
        # tau_xx = 2*mu*(du/dx - 1/3*(du/dx+dv/dy)) = 2*0.1*(0 - 0) = 0
        # tau_yy = 2*mu*(dv/dy - 1/3*(du/dx+dv/dy)) = 0
        
        grad_u_val = [0.0, 1.0]
        grad_u_wp = wp.array([grad_u_val], dtype=wp.vec2, device=device)
        
        wp.launch(
            kernel=kernel_viscous_flux,
            dim=1,
            inputs=[q_wp, grad_u_wp, grad_v_wp, grad_T_wp, flux_x_out, flux_y_out, params],
            device=device
        )
        
        # Fv = [0, tau_xx, tau_xy, u*tau_xx + v*tau_xy - qx]
        # Fv = [0, 0, 0.1, 1*0 + 0*0.1 - 0] = [0, 0, 0.1, 0]
        expected_x = np.array([0.0, 0.0, 0.1, 0.0])
        np.testing.assert_allclose(flux_x_out.numpy()[0], expected_x, atol=1e-7)
        
        # Gv = [0, tau_yx, tau_yy, u*tau_yx + v*tau_yy - qy]
        # Gv = [0, 0.1, 0, 1*0.1 + 0*0 - 0] = [0, 0.1, 0, 0.1]
        expected_y = np.array([0.0, 0.1, 0.0, 0.1])
        np.testing.assert_allclose(flux_y_out.numpy()[0], expected_y, atol=1e-7)
