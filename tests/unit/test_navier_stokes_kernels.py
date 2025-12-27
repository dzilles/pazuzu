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
    @pytest.fixture(autouse=True)
    def setup(self):
        wp.init()
        self.device = "cpu"
        
    def test_primitive_gradients_constant(self):
        """Test that a constant state results in zero gradients."""
        Np = 4
        num_elements = 1
        
        q_host = np.zeros((num_elements, Np, 4), dtype=np.float32)
        q_host[:, :, 0] = 1.0 # rho
        q_host[:, :, 1] = 0.5 # rho*u => u=0.5
        q_host[:, :, 2] = 0.0 # rho*v => v=0.0
        q_host[:, :, 3] = 2.5 # E
        
        q = wp.array(q_host, dtype=wp.vec4, device=self.device)
        grad_u = wp.zeros((num_elements, Np), dtype=wp.vec2, device=self.device)
        grad_v = wp.zeros((num_elements, Np), dtype=wp.vec2, device=self.device)
        grad_T = wp.zeros((num_elements, Np), dtype=wp.vec2, device=self.device)
        
        # Identity differentiation matrices for constant test
        Dr = wp.array(np.zeros((Np, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        Ds = wp.array(np.zeros((Np, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        
        rx = wp.array(np.ones((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        ry = wp.array(np.zeros((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        sx = wp.array(np.zeros((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        sy = wp.array(np.ones((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        
        params = EquationParams32()
        params.gamma = 1.4
        params.rho_floor = 1e-5
        params.p_floor = 1e-5
        params.one = 1.0
        params.cp = 1004.5
        params.gas_constant = 287.0
        
        wp.launch(
            kernel=nsk.compute_primitive_gradients_volume,
            dim=(num_elements, Np),
            inputs=[q, grad_u, grad_v, grad_T, Dr, Ds, rx, ry, sx, sy, Np, params],
            device=self.device
        )
        
        res_u = grad_u.numpy()
        res_v = grad_v.numpy()
        res_T = grad_T.numpy()
        
        assert np.allclose(res_u, 0.0)
        assert np.allclose(res_v, 0.0)
        assert np.allclose(res_T, 0.0)

    def test_viscous_flux_calculation(self):
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
        q_wp = wp.array([q_val], dtype=wp.vec4, device=self.device)
        grad_u_wp = wp.zeros(1, dtype=wp.vec2, device=self.device)
        grad_v_wp = wp.zeros(1, dtype=wp.vec2, device=self.device)
        grad_T_wp = wp.zeros(1, dtype=wp.vec2, device=self.device)
        
        flux_x_out = wp.zeros(1, dtype=wp.vec4, device=self.device)
        flux_y_out = wp.zeros(1, dtype=wp.vec4, device=self.device)
        
        wp.launch(
            kernel=kernel_viscous_flux,
            dim=1,
            inputs=[q_wp, grad_u_wp, grad_v_wp, grad_T_wp, flux_x_out, flux_y_out, params],
            device=self.device
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
        grad_u_wp = wp.array([grad_u_val], dtype=wp.vec2, device=self.device)
        
        wp.launch(
            kernel=kernel_viscous_flux,
            dim=1,
            inputs=[q_wp, grad_u_wp, grad_v_wp, grad_T_wp, flux_x_out, flux_y_out, params],
            device=self.device
        )
        
        # Fv = [0, tau_xx, tau_xy, u*tau_xx + v*tau_xy - qx]
        # Fv = [0, 0, 0.1, 1*0 + 0*0.1 - 0] = [0, 0, 0.1, 0]
        expected_x = np.array([0.0, 0.0, 0.1, 0.0])
        np.testing.assert_allclose(flux_x_out.numpy()[0], expected_x, atol=1e-7)
        
        # Gv = [0, tau_yx, tau_yy, u*tau_yx + v*tau_yy - qy]
        # Gv = [0, 0.1, 0, 1*0.1 + 0*0 - 0] = [0, 0.1, 0, 0.1]
        expected_y = np.array([0.0, 0.1, 0.0, 0.1])
        np.testing.assert_allclose(flux_y_out.numpy()[0], expected_y, atol=1e-7)
