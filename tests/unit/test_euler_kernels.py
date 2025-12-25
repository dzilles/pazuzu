import unittest
import numpy as np
import warp as wp
import sys
import os

# Add project root (parent of src) to path to allow 'src.' imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.kernels import euler_kernels

# Helper kernels to test wp.func
@wp.kernel
def kernel_pressure(
    q: wp.array(dtype=wp.vec4),
    p_out: wp.array(dtype=float),
    gamma: float,
    rho_floor: float,
    p_floor: float
):
    tid = wp.tid()
    p_out[tid] = euler_kernels.pressure(q[tid], gamma, rho_floor, p_floor, 0.5, 1.0)

@wp.kernel
def kernel_flux_x(
    q: wp.array(dtype=wp.vec4),
    flux_out: wp.array(dtype=wp.vec4),
    gamma: float,
    rho_floor: float,
    p_floor: float
):
    tid = wp.tid()
    flux_out[tid] = euler_kernels.flux_x(q[tid], gamma, rho_floor, p_floor, 0.5, 1.0)

@wp.kernel
def kernel_rusanov_flux(
    q_l: wp.array(dtype=wp.vec4),
    q_r: wp.array(dtype=wp.vec4),
    nx: float,
    ny: float,
    flux_out: wp.array(dtype=wp.vec4),
    gamma: float,
    rho_floor: float,
    p_floor: float
):
    tid = wp.tid()
    flux_out[tid] = euler_kernels.rusanov_flux(q_l[tid], q_r[tid], nx, ny, gamma, rho_floor, p_floor, 0.5, 1.0)

class TestEulerKernels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        wp.init()
        cls.device = "cpu"
        cls.gamma = 1.4
        cls.rho_floor = 1e-6
        cls.p_floor = 1e-6

    def test_pressure_calculation(self):
        # rho=1, u=1, v=0, p=1 => E = p/(gamma-1) + 0.5*rho*u^2 = 1/0.4 + 0.5 = 2.5 + 0.5 = 3.0
        q_val = [1.0, 1.0, 0.0, 3.0]
        q_wp = wp.array([q_val], dtype=wp.vec4, device=self.device)
        p_out_wp = wp.zeros(1, dtype=float, device=self.device)
        
        wp.launch(
            kernel=kernel_pressure,
            dim=1,
            inputs=[q_wp, p_out_wp, self.gamma, self.rho_floor, self.p_floor],
            device=self.device
        )
        
        self.assertAlmostEqual(p_out_wp.numpy()[0], 1.0, places=6)

    def test_flux_x_calculation(self):
        # rho=1, u=2, v=0, p=1
        # q = [1, 2, 0, E]
        # E = 1/0.4 + 0.5*1*2^2 = 2.5 + 2 = 4.5
        q_val = [1.0, 2.0, 0.0, 4.5]
        q_wp = wp.array([q_val], dtype=wp.vec4, device=self.device)
        flux_out_wp = wp.zeros(1, dtype=wp.vec4, device=self.device)
        
        wp.launch(
            kernel=kernel_flux_x,
            dim=1,
            inputs=[q_wp, flux_out_wp, self.gamma, self.rho_floor, self.p_floor],
            device=self.device
        )
        
        # Flux F = [rho*u, rho*u^2 + p, rho*u*v, (E+p)*u]
        # F = [2, 1*4 + 1, 0, (4.5+1)*2] = [2, 5, 0, 11]
        expected = np.array([2.0, 5.0, 0.0, 11.0])
        np.testing.assert_allclose(flux_out_wp.numpy()[0], expected, atol=1e-6)

    def test_rusanov_flux_symmetric(self):
        # Symmetric left and right states should give average flux if wave speed is zero? 
        # No, if states are same, flux should be physical flux.
        q_val = [1.0, 0.0, 0.0, 2.5] # rho=1, p=1, u=0, v=0
        q_wp = wp.array([q_val], dtype=wp.vec4, device=self.device)
        flux_out_wp = wp.zeros(1, dtype=wp.vec4, device=self.device)
        
        wp.launch(
            kernel=kernel_rusanov_flux,
            dim=1,
            inputs=[q_wp, q_wp, 1.0, 0.0, flux_out_wp, self.gamma, self.rho_floor, self.p_floor],
            device=self.device
        )
        
        # Expected flux F_x = [0, 1, 0, 0]
        expected = np.array([0.0, 1.0, 0.0, 0.0])
        np.testing.assert_allclose(flux_out_wp.numpy()[0], expected, atol=1e-6)

if __name__ == '__main__':
    unittest.main()
