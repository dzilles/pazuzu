import unittest
import numpy as np
import warp as wp
import sys
import os
from typing import Any

# Add project root (parent of src) to path to allow 'src.' imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.kernels import euler_kernels
from src.kernels.structs import EquationParams32

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
    
    # Compute normal fluxes needed for the new rusanov_flux signature
    fln = euler_kernels.flux_x(ql, params) * nx + euler_kernels.flux_y(ql, params) * ny
    frn = euler_kernels.flux_x(qr, params) * nx + euler_kernels.flux_y(qr, params) * ny
    
    flux_out[tid] = euler_kernels.rusanov_flux(ql, qr, fln, frn, nx, ny, params)

class TestEulerKernels(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        wp.init()
        cls.device = "cpu"
        cls.gamma = 1.4
        cls.rho_floor = 1e-6
        cls.p_floor = 1e-6
        
        cls.params = EquationParams32()
        cls.params.gamma = cls.gamma
        cls.params.rho_floor = cls.rho_floor
        cls.params.p_floor = cls.p_floor
        cls.params.half = 0.5
        cls.params.one = 1.0

    def test_pressure_calculation(self):
        # rho=1, u=1, v=0, p=1 => E = p/(gamma-1) + 0.5*rho*u^2 = 1/0.4 + 0.5 = 2.5 + 0.5 = 3.0
        q_val = [1.0, 1.0, 0.0, 3.0]
        q_wp = wp.array([q_val], dtype=wp.vec4, device=self.device)
        p_out_wp = wp.zeros(1, dtype=float, device=self.device)
        
        wp.launch(
            kernel=kernel_pressure,
            dim=1,
            inputs=[q_wp, p_out_wp, self.params],
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
            inputs=[q_wp, flux_out_wp, self.params],
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
            inputs=[q_wp, q_wp, 1.0, 0.0, flux_out_wp, self.params],
            device=self.device
        )
        
        # Expected flux F_x = [0, 1, 0, 0]
        expected = np.array([0.0, 1.0, 0.0, 0.0])
        np.testing.assert_allclose(flux_out_wp.numpy()[0], expected, atol=1e-6)

if __name__ == '__main__':
    unittest.main()
