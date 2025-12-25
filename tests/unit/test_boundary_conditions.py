import unittest
import numpy as np
import warp as wp
import sys
import os
from typing import Any

# Add project root (parent of src) to path to allow 'src.' imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.kernels import boundary_conditions as bc
from src.kernels.structs import EquationParams32

@wp.kernel
def kernel_slip_wall(
    q_inner: wp.array(dtype=wp.vec4),
    nx: float,
    ny: float,
    q_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    q_out[tid] = bc.apply_slip_wall(q_inner[tid], nx, ny, params)

@wp.kernel
def kernel_inlet(
    q_inner: wp.array(dtype=wp.vec4),
    nx: float,
    ny: float,
    t: float,
    ramp_time: float,
    q_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    q_out[tid] = bc.apply_inlet(q_inner[tid], nx, ny, t, ramp_time, params)

@wp.kernel
def kernel_outlet(
    q_inner: wp.array(dtype=wp.vec4),
    nx: float,
    ny: float,
    q_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    q_out[tid] = bc.apply_outlet(q_inner[tid], nx, ny, params)

class TestBoundaryConditions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        wp.init()
        cls.device = "cpu"
        cls.gamma = 1.4
        
        cls.params = EquationParams32()
        cls.params.gamma = 1.4
        cls.params.rho_floor = 1e-6
        cls.params.p_floor = 1e-6
        cls.params.half = 0.5
        cls.params.one = 1.0

    def test_slip_wall(self):
        # Normal in x: (1, 0)
        # Inner velocity (1, 1) => rhou=1, rhov=1
        # Expect reflected rhou = -1, rhov = 1
        q_inner = [1.0, 1.0, 1.0, 3.0]
        q_inner_wp = wp.array([q_inner], dtype=wp.vec4, device=self.device)
        q_out_wp = wp.zeros(1, dtype=wp.vec4, device=self.device)
        
        wp.launch(
            kernel=kernel_slip_wall,
            dim=1,
            inputs=[q_inner_wp, 1.0, 0.0, q_out_wp, self.params],
            device=self.device
        )
        
        expected = np.array([1.0, -1.0, 1.0, 3.0])
        np.testing.assert_allclose(q_out_wp.numpy()[0], expected, atol=1e-6)

    def test_inlet_full(self):
        # t=2.0, ramp_time=1.0 => target u should be 1.0
        # If we use supersonic inflow (un_i < -c_i), we should get exactly target state
        # rho=1, u=-2, v=0, p=1 => un_i = -2, c_i = sqrt(1.4) ~ 1.18. 
        # un_i < -c_i is True.
        # target (freestream): rho=1, u=1, v=0, p=1 => E = 3.0
        q_inner = [1.0, -2.0, 0.0, 4.0] 
        q_inner_wp = wp.array([q_inner], dtype=wp.vec4, device=self.device)
        q_out_wp = wp.zeros(1, dtype=wp.vec4, device=self.device)
        
        wp.launch(
            kernel=kernel_inlet,
            dim=1,
            inputs=[q_inner_wp, 1.0, 0.0, 2.0, 1.0, q_out_wp, self.params],
            device=self.device
        )
        
        # rho=1, u=1, v=0, p=1 => E = 2.5 + 0.5 = 3.0
        expected = np.array([1.0, 1.0, 0.0, 3.0])
        np.testing.assert_allclose(q_out_wp.numpy()[0], expected, atol=1e-6)

    def test_outlet_subsonic(self):
        # Subsonic outflow
        # rho=1, u=0.5, v=0, p=1 => un_i = 0.5, c_i = 1.18. Subsonic.
        # target: p=1, same rho, u, v
        # j_inner = 0.5 + 2*1.18/0.4 = 0.5 + 5.9 = 6.4
        # j_target = 0.5 - 2*1.18/0.4 = 0.5 - 5.9 = -5.4
        # un_b = 0.5 * (6.4 - 5.4) = 0.5
        # c_b = 0.25 * 0.4 * (6.4 + 5.4) = 0.1 * 11.8 = 1.18
        # Since un_i > 0, we use s_inner. s_inner = 1/1^1.4 = 1.
        # Reconstruct: rho_b = (1.18^2 / (1.4*1))^(1/0.4) = (1.4/1.4)^2.5 = 1.
        # p_b = 1 * 1^1.4 = 1.
        # So we should get back roughly the same state.
        
        q_inner = [1.0, 0.5, 0.0, 2.5 + 0.125] # rho=1, u=0.5, v=0, p=1
        q_inner_wp = wp.array([q_inner], dtype=wp.vec4, device=self.device)
        q_out_wp = wp.zeros(1, dtype=wp.vec4, device=self.device)
        
        wp.launch(
            kernel=kernel_outlet,
            dim=1,
            inputs=[q_inner_wp, 1.0, 0.0, q_out_wp, self.params],
            device=self.device
        )
        
        expected = np.array([1.0, 0.5, 0.0, 2.625])
        np.testing.assert_allclose(q_out_wp.numpy()[0], expected, atol=1e-6)

if __name__ == '__main__':
    unittest.main()
