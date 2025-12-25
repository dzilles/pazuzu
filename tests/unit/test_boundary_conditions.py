import unittest
import numpy as np
import warp as wp
import sys
import os
from typing import Any

# Add project root (parent of src) to path to allow 'src.' imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.kernels import boundary_conditions as bc
from src.kernels.structs import EquationParams32, BoundaryState32

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
    bc_data: wp.array(dtype=BoundaryState32),
    nx: float,
    ny: float,
    t: float,
    ramp_time: float,
    q_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    bc_state = bc_data[0]
    q_out[tid] = bc.apply_inlet(q_inner[tid], nx, ny, t, ramp_time, bc_state, params)

@wp.kernel
def kernel_outlet(
    q_inner: wp.array(dtype=wp.vec4),
    bc_data: wp.array(dtype=BoundaryState32),
    nx: float,
    ny: float,
    q_out: wp.array(dtype=wp.vec4),
    params: Any
):
    tid = wp.tid()
    bc_state = bc_data[0]
    q_out[tid] = bc.apply_outlet(q_inner[tid], nx, ny, bc_state, params)

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
        cls.params.rho_inf = 1.0
        cls.params.u_inf = 0.0
        cls.params.v_inf = 0.0
        cls.params.p_inf = 1.0

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
        # target (freestream): rho=1, u=1, v=0, p=1 => E = 3.0
        q_inner = [1.0, -2.0, 0.0, 4.0] 
        q_inner_wp = wp.array([q_inner], dtype=wp.vec4, device=self.device)
        q_out_wp = wp.zeros(1, dtype=wp.vec4, device=self.device)
        
        # BC Data: type=INLET, v0=rho=1, v1=u=1, v2=v=0, v3=p=1
        bc_data_np = np.zeros(1, dtype=BoundaryState32.numpy_dtype())
        bc_data_np[0]['type'] = bc.BC_INLET
        bc_data_np[0]['v0'] = 1.0
        bc_data_np[0]['v1'] = 1.0
        bc_data_np[0]['v2'] = 0.0
        bc_data_np[0]['v3'] = 1.0
        bc_data_wp = wp.array(bc_data_np, dtype=BoundaryState32, device=self.device)

        wp.launch(
            kernel=kernel_inlet,
            dim=1,
            inputs=[q_inner_wp, bc_data_wp, 1.0, 0.0, 2.0, 1.0, q_out_wp, self.params],
            device=self.device
        )
        
        # rho=1, u=1, v=0, p=1 => E = 2.5 + 0.5 = 3.0
        expected = np.array([1.0, 1.0, 0.0, 3.0])
        np.testing.assert_allclose(q_out_wp.numpy()[0], expected, atol=1e-6)

    def test_outlet_subsonic(self):
        # Subsonic outflow
        # rho=1, u=0.5, v=0, p=1 => un_i = 0.5, c_i = 1.18. Subsonic.
        # target: p=1, same rho, u, v
        q_inner = [1.0, 0.5, 0.0, 2.5 + 0.125] # rho=1, u=0.5, v=0, p=1
        q_inner_wp = wp.array([q_inner], dtype=wp.vec4, device=self.device)
        q_out_wp = wp.zeros(1, dtype=wp.vec4, device=self.device)
        
        # BC Data: type=OUTLET, v0=p_back=1
        bc_data_np = np.zeros(1, dtype=BoundaryState32.numpy_dtype())
        bc_data_np[0]['type'] = bc.BC_OUTLET
        bc_data_np[0]['v0'] = 1.0
        bc_data_wp = wp.array(bc_data_np, dtype=BoundaryState32, device=self.device)

        wp.launch(
            kernel=kernel_outlet,
            dim=1,
            inputs=[q_inner_wp, bc_data_wp, 1.0, 0.0, q_out_wp, self.params],
            device=self.device
        )
        
        expected = np.array([1.0, 0.5, 0.0, 2.625])
        np.testing.assert_allclose(q_out_wp.numpy()[0], expected, atol=1e-6)

if __name__ == '__main__':
    unittest.main()
