import numpy as np
import pytest
import warp as wp
import os
import sys
from typing import Any

# Add project root (parent of src) to path to allow 'src.' imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.kernels import boundary_conditions as bc
from src.kernels.structs import EquationParams32, BoundaryState32

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
    p.rho_inf = 1.0
    p.u_inf = 0.0
    p.v_inf = 0.0
    p.p_inf = 1.0
    return p

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

@pytest.mark.parametrize("nx, ny, u_in, v_in", [
    (1.0, 0.0, 1.0, 1.0),   # Normal in x, velocity (1,1) -> reflected (-1,1)
    (0.0, 1.0, 1.0, 1.0),   # Normal in y, velocity (1,1) -> reflected (1,-1)
    (0.7071, 0.7071, 1.0, 0.0), # Diagonal normal
])
def test_slip_wall(nx, ny, u_in, v_in, device, params):
    # rho=1, p=1 => E = 2.5 + 0.5*(u^2 + v^2)
    rho = 1.0
    E = 2.5 + 0.5 * rho * (u_in**2 + v_in**2)
    q_inner = [rho, rho * u_in, rho * v_in, E]
    q_inner_wp = wp.array([q_inner], dtype=wp.vec4, device=device)
    q_out_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    
    wp.launch(
        kernel=kernel_slip_wall,
        dim=1,
        inputs=[q_inner_wp, nx, ny, q_out_wp, params],
        device=device
    )
    
    q_out = q_out_wp.numpy()[0]
    u_out = q_out[1] / q_out[0]
    v_out = q_out[2] / q_out[0]
    
    # Velocity reflection: V_out = V_in - 2 * (V_in . n) * n
    v_dot_n = u_in * nx + v_in * ny
    u_expected = u_in - 2.0 * v_dot_n * nx
    v_expected = v_in - 2.0 * v_dot_n * ny
    
    np.testing.assert_allclose(u_out, u_expected, atol=1e-6)
    np.testing.assert_allclose(v_out, v_expected, atol=1e-6)
    assert pytest.approx(q_out[0]) == q_inner[0] # Density should be same
    assert pytest.approx(q_out[3]) == q_inner[3] # Energy should be same

@pytest.mark.parametrize("rho_target, u_target, v_target, p_target", [
    (1.0, -2.0, 0.0, 1.0), # Supersonic inflow (un < -c)
    (1.2, -3.0, 0.2, 1.5), # Supersonic inflow
])
def test_inlet_full(rho_target, u_target, v_target, p_target, device, params):
    # t=2.0, ramp_time=1.0 => fully ramped
    # Use a supersonic inflow state for q_inner too, or just anything that makes un < -c
    q_inner = [1.0, -2.0, 0.0, 5.0] 
    q_inner_wp = wp.array([q_inner], dtype=wp.vec4, device=device)
    q_out_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    
    bc_data_np = np.zeros(1, dtype=BoundaryState32.numpy_dtype())
    bc_data_np[0]['type'] = bc.BC_INLET
    bc_data_np[0]['v0'] = rho_target
    bc_data_np[0]['v1'] = u_target
    bc_data_np[0]['v2'] = v_target
    bc_data_np[0]['v3'] = p_target
    bc_data_wp = wp.array(bc_data_np, dtype=BoundaryState32, device=device)

    wp.launch(
        kernel=kernel_inlet,
        dim=1,
        inputs=[q_inner_wp, bc_data_wp, 1.0, 0.0, 2.0, 1.0, q_out_wp, params],
        device=device
    )
    
    E_target = p_target / (params.gamma - 1.0) + 0.5 * rho_target * (u_target**2 + v_target**2)
    expected = np.array([rho_target, rho_target * u_target, rho_target * v_target, E_target])
    np.testing.assert_allclose(q_out_wp.numpy()[0], expected, atol=1e-6)

def test_outlet_subsonic(device, params):
    # Subsonic outflow
    # rho=1, u=0.5, v=0, p=1 => un_i = 0.5, c_i = 1.18. Subsonic.
    # target: p=1, same rho, u, v
    q_inner = [1.0, 0.5, 0.0, 2.5 + 0.125] # rho=1, u=0.5, v=0, p=1
    q_inner_wp = wp.array([q_inner], dtype=wp.vec4, device=device)
    q_out_wp = wp.zeros(1, dtype=wp.vec4, device=device)
    
    # BC Data: type=OUTLET, v0=p_back=1
    bc_data_np = np.zeros(1, dtype=BoundaryState32.numpy_dtype())
    bc_data_np[0]['type'] = bc.BC_OUTLET
    bc_data_np[0]['v0'] = 1.0
    bc_data_wp = wp.array(bc_data_np, dtype=BoundaryState32, device=device)

    wp.launch(
        kernel=kernel_outlet,
        dim=1,
        inputs=[q_inner_wp, bc_data_wp, 1.0, 0.0, q_out_wp, params],
        device=device
    )
    
    expected = np.array([1.0, 0.5, 0.0, 2.625])
    np.testing.assert_allclose(q_out_wp.numpy()[0], expected, atol=1e-6)