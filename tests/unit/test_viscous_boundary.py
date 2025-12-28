import pytest
import numpy as np
import warp as wp
import os
import sys
from typing import Any

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.kernels import navier_stokes_kernels as nsk
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
    p.mu = 0.1
    p.cp = 1000.0
    p.prandtl = 0.72
    p.gas_constant = 287.0
    return p

@wp.kernel
def launch_viscous_surface_reflection(
    q: wp.array(dtype=Any, ndim=2),
    grad_u: wp.array(dtype=Any, ndim=2),
    grad_v: wp.array(dtype=Any, ndim=2),
    grad_T: wp.array(dtype=Any, ndim=2),
    rhs: wp.array(dtype=Any, ndim=2),
    connectivity: wp.array(dtype=wp.int32, ndim=2),
    neighbor_face_indices: wp.array(dtype=wp.int32, ndim=2),
    face_map: wp.array(dtype=wp.int32, ndim=2),
    LIFT: wp.array(dtype=Any, ndim=2),
    face_geo_factors: wp.array(dtype=Any, ndim=3), 
    J: wp.array(dtype=Any, ndim=2),     
    bc_mask: wp.array(dtype=wp.int32, ndim=2),
    bc_data: wp.array(dtype=Any, ndim=1),
    coord_x: wp.array(dtype=Any, ndim=2),
    coord_y: wp.array(dtype=Any, ndim=2),
    Nfp: wp.int32,
    params: Any
):
    e = wp.tid()
    nsk.compute_viscous_surface_term(
        e, q, grad_u, grad_v, grad_T, rhs,
        connectivity, neighbor_face_indices, face_map,
        LIFT, face_geo_factors, J,
        bc_mask, bc_data, coord_x, coord_y,
        Nfp, wp.float32(0.0), wp.float32(1.0), params
    )

def test_no_slip_adiabatic_heat_flux(device, params):
    """
    Verifies that for a No-Slip Wall with a temperature gradient normal to the wall,
    the reflected ghost gradient gTo ensures zero net heat flux contribution from the boundary.
    Specifically, Fv_star . n should have zero contribution from the temperature gradient.
    """
    num_elements = 1
    Np = 1
    Nfp = 1
    
    q_host = np.zeros((num_elements, Np, 4), dtype=np.float32)
    q_host[0, 0] = [1.0, 0.0, 0.0, 2.5] # rho=1, u=0, v=0, p=1
    
    grad_u_host = np.zeros((num_elements, Np, 2), dtype=np.float32)
    grad_v_host = np.zeros((num_elements, Np, 2), dtype=np.float32)
    grad_T_host = np.zeros((num_elements, Np, 2), dtype=np.float32)
    grad_T_host[0, 0] = [1.0, 0.0]
    
    q = wp.array(q_host, dtype=wp.vec4, device=device)
    grad_u = wp.array(grad_u_host, dtype=wp.vec2, device=device)
    grad_v = wp.array(grad_v_host, dtype=wp.vec2, device=device)
    grad_T = wp.array(grad_T_host, dtype=wp.vec2, device=device)
    rhs = wp.zeros((num_elements, Np), dtype=wp.vec4, device=device)
    
    connectivity = wp.array(-np.ones((num_elements, 4), dtype=np.int32), dtype=wp.int32, device=device)
    neighbor_face_indices = wp.zeros((num_elements, 4), dtype=wp.int32, device=device)
    face_map = wp.array(np.zeros((4, Nfp), dtype=np.int32), dtype=wp.int32, device=device)
    
    LIFT = wp.array(np.ones((Np, 4*Nfp), dtype=np.float32), dtype=wp.float32, device=device)
    
    # Face 0: nx=1, ny=0, surf_J=1
    face_geo_factors_host = np.zeros((num_elements, 4, 3), dtype=np.float32)
    face_geo_factors_host[0, 0] = [1.0, 0.0, 1.0]
    face_geo_factors = wp.array(face_geo_factors_host, dtype=wp.float32, device=device)
    
    J = wp.array(np.ones((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=device)
    
    bc_mask_host = np.zeros((num_elements, 4), dtype=np.int32)
    bc_mask_host[0, 0] = bc.BC_NO_SLIP_WALL 
    bc_mask = wp.array(bc_mask_host, dtype=wp.int32, device=device)
    
    bc_data_np = np.zeros(1, dtype=BoundaryState32.numpy_dtype())
    bc_data_np[0]['type'] = bc.BC_NO_SLIP_WALL
    bc_data = wp.array(bc_data_np, dtype=BoundaryState32, device=device)
    
    coord_x = wp.zeros((num_elements, Np), dtype=wp.float32, device=device)
    coord_y = wp.zeros((num_elements, Np), dtype=wp.float32, device=device)
    
    wp.launch(
        kernel=launch_viscous_surface_reflection,
        dim=num_elements,
        inputs=[q, grad_u, grad_v, grad_T, rhs, connectivity, neighbor_face_indices, face_map,
                LIFT, face_geo_factors, J, bc_mask, bc_data, coord_x, coord_y, Nfp, params],
        device=device
    )
    
    res_rhs = rhs.numpy()[0, 0]
    # k = mu * cp / Pr = 0.1 * 1000 / 0.72 = 138.888...
    k_expected = params.mu * params.cp / params.prandtl
    assert res_rhs[3] == pytest.approx(-k_expected, rel=1e-5)

def test_no_slip_velocity_gradient_reflection(device, params):
    """
    Verifies that for a No-Slip Wall with a velocity gradient (shear),
    the reflected ghost gradients guo, gvo ensure consistent viscous stresses.
    For Dirichlet=0, we preserve the NORMAL component and flip the TANGENTIAL component.
    g_out = 2*(g_in . n) * n - g_in
    """
    num_elements = 1
    Np = 1; Nfp = 1
    
    # Internal state: rho=1, u=0, v=0, p=1
    q_host = np.zeros((num_elements, Np, 4), dtype=np.float32)
    q_host[0, 0] = [1.0, 0.0, 0.0, 2.5]
    
    # Case: Shear flow near a wall at x=0 (nx=1, ny=0).
    # u = y  => grad_u = [du/dx, du/dy] = [0, 1]
    # v = 0  => grad_v = [0, 0]
    # grad_T = [0, 0]
    grad_u_host = np.zeros((num_elements, Np, 2), dtype=np.float32)
    grad_u_host[0, 0] = [0.0, 1.0] 
    grad_v_host = np.zeros((num_elements, Np, 2), dtype=np.float32)
    grad_T_host = np.zeros((num_elements, Np, 2), dtype=np.float32)
    
    q = wp.array(q_host, dtype=wp.vec4, device=device)
    grad_u = wp.array(grad_u_host, dtype=wp.vec2, device=device)
    grad_v = wp.array(grad_v_host, dtype=wp.vec2, device=device)
    grad_T = wp.array(grad_T_host, dtype=wp.vec2, device=device)
    rhs = wp.zeros((num_elements, Np), dtype=wp.vec4, device=device)
    
    connectivity = wp.array(-np.ones((num_elements, 4), dtype=np.int32), dtype=wp.int32, device=device)
    neighbor_face_indices = wp.zeros((num_elements, 4), dtype=wp.int32, device=device)
    face_map = wp.array(np.zeros((4, Nfp), dtype=np.int32), dtype=wp.int32, device=device)
    LIFT = wp.array(np.ones((Np, 4*Nfp), dtype=np.float32), dtype=wp.float32, device=device)
    
    # Face 0: nx=1, ny=0, surf_J=1
    face_geo_factors_host = np.zeros((num_elements, 4, 3), dtype=np.float32)
    face_geo_factors_host[0, 0] = [1.0, 0.0, 1.0]
    face_geo_factors = wp.array(face_geo_factors_host, dtype=wp.float32, device=device)
    J = wp.array(np.ones((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=device)
    
    bc_mask_host = np.zeros((num_elements, 4), dtype=np.int32)
    bc_mask_host[0, 0] = bc.BC_NO_SLIP_WALL 
    bc_mask = wp.array(bc_mask_host, dtype=wp.int32, device=device)
    
    bc_data_np = np.zeros(1, dtype=BoundaryState32.numpy_dtype())
    bc_data_np[0]['type'] = bc.BC_NO_SLIP_WALL
    bc_data = wp.array(bc_data_np, dtype=BoundaryState32, device=device)
    
    coord_x = wp.zeros((num_elements, Np), dtype=wp.float32, device=device)
    coord_y = wp.zeros((num_elements, Np), dtype=wp.float32, device=device)
    
    wp.launch(
        kernel=launch_viscous_surface_reflection,
        dim=num_elements,
        inputs=[q, grad_u, grad_v, grad_T, rhs, connectivity, neighbor_face_indices, face_map,
                LIFT, face_geo_factors, J, bc_mask, bc_data, coord_x, coord_y, Nfp, params],
        device=device
    )
    
    res_rhs = rhs.numpy()[0, 0]
    # Check Momentum-V (index 2)
    assert res_rhs[2] == pytest.approx(-0.1, rel=1e-5)
    # Energy contribution should be zero
    assert res_rhs[3] == pytest.approx(0.0, abs=1e-7)