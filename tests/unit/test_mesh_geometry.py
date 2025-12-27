import numpy as np
import pytest
import warp as wp
from src.geometry.mesh import Mesh
from src.core.basis import Basis
from tests.test_utils import create_dummy_mesh

@pytest.fixture(autouse=True)
def init_warp():
    """Ensure Warp is initialized for tests."""
    wp.init()

def test_metric_calculation():
    """
    Test metric calculation for a single rectangular element.
    Rectangle: (0,0) to (2,1).
    Area = 2.0.
    Mapping from [-1, 1]x[-1, 1]:
    x = xi + 1  => dx/dxi = 1, dx/deta = 0
    y = 0.5*eta + 0.5 => dy/dxi = 0, dy/deta = 0.5
    J = dx/dxi * dy/deta - dx/deta * dy/dxi = 1 * 0.5 - 0 = 0.5.
    """
    device = "cpu"
    # Create 1x1 element mesh covering [0,2]x[0,1]
    mesh = create_dummy_mesh(nx=1, ny=1, x_min=0, x_max=2, y_min=0, y_max=1, device=device)
    
    # Use P=1 basis
    basis = Basis(polynomial_degree=1, device=device)
    mesh.compute_geometry(basis)
    
    # 1. Check Jacobian
    # J_host is (num_elements, Np). For P=1, Np=4.
    expected_J = 0.5
    np.testing.assert_allclose(mesh.J_host[0, :], expected_J, atol=1e-7)
    
    # 2. Check Volume (Area)
    # Area = sum(weights * J). sum(weights) = 2*2 = 4 for 2D. 4 * 0.5 = 2.0.
    assert pytest.approx(mesh.vol_host[0]) == 2.0
    
    # 3. Check Normals (face_geo_factors: nx, ny, J_surf)
    # Face 0: Bottom (y=0). Normal points down: (0, -1). Length = dx = 2.
    # Face 1: Right (x=2). Normal points right: (1, 0). Length = dy = 1.
    # Face 2: Top (y=1). Normal points up: (0, 1). Length = dx = 2.
    # Face 3: Left (x=0). Normal points left: (-1, 0). Length = dy = 1.
    
    # Face 0
    np.testing.assert_allclose(mesh.face_geo_factors_host[0, 0, :2], [0.0, -1.0], atol=1e-7)
    np.testing.assert_allclose(mesh.face_geo_factors_host[0, 0, 2], 1.0, atol=1e-7)
    
    # Face 1
    np.testing.assert_allclose(mesh.face_geo_factors_host[0, 1, :2], [1.0, 0.0], atol=1e-7)
    np.testing.assert_allclose(mesh.face_geo_factors_host[0, 1, 2], 0.5, atol=1e-7)
    
    # Face 2
    np.testing.assert_allclose(mesh.face_geo_factors_host[0, 2, :2], [0.0, 1.0], atol=1e-7)
    np.testing.assert_allclose(mesh.face_geo_factors_host[0, 2, 2], 1.0, atol=1e-7)
    
    # Face 3
    np.testing.assert_allclose(mesh.face_geo_factors_host[0, 3, :2], [-1.0, 0.0], atol=1e-7)
    np.testing.assert_allclose(mesh.face_geo_factors_host[0, 3, 2], 0.5, atol=1e-7)

def test_periodic_linking():
    """
    Test periodic boundary condition linking for a 3x1 strip.
    """
    device = "cpu"
    # 3 elements in X, 1 in Y.
    # Elements: 0, 1, 2.
    # Connectivity of 0: [None, 1, None, Tag3]
    # Connectivity of 2: [None, Tag4, None, 1]
    mesh = create_dummy_mesh(nx=3, ny=1, x_min=0, x_max=3, y_min=0, y_max=1, device=device)
    
    # Tags: Left=3, Right=4.
    tag_left = 3
    tag_right = 4
    
    # Verify initial state
    assert mesh.boundary_tags_host[0, 3] == tag_left
    assert mesh.boundary_tags_host[2, 1] == tag_right
    assert mesh.connectivity_host[0, 3, 0] == -1
    assert mesh.connectivity_host[2, 1, 0] == -1
    
    # Apply periodicity in X
    mesh.apply_periodic_condition(tag1=tag_left, tag2=tag_right, axis='x')
    
    # Verify linking
    # Element 0, Face 3 (Left) should link to Element 2, Face 1 (Right)
    assert mesh.connectivity_host[0, 3, 0] == 2
    assert mesh.connectivity_host[0, 3, 1] == 1
    
    # Element 2, Face 1 (Right) should link to Element 0, Face 3 (Left)
    assert mesh.connectivity_host[2, 1, 0] == 0
    assert mesh.connectivity_host[2, 1, 1] == 3
    
    # Verify tags are cleared
    assert mesh.boundary_tags_host[0, 3] == 0
    assert mesh.boundary_tags_host[2, 1] == 0
