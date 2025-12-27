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
    if wp.get_device().is_cuda:
        # For unit tests on CI or local machines, we often want to force CPU
        # but let's just use what's available and hope it works.
        pass

def test_mesh_inversion_detection():
    """Verifies that the Mesh class detects inverted elements during geometry computation."""
    device = "cpu"
    
    # Create a simple 1x1 Cartesian mesh using utility
    mesh = create_dummy_mesh(nx=1, ny=1, x_min=0, x_max=1, y_min=0, y_max=1, device=device)
    
    # Manually invert the element by swapping two vertices
    # Original CCW: [0,0], [1,0], [1,1], [0,1]
    # Inverted (CW): [0,0], [0,1], [1,1], [1,0]
    mesh.vertices_host[0, 1, :] = [0, 1]
    mesh.vertices_host[0, 3, :] = [1, 0]
    
    # Initialize a basis (P=1)
    basis = Basis(polynomial_degree=1, device=device)
    
    # Attempting to compute geometry should raise ValueError
    with pytest.raises(ValueError, match="Mesh contains inverted element at index 0"):
        mesh.compute_geometry(basis)

def test_valid_mesh_computation():
    """Ensures a valid mesh does not raise errors during geometry computation."""
    device = "cpu"
    
    # Create a simple 2x2 Cartesian mesh using utility
    mesh = create_dummy_mesh(nx=2, ny=2, x_min=0, x_max=1, y_min=0, y_max=1, device=device)
    
    # Initialize a basis (P=1)
    basis = Basis(polynomial_degree=1, device=device)
    
    # This should succeed
    mesh.compute_geometry(basis)
    
    # Check that Jacobians are positive
    assert np.all(mesh.J_host > 0.0)
