import pytest
import numpy as np
import meshio
from src.geometry.mesh import Mesh
import warp as wp

@pytest.fixture(scope="module")
def init_warp():
    wp.init()

def test_dx_calculation_from_file(tmp_path, init_warp):
    """
    Verifies that Mesh loads correctly and calculates dx based on the minimum element area.
    """
    # Create vertices for two quads of different sizes
    # Quad 1: 1x1 square (Area = 1.0)
    # Vertices: (0,0), (1,0), (1,1), (0,1)
    
    # Quad 2: 2x2 square (Area = 4.0)
    # Vertices: (1,0), (3,0), (3,2), (1,2)
    # Note: They share the node (1,0) and the vertical line x=1 partially (but nodes don't align perfectly for a conformal mesh)
    # For this test, we don't care about valid connectivity/conformal mesh, just the area calculation.
    
    points = np.array([
        [0.0, 0.0, 0.0], # 0
        [1.0, 0.0, 0.0], # 1
        [1.0, 1.0, 0.0], # 2
        [0.0, 1.0, 0.0], # 3
        [3.0, 0.0, 0.0], # 4
        [3.0, 2.0, 0.0], # 5
        [1.0, 2.0, 0.0], # 6
    ], dtype=np.float64)
    
    cells = [
        ("quad", np.array([
            [0, 1, 2, 3],   # Area 1
            [1, 4, 5, 6]    # Area 4
        ]))
    ]
    
    # Create meshio object
    mesh = meshio.Mesh(
        points,
        cells
    )
    
    filename = str(tmp_path / "test_variable_area.msh")
    mesh.write(filename)
    
    # Load Mesh using our class
    # Suppress output to keep test clean
    m = Mesh(filename=filename, device="cpu")
    
    # Check areas
    # We expect the code to calculate areas 1.0 and 4.0.
    # The new logic should take min_area = 1.0.
    # So dx = sqrt(1.0) = 1.0.
    
    # (Old logic would give sqrt((1+4)/2) = sqrt(2.5) approx 1.5811)
    
    assert m.dx == pytest.approx(1.0, rel=1e-5), f"Expected dx=1.0 (sqrt of min area), but got {m.dx}"
