import unittest
import numpy as np
from unittest.mock import MagicMock
from src.core.boundary_condition_manager import BoundaryConditionManager
from src.geometry.mesh import Mesh
from src.kernels import boundary_conditions as bc

class MockConfig:
    def __init__(self, boundaries):
        self.boundaries = boundaries

class TestBoundaryConditionManager(unittest.TestCase):
    def test_cartesian_mapping(self):
        # Create a small Cartesian mesh (internally generated)
        # nx=2, ny=2. 
        # Boundary tags should be: Bottom=1, Top=2, Left=3, Right=4
        mesh = Mesh(nx=2, ny=2, device="cpu")
        
        # Config using descriptive names
        # This is what we WANT to work
        boundaries = {
            "Bottom": {"type": "slip_wall"},
            "Top": {"type": "slip_wall"},
            "Left": {"type": "inlet"},
            "Right": {"type": "outlet"}
        }
        config = MockConfig(boundaries)
        
        manager = BoundaryConditionManager(mesh, config)
        bc_mask = manager.setup_boundary_conditions()
        
        # Check specific elements
        # Element 0 (Bottom-Left): Left face (3) should be INLET, Bottom face (0) should be WALL
        # bc_mask is (NumElements, 4)
        
        # Current implementation of Mesh Cartesian connectivity:
        # Faces: 0: Bottom, 1: Right, 2: Top, 3: Left
        
        # Element 0 (0,0): 
        #   Face 0 (Bottom) -> Tag 1
        #   Face 3 (Left)   -> Tag 3
        #   Face 1 (Right)  -> Internal (neighbor)
        #   Face 2 (Top)    -> Internal (neighbor)
        
        # If mapping fails, it defaults to BC_FARFIELD (0)
        
        # We expect:
        #   Tag 1 (Bottom) -> BC_WALL
        #   Tag 3 (Left)   -> BC_INLET
        
        # Verify Element 0
        self.assertEqual(bc_mask[0, 0], bc.BC_WALL, "Bottom face of Element 0 should be WALL")
        self.assertEqual(bc_mask[0, 3], bc.BC_INLET, "Left face of Element 0 should be INLET")
        
        # Element 3 (Top-Right) (ix=1, iy=1) -> ID = 1*2 + 1 = 3
        #   Face 1 (Right) -> Tag 4
        #   Face 2 (Top)   -> Tag 2
        
        self.assertEqual(bc_mask[3, 1], bc.BC_OUTLET, "Right face of Element 3 should be OUTLET")
        self.assertEqual(bc_mask[3, 2], bc.BC_WALL, "Top face of Element 3 should be WALL")

if __name__ == '__main__':
    unittest.main()
