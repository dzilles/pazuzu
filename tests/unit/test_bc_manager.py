import unittest
import numpy as np
from unittest.mock import MagicMock
from src.core.boundary_condition_manager import BoundaryConditionManager
from src.geometry.mesh import Mesh
from src.kernels import boundary_conditions as bc

class MockPhysics:
    def __init__(self):
        self.rho_inf = 1.0
        self.u_inf = 0.0
        self.v_inf = 0.0
        self.p_inf = 1.0

class MockConfig:
    def __init__(self, boundaries):
        self.boundaries = boundaries
        self.physics = MockPhysics()

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
            "Left": {
                "type": "inlet",
                "params": {"rho": 1.0, "u": 1.0, "v": 0.0, "p": 1.0}
            },
            "Right": {
                "type": "outlet",
                "params": {"p_back": 1.0}
            }
        }
        config = MockConfig(boundaries)
        
        manager = BoundaryConditionManager(mesh, config)
        bc_mask, bc_data = manager.setup_boundary_conditions()
        
        # Check Element 0 (Bottom-Left)
        # Face 0 (Bottom) -> Tag 1 -> slip_wall
        # Face 3 (Left)   -> Tag 3 -> inlet
        
        idx_bottom = bc_mask[0, 0]
        idx_left = bc_mask[0, 3]
        
        self.assertGreaterEqual(idx_bottom, 0)
        self.assertGreaterEqual(idx_left, 0)
        
        self.assertEqual(bc_data[idx_bottom]['type'], bc.BC_WALL)
        self.assertEqual(bc_data[idx_left]['type'], bc.BC_INLET)
        
        # Element 3 (Top-Right)
        # Face 1 (Right) -> Tag 4 -> outlet
        # Face 2 (Top)   -> Tag 2 -> slip_wall
        
        idx_right = bc_mask[3, 1]
        idx_top = bc_mask[3, 2]
        
        self.assertEqual(bc_data[idx_right]['type'], bc.BC_OUTLET)
        self.assertEqual(bc_data[idx_top]['type'], bc.BC_WALL)

if __name__ == '__main__':
    unittest.main()
