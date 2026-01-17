import pytest
import numpy as np
import warp as wp
import meshio
from src.core.config import PazuzuConfig
from solver import PazuzuSolver

# Ensure Warp is initialized
wp.init()

class TestIBMIntegration:
    def test_analytical_cylinder_integration(self):
        """
        Verify that setting mode="analytical" correctly populates state.phi.
        """
        # 1. Config
        # Cylinder at (0,0), Radius 0.5
        # Domain [-1, 1] x [-1, 1]
        config_dict = {
            "case_name": "test_analytical_ibm",
            "io": {"write_interval": 10},
            "mesh": {"x_min": -1.0, "x_max": 1.0, "y_min": -1.0, "y_max": 1.0, "initial_depth": 2},
            "ibm": {
                "enabled": True,
                "mode": "analytical",
                "geometric_params": {"center_x": 0.0, "center_y": 0.0, "radius": 0.5}
            }
        }
        config = PazuzuConfig(**config_dict)
        
        # 2. Init Solver
        solver = PazuzuSolver(config)
        
        # 3. Verify Phi
        # We need to find points closest to (0,0) and (1,0)
        # Or just inspect the whole array.
        phi = solver.state.phi.numpy()
        x = solver.state.x.numpy()
        y = solver.state.y.numpy()
        
        # Check active blocks
        # num_blocks is likely an int in Quadtree implementation
        num_blocks = solver.quadtree.num_blocks
        if hasattr(num_blocks, "numpy"):
            num_blocks = int(num_blocks.numpy()[0])
            
        active_indices = solver.state.active_block_indices.numpy()[:num_blocks]
        
        found_inside = False
        found_outside = False
        
        for idx in active_indices:
            # Check nodes in this block
            for n in range(solver.basis.Np):
                px = x[idx, n]
                py = y[idx, n]
                val = phi[idx, n]
                
                dist = np.sqrt(px**2 + py**2)
                expected = dist - 0.5
                
                assert np.isclose(val, expected, atol=1e-5), f"Phi mismatch at ({px},{py})"
                
                if dist < 0.5:
                    found_inside = True
                else:
                    found_outside = True
                    
        assert found_inside, "No points found inside cylinder"
        assert found_outside, "No points found outside cylinder"

    def test_stl_file_integration_full_run(self, tmp_path):
        """
        Verify that the solver can load an STL, generate SDF, and run a time step.
        """
        # 1. Create Dummy STL (Cube/Square)
        stl_path = tmp_path / "dummy_cube.stl"
        points = np.array([
            [-0.5, -0.5, -0.5], [ 0.5, -0.5, -0.5], [ 0.5,  0.5, -0.5], [-0.5,  0.5, -0.5],
            [-0.5, -0.5,  0.5], [ 0.5, -0.5,  0.5], [ 0.5,  0.5,  0.5], [-0.5,  0.5,  0.5]
        ])
        # Simple cube connectivity
        # Ensure OUTWARD normals (CCW looking from outside)
        cells = [("triangle", np.array([
            [0, 2, 1], [0, 3, 2], # Bottom (Normal -Z)
            [4, 5, 6], [4, 6, 7], # Top (Normal +Z)
            [0, 1, 5], [0, 5, 4], # Front (Normal -Y)
            [1, 2, 6], [1, 6, 5], # Right
            [2, 3, 7], [2, 7, 6], # Back
            [3, 0, 4], [3, 4, 7]  # Left
        ]))]
        meshio.write(stl_path, meshio.Mesh(points, cells))
        
        # 2. Config
        config_dict = {
            "case_name": "test_stl_ibm",
            "io": {"write_interval": 1},
            "simulation": {"max_steps": 1, "t_final": 0.1},
            "mesh": {"x_min": -1.0, "x_max": 1.0, "y_min": -1.0, "y_max": 1.0, "initial_depth": 1},
            "ibm": {
                "enabled": True,
                "mode": "stl_file",
                "stl_path": str(stl_path)
            },
            # Initial condition: Uniform flow to test forcing?
            # Or vortex. Let's stick to default vortex.
            "initial_condition": {"name": "vortex"}
        }
        config = PazuzuConfig(**config_dict)
        
        # 3. Init Solver
        solver = PazuzuSolver(config)
        
        # Assertions Pre-Run
        assert solver.ibm.mesh is not None
        
        # Check phi
        phi = solver.state.phi.numpy()
        
        num_blocks = solver.quadtree.num_blocks
        if hasattr(num_blocks, "numpy"):
            num_blocks = int(num_blocks.numpy()[0])
            
        active_indices = solver.state.active_block_indices.numpy()[:num_blocks]
        
        # Check center (0,0) - should be inside cube (-0.5 to 0.5)
        # Dist to face is 0.5. Inside -> -0.5
        # We need to find the node closest to 0,0
        min_dist_origin = 1e9
        phi_at_origin = 0.0
        
        x = solver.state.x.numpy()
        y = solver.state.y.numpy()
        
        for idx in active_indices:
            for n in range(solver.basis.Np):
                px = x[idx, n]
                py = y[idx, n]
                d = np.sqrt(px**2 + py**2)
                if d < min_dist_origin:
                    min_dist_origin = d
                    phi_at_origin = phi[idx, n]
        
        print(f"Phi at approx origin (dist {min_dist_origin}): {phi_at_origin}")
        # Should be negative
        assert phi_at_origin < 0.0
        
        # 4. Action: Run
        print("Running solver step...")
        solver.run()
        
        # 5. Assertions Post-Run
        assert solver.state.step == 1
        print("Solver completed step 1 successfully.")
