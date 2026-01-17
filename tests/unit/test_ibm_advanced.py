import warp as wp
import numpy as np
import pytest
import os
from src.core.config import IBMConfig
from src.geometry.ibm import IBMManager
from src.kernels.ibm_kernels import generate_sdf_cylinder, generate_sdf_from_mesh

# Ensure Warp is initialized
wp.init()

class TestIBMAdvanced:
    def test_sdf_gradient_accuracy(self):
        """
        Verifies that the gradient of the SDF (phi) accurately recovers the surface normal.
        n = grad(phi). For a cylinder, normal is radial.
        """
        N = 20
        # Grid range [-2, 2]
        x_lin = np.linspace(-2.0, 2.0, N)
        y_lin = np.linspace(-2.0, 2.0, N)
        X, Y = np.meshgrid(x_lin, y_lin)
        
        # Flatten for Warp
        x_flat = X.flatten().astype(np.float32)
        y_flat = Y.flatten().astype(np.float32)
        num_points = x_flat.size
        
        # Reshape to (1, num_points) for 2D kernel compatibility
        x_wp = wp.array(x_flat.reshape(1, -1), dtype=float, device="cpu")
        y_wp = wp.array(y_flat.reshape(1, -1), dtype=float, device="cpu")
        phi_wp = wp.zeros((1, num_points), dtype=float, device="cpu")
        active_indices = wp.array([0], dtype=wp.int32, device="cpu")
        
        center_x, center_y, radius = 0.0, 0.0, 1.0
        
        # 1. Generate SDF
        wp.launch(
            kernel=generate_sdf_cylinder,
            dim=(1, num_points),
            inputs=[x_wp, y_wp, phi_wp, active_indices, center_x, center_y, radius],
            device="cpu"
        )
        
        phi = phi_wp.numpy()[0].reshape(N, N)
        dx = x_lin[1] - x_lin[0]
        dy = y_lin[1] - y_lin[0]
        
        # 2. Compute Numerical Gradient (Central Difference)
        grad = np.gradient(phi, dy, dx) 
        grad_x = grad[1]
        grad_y = grad[0]
        
        # 3. Compare with Analytical Normal
        total_error = 0.0
        count = 0
        
        for i in range(N):
            for j in range(N):
                px = X[i, j]
                py = Y[i, j]
                dist = np.sqrt(px**2 + py**2)
                
                if dist > 0.1 and 1 < i < N-2 and 1 < j < N-2:
                    nx_ana = px / dist
                    ny_ana = py / dist
                    nx_num = grad_x[i, j]
                    ny_num = grad_y[i, j]
                    
                    error = np.sqrt((nx_ana - nx_num)**2 + (ny_ana - ny_num)**2)
                    total_error += error
                    count += 1
        
        avg_error = total_error / count if count > 0 else 0.0
        print(f"Average Gradient Error: {avg_error}")
        assert avg_error < 2.0e-2

    def test_sdf_inversion_flag(self):
        """
        Tests that the 'invert' flag correctly flips the sign of the SDF.
        """
        points_data = np.array([
            [-1, -1, -1], [ 1, -1, -1], [ 1,  1, -1], [-1,  1, -1],
            [-1, -1,  1], [ 1, -1,  1], [ 1,  1,  1], [-1,  1,  1]
        ], dtype=np.float32)
        
        indices_data = np.array([
            0, 2, 1, 0, 3, 2,
            4, 5, 6, 4, 6, 7,
            0, 1, 5, 0, 5, 4,
            1, 2, 6, 1, 6, 5,
            2, 3, 7, 2, 7, 6,
            3, 0, 4, 3, 4, 7
        ], dtype=np.int32)
        
        points = wp.array(points_data, dtype=wp.vec3, device="cpu")
        indices = wp.array(indices_data, dtype=wp.int32, device="cpu")
        mesh = wp.Mesh(points, indices, points)
        
        # Query Points
        Np = 2
        x = wp.array([[0.0, 2.0]], dtype=float, device="cpu")
        y = wp.array([[0.0, 0.0]], dtype=float, device="cpu")
        phi_0 = wp.zeros((1, 2), dtype=float, device="cpu")
        phi_1 = wp.zeros((1, 2), dtype=float, device="cpu")
        active_indices = wp.array([0], dtype=wp.int32, device="cpu")
        
        # Run Invert=0
        wp.launch(
            kernel=generate_sdf_from_mesh,
            dim=(1, Np),
            inputs=[x, y, phi_0, active_indices, mesh.id, 10.0, 0],
            device="cpu"
        )
        
        # Run Invert=1
        wp.launch(
            kernel=generate_sdf_from_mesh,
            dim=(1, Np),
            inputs=[x, y, phi_1, active_indices, mesh.id, 10.0, 1],
            device="cpu"
        )
        
        p0 = phi_0.numpy()[0]
        p1 = phi_1.numpy()[0]
        
        assert np.allclose(p0, -p1)

    def test_mesh_boundary_robustness(self):
        """
        Tests SDF robustness when querying points exactly on vertices/edges.
        """
        points_data = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ], dtype=np.float32)
        indices_data = np.array([0, 1, 2], dtype=np.int32)
        
        points = wp.array(points_data, dtype=wp.vec3, device="cpu")
        indices = wp.array(indices_data, dtype=wp.int32, device="cpu")
        mesh = wp.Mesh(points, indices, points)
        
        # Query exactly on vertices
        Np = 3
        x = wp.array([[0.0, 1.0, 0.0]], dtype=float, device="cpu")
        y = wp.array([[0.0, 0.0, 1.0]], dtype=float, device="cpu")
        phi = wp.zeros((1, 3), dtype=float, device="cpu")
        active_indices = wp.array([0], dtype=wp.int32, device="cpu")
        
        wp.launch(
            kernel=generate_sdf_from_mesh,
            dim=(1, Np),
            inputs=[x, y, phi, active_indices, mesh.id, 10.0, 0],
            device="cpu"
        )
        
        res = phi.numpy()[0]
        assert np.allclose(res, 0.0, atol=1e-6)
        assert not np.any(np.isnan(res))

    def test_config_validation(self):
        """
        Verifies IBMManager raises FileNotFoundError for missing STL.
        """
        config = IBMConfig(
            enabled=True,
            mode="stl_file",
            stl_path="non_existent_ghost_file.stl"
        )
        
        with pytest.raises(FileNotFoundError):
            IBMManager(config, device="cpu")