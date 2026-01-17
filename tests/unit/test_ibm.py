import warp as wp
import numpy as np
import pytest
import os
from src.core.config import PazuzuConfig, IBMConfig
from src.core.simulation_state import SimulationState
from src.kernels.ibm_kernels import generate_sdf_cylinder, generate_sdf_from_mesh
from src.geometry.ibm import IBMManager

# Ensure Warp is initialized
wp.init()

class TestIBMConfig:
    def test_defaults(self):
        """Test that IBM is disabled by default and has correct default params."""
        config = PazuzuConfig(
            io={"output_dir": "test", "write_interval": 10}
        )
        assert config.ibm.enabled is False
        assert config.ibm.mode == "analytical"
        assert config.ibm.stl_path is None

    def test_custom_config(self):
        """Test setting IBM parameters."""
        ibm_conf = IBMConfig(
            enabled=True,
            mode="stl_file",
            stl_path="assets/test.stl",
            invert_inside_outside=True
        )
        config = PazuzuConfig(
            io={"output_dir": "test", "write_interval": 10},
            ibm=ibm_conf
        )
        assert config.ibm.enabled is True
        assert config.ibm.stl_path == "assets/test.stl"

class TestIBMState:
    def test_phi_allocation(self):
        """Verify SimulationState allocates the phi buffer."""
        # Setup minimal state
        Np = 4
        max_blocks = 2
        
        state = SimulationState(
            Np=Np,
            dtype=wp.vec4,
            device="cpu",
            max_blocks=max_blocks,
            scalar_dtype=wp.float32
        )
        
        # Check phi existence and shape
        assert hasattr(state, 'phi'), "SimulationState missing 'phi' attribute"
        assert state.phi.shape == (max_blocks, Np)
        assert state.phi.dtype == wp.float32
        
        # Check initialization to zero
        phi_host = state.phi.numpy()
        assert np.allclose(phi_host, 0.0)

class TestIBMKernels:
    def test_sdf_analytical_cylinder(self):
        """Test the analytical cylinder SDF kernel (2D array version)."""
        # 1. Setup Grid Points (3 points: Inside, On Surface, Outside)
        # Cylinder at (0,0) radius 1.0
        # P1: (0,0) -> dist = -1.0
        # P2: (1,0) -> dist = 0.0
        # P3: (2,0) -> dist = 1.0
        
        Np = 3
        # Reshape to (1, Np) for 2D kernel compatibility
        x = wp.array([[0.0, 1.0, 2.0]], dtype=wp.float32, device="cpu")
        y = wp.array([[0.0, 0.0, 0.0]], dtype=wp.float32, device="cpu")
        phi = wp.zeros((1, 3), dtype=wp.float32, device="cpu")
        active_indices = wp.array([0], dtype=wp.int32, device="cpu")
        
        center_x = 0.0
        center_y = 0.0
        radius = 1.0
        
        # 2. Launch Kernel (2D launch: 1 block, Np nodes)
        wp.launch(
            kernel=generate_sdf_cylinder,
            dim=(1, Np),
            inputs=[x, y, phi, active_indices, center_x, center_y, radius],
            device="cpu"
        )
        
        # 3. Verify
        res = phi.numpy()[0]
        assert np.isclose(res[0], -1.0)
        assert np.isclose(res[1], 0.0)
        assert np.isclose(res[2], 1.0)

    def test_sdf_from_mesh_cube(self):
        """
        Tests the Warp Mesh SDF query using a closed Cube mesh (2D array version).
        Cube from [-1,-1,-1] to [1,1,1].
        """
        # 1. Construct a Cube Mesh manually (12 triangles)
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
        
        # 2. Define Query Points (Z=0 plane)
        Np = 2
        x = wp.array([[0.0, 2.0]], dtype=wp.float32, device="cpu")
        y = wp.array([[0.0, 0.0]], dtype=wp.float32, device="cpu")
        phi = wp.zeros((1, 2), dtype=wp.float32, device="cpu")
        active_indices = wp.array([0], dtype=wp.int32, device="cpu")
        
        # 3. Launch Kernel
        wp.launch(
            kernel=generate_sdf_from_mesh,
            dim=(1, Np),
            inputs=[x, y, phi, active_indices, mesh.id, 10.0, 0],
            device="cpu"
        )
        
        # 4. Verify
        res = phi.numpy()[0]
        print(f"Cube Mesh SDF Results: {res}")
        
        assert np.isclose(res[0], -1.0), "Point (0,0) should be inside cube (dist -1.0)"
        assert np.isclose(res[1], 1.0), "Point (2,0) should be outside cube (dist 1.0)"

class TestIBMManager:
    def test_manager_analytical(self, tmp_path):
        """Test initializing manager in analytical mode."""
        # Create dummy config
        config_data = IBMConfig(enabled=True, mode="analytical")
        
        # Init Manager
        manager = IBMManager(config_data, device="cpu")
        
        # Ensure no mesh is loaded
        assert manager.mesh is None

    def test_manager_stl_loading(self, tmp_path):
        """
        Test full pipeline: Write dummy STL -> IBMManager load -> Verify Mesh.
        """
        import meshio
        
        # 1. Create a dummy STL file using meshio
        stl_path = tmp_path / "dummy.stl"
        
        points = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0]
        ])
        cells = [("triangle", np.array([[0, 1, 2]]))]
        
        meshio.write(stl_path, meshio.Mesh(points, cells))
        
        # 2. Config pointing to this file
        config_data = IBMConfig(
            enabled=True, 
            mode="stl_file", 
            stl_path=str(stl_path)
        )
        
        # 3. Init Manager
        manager = IBMManager(config_data, device="cpu")
        
        # 4. Verify
        assert manager.mesh is not None
        assert isinstance(manager.mesh, wp.Mesh)