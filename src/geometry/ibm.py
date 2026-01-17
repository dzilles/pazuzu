import os
from typing import Optional

import meshio
import numpy as np
import warp as wp

from src.core.config import IBMConfig


class IBMManager:
    """
    Manages geometry for the Immersed Boundary Method (IBM).
    Handles loading of STL files and construction of acceleration structures (BVH) via Warp.
    """

    def __init__(self, config: IBMConfig, device: str):
        self.config = config
        self.device = device
        self.mesh: Optional[wp.Mesh] = None

        if self.config.enabled:
            self.load_geometry()

    def load_geometry(self):
        """
        Loads the geometry based on the configuration mode.
        """
        if self.config.mode == "stl_file":
            self._load_stl()
        elif self.config.mode == "analytical":
            # Placeholder for analytical shapes (e.g., cylinder, sphere)
            # These will eventually be handled by specific SDF kernels
            print(f"IBMManager: Analytical mode selected. Params: {self.config.geometric_params}")
        else:
            raise ValueError(f"Unknown IBM mode: {self.config.mode}")

    def _load_stl(self):
        """
        Loads an STL file using meshio, extracts vertices and faces,
        and builds a Warp Mesh object (BVH).
        """
        stl_path = self.config.stl_path
        if not stl_path:
            raise ValueError("IBM mode is 'stl_file' but 'stl_path' is not provided in config.")

        if not os.path.exists(stl_path):
            raise FileNotFoundError(f"STL file not found: {stl_path}")

        print(f"IBMManager: Loading STL geometry from {stl_path}...")

        try:
            # Load mesh using meshio
            mesh = meshio.read(stl_path)
        except Exception as e:
            raise RuntimeError(f"Failed to read STL file with meshio: {e}")

        # Extract vertices (ensure they are float32 for Warp)
        # mesh.points is typically (N, 3)
        vertices_np = mesh.points.astype(np.float32)

        # Extract faces (triangles)
        faces_np = None
        for cell_block in mesh.cells:
            if cell_block.type == "triangle":
                # cell_block.data is (num_triangles, 3)
                faces_np = cell_block.data.astype(np.int32)
                break
        
        if faces_np is None:
            raise ValueError(f"No triangular faces found in {stl_path}. Ensure it is a valid binary/ascii STL.")

        # Warp Mesh expects flattened indices for triangle lists: (v0, v1, v2, v0, v1, v2, ...)
        indices_flat = faces_np.flatten()

        print(f"IBMManager: Constructing Warp Mesh (BVH) with {len(vertices_np)} vertices and {len(faces_np)} faces...")

        # Create Warp arrays on the specified device
        # wp.vec3 corresponds to float32[3]
        wp_points = wp.array(vertices_np, dtype=wp.vec3, device=self.device)
        wp_indices = wp.array(indices_flat, dtype=wp.int32, device=self.device)

        # Create the Warp Mesh
        # This automatically builds the BVH for spatial queries
        self.mesh = wp.Mesh(
            points=wp_points,
            indices=wp_indices
        )

        print("IBMManager: Geometry loaded successfully.")
