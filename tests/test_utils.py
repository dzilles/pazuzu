import numpy as np
from src.geometry.mesh import Mesh

def create_dummy_mesh(nx=2, ny=2, x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0, device="cpu"):
    """
    Creates a structured Cartesian Mesh for testing purposes.
    This replaces the internal generation logic previously in Mesh class.
    
    Args:
        nx (int): Number of elements in X.
        ny (int): Number of elements in Y.
        x_min, x_max, y_min, y_max: Domain bounds.
        device (str): Compute device.
        
    Returns:
        Mesh: An initialized Mesh object.
    """
    num_elements = nx * ny
    dx = (x_max - x_min) / nx
    dy = (y_max - y_min) / ny

    vertices = np.zeros((num_elements, 4, 2), dtype=np.float32)
    connectivity = -np.ones((num_elements, 4, 2), dtype=np.int32)
    boundary_tags = np.zeros((num_elements, 4), dtype=np.int32)

    # 1. Vertices (Counter-Clockwise)
    for j in range(ny):
        for i in range(nx):
            element_id = j * nx + i
            x0 = x_min + i * dx
            y0 = y_min + j * dy
            
            vertices[element_id, 0, :] = [x0, y0]
            vertices[element_id, 1, :] = [x0 + dx, y0]
            vertices[element_id, 2, :] = [x0 + dx, y0 + dy]
            vertices[element_id, 3, :] = [x0, y0 + dy]

    # 2. Connectivity and Tags
    # Faces: 0: Bottom, 1: Right, 2: Top, 3: Left
    for j in range(ny):
        for i in range(nx):
            element_id = j * nx + i

            # Left (Face 3)
            if i > 0:
                ni, nj = i - 1, j
                connectivity[element_id, 3] = [nj * nx + ni, 1]
            else:
                boundary_tags[element_id, 3] = 3 # Left Tag (Inlet)

            # Right (Face 1)
            if i < nx - 1:
                ni, nj = i + 1, j
                connectivity[element_id, 1] = [nj * nx + ni, 3]
            else:
                boundary_tags[element_id, 1] = 4 # Right Tag (Outlet)

            # Bottom (Face 0)
            if j > 0:
                ni, nj = i, j - 1
                connectivity[element_id, 0] = [nj * nx + ni, 2]
            else:
                boundary_tags[element_id, 0] = 1 # Bottom Tag

            # Top (Face 2)
            if j < ny - 1:
                ni, nj = i, j + 1
                connectivity[element_id, 2] = [nj * nx + ni, 0]
            else:
                boundary_tags[element_id, 2] = 2 # Top Tag

    physical_groups = {
        "Bottom": 1,
        "Top": 2,
        "Left": 3,
        "Right": 4
    }

    # Initialize the Mesh object with explicit data
    mesh = Mesh(
        vertices=vertices,
        connectivity=connectivity,
        boundary_tags=boundary_tags,
        physical_groups=physical_groups,
        device=device
    )
    
    # Manually attach attributes that were previously set for Cartesian meshes
    mesh.nx = nx
    mesh.ny = ny
    mesh.dx = dx
    mesh.dy = dy
    mesh.x_min = x_min
    mesh.x_max = x_max
    mesh.y_min = y_min
    mesh.y_max = y_max
    
    return mesh
