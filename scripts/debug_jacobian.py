import sys
import os
import numpy as np
import warp as wp

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../")))

from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.physics.euler_2d import Euler2DSolver

def check_jacobians(mesh_file):
    print(f"Checking Jacobians for {mesh_file}...")
    
    # Initialize Warp (CPU is enough for checking)
    wp.init()
    
    # Load Mesh
    try:
        mesh = Mesh(filename=mesh_file, device="cpu")
    except Exception as e:
        print(f"Failed to load mesh: {e}")
        return

    # Basis (Degree 1 is sufficient for geometry check)
    basis = Basis(polynomial_degree=1, device="cpu")
    
    # Initialize Solver (this computes metrics)
    # We pass an empty config as we don't need full simulation setup
    solver = Euler2DSolver(mesh, basis, config={})
    
    # Check J
    J = solver.mesh.J.numpy()
    
    min_J = np.min(J)
    max_J = np.max(J)
    avg_J = np.mean(J)
    
    print(f"Jacobian Stats:")
    print(f"  Min J: {min_J}")
    print(f"  Max J: {max_J}")
    print(f"  Avg J: {avg_J}")
    
    negative_indices = np.where(J <= 0)[0]
    num_negative = len(negative_indices)
    
    if num_negative > 0:
        print(f"!!! FOUND {num_negative} ELEMENTS WITH NON-POSITIVE JACOBIAN !!!")
        print("First 5 bad elements:")
        for idx in negative_indices[:5]:
            print(f"  Elem {idx}: J={J[idx]}")
            verts = mesh.vertices.numpy()[idx]
            print(f"    Vertices: {verts}")
    else:
        print("All Jacobians are positive. Mesh orientation seems OK.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        mesh_file = sys.argv[1]
    else:
        mesh_file = "tests/verification/cylinder_flow/cylinder.msh"
        
    check_jacobians(mesh_file)
