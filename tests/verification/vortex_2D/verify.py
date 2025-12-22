import sys
import os
import numpy as np
import h5py
import subprocess

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))

from src.physics.initial_conditions import vortex

def run_simulation():
    config_path = os.path.join(os.path.dirname(__file__), "vortex.yaml")
    # Run from root
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
    
    cmd = [sys.executable, "run_simulation.py", config_path]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=root_dir, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("Simulation failed!")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        sys.exit(1)
    else:
        print("Simulation finished successfully.")

def verify():
    output_dir = os.path.join(os.path.dirname(__file__), "output")
    h5_file = os.path.join(output_dir, "results.h5")
    
    if not os.path.exists(h5_file):
        print(f"Error: {h5_file} not found.")
        sys.exit(1)
        
    with h5py.File(h5_file, 'r') as f:
        # Load Mesh
        points = f['mesh/points'][:]
        conn = f['mesh/connectivity'][:]
        
        # Calculate Centroids
        # conn is (N_elems, 4)
        # points is (N_points, 3)
        
        # Get coordinates of all 4 nodes for all elements
        # Shape (N_elems, 4, 3)
        element_coords = points[conn] 
        centroids = np.mean(element_coords, axis=1) # (N_elems, 3)
        
        # Load Final Step Data
        # Find last step (numerical sort)
        steps = [k for k in f['data'].keys() if k.startswith('step_')]
        if not steps:
            print("No steps found in output file.")
            sys.exit(1)
            
        # Sort by integer suffix
        steps.sort(key=lambda x: int(x.split('_')[1]))
        
        last_step = steps[-1]
        time = f['data'][last_step].attrs['time']
        print(f"Verifying step: {last_step} at t={time:.4f}")
        
        rho_num = f['data'][last_step]['rho'][:]
        
        # Calculate Analytical Solution at Centroids
        # Centroids are in 3D (z=0), we take x,y
        x = centroids[:, 0]
        y = centroids[:, 1]
        
        # Vectorized call to vortex
        rho_ana, _, _, _ = vortex(x, y, t=time)
        
        # Error (L2 norm relative)
        diff = rho_num - rho_ana
        l2_error = np.linalg.norm(diff) / np.linalg.norm(rho_ana)
        
        print(f"L2 Error (Density): {l2_error:.4e}")
        
        # Threshold
        # Since we use cell averages on a 20x20 mesh compared to point-wise analytical,
        # the error might be dominated by the averaging/sampling error rather than solver error.
        # But it should be reasonably small. 
        threshold = 5.0e-2 # 5% tolerance
        if l2_error > threshold:
            print(f"Test FAILED: Error {l2_error:.4e} > {threshold}")
            sys.exit(1)
        else:
            print("Test PASSED.")

if __name__ == "__main__":
    run_simulation()
    verify()
