import sys
import os
import numpy as np
import h5py
import subprocess
import matplotlib.pyplot as plt
import pytest

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))

def run_simulation():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "output")
    mesh_file = os.path.join(output_dir, "pulse.msh")
    mesh_script = os.path.join(script_dir, "generate_mesh.py")
    
    print(f"Generating mesh: {mesh_file}")
    try:
        subprocess.run([sys.executable, mesh_script], check=True, cwd=script_dir)
    except subprocess.CalledProcessError as e:
        print(f"Mesh generation failed: {e}")
        sys.exit(1)
    
    config_path = os.path.join(script_dir, "pulse.yaml")
    root_dir = os.path.abspath(os.path.join(script_dir, "../../../"))
    
    cmd = [sys.executable, "run_simulation.py", config_path]
    print(f"Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, cwd=root_dir, check=True, timeout=600)
        print("Simulation finished successfully.")
    except subprocess.TimeoutExpired:
        print("❌ Simulation timed out after 10 minutes!")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Simulation failed: {e}")
        sys.exit(1)

def verify():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "output")
    h5_file = os.path.join(output_dir, "results.h5")
    
    if not os.path.exists(h5_file):
        print(f"Error: {h5_file} not found.")
        sys.exit(1)
        
    with h5py.File(h5_file, 'r') as f:
        # Load Mesh
        points = f['mesh/points'][:]
        conn = f['mesh/connectivity'][:]
        
        # Load Final Step Data
        steps = [k for k in f['data'].keys() if k.startswith('step_')]
        if not steps:
            print("No steps found in output file.")
            sys.exit(1)
            
        steps.sort(key=lambda x: int(x.split('_')[1]))
        last_step = steps[-1]
        time = f['data'][last_step].attrs['time']
        print(f"Verifying step: {last_step} at t={time:.4f}")
        
        p = f['data'][last_step]['p'][:]
        
        # Determine coordinates for plotting
        if p.shape[0] == points.shape[0]:
            x = points[:, 0]
            y = points[:, 1]
            mode = "nodal"
        else:
            element_coords = points[conn]
            centroids = np.mean(element_coords, axis=1)
            x = centroids[:, 0]
            y = centroids[:, 1]
            mode = "cell"
            
        # Target pressure
        p_inf = 1.0
        diff = p - p_inf
        max_abs_error = np.max(np.abs(diff))
        l2_error = np.sqrt(np.mean(diff**2))
        
        print(f"Max Absolute Error in Pressure: {max_abs_error:.6e}")
        print(f"L2 Error in Pressure:           {l2_error:.6e}")
        
        # Pass Condition
        threshold = 1.0e-1
        if max_abs_error < threshold:
            print(f"✅ PASSED: Error is below threshold {threshold}")
        else:
            print(f"❌ FAILED: Error exceeds threshold {threshold}")
            
        # --- Visualization ---
        plt.figure(figsize=(8, 6))
        if mode == "nodal":
            plt.tricontourf(x, y, p, levels=50, cmap='viridis')
        else:
            plt.scatter(x, y, c=p, cmap='viridis', s=10)
            
        plt.colorbar(label='Pressure')
        plt.title(f'Final Pressure Field at t={time:.2f}\nMax Abs Error = {max_abs_error:.2e}')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        
        plot_file = os.path.join(output_dir, "final_pressure.png")
        plt.savefig(plot_file)
        print(f"Plot saved to {plot_file}")
        plt.close()
        
        if max_abs_error >= threshold:
            sys.exit(1)

@pytest.mark.slow
def test_verification():
    run_simulation()
    verify()

if __name__ == "__main__":
    run_simulation()
    verify()
