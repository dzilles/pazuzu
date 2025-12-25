import os
import sys
import subprocess
import numpy as np
import h5py
import pytest

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))

def run_test():
    # Setup paths
    test_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(test_dir, "../../../"))
    
    mesh_script = os.path.join(test_dir, "generate_channel.py")
    config_file = os.path.join(test_dir, "channel.yaml")
    run_script = os.path.join(project_root, "run_simulation.py")
    
    output_dir = os.path.join(test_dir, "output")
    results_file = os.path.join(output_dir, "results.h5")

    # 1. Generate Mesh
    print("Generating mesh...")
    try:
        subprocess.run([sys.executable, mesh_script], check=True, cwd=test_dir)
    except subprocess.CalledProcessError as e:
        print(f"Mesh generation failed: {e}")
        sys.exit(1)

    # 2. Run Simulation
    print(f"Running simulation with {config_file}...")
    try:
        # Check if config exists
        if not os.path.exists(config_file):
            print(f"Config file not found: {config_file}")
            sys.exit(1)

        subprocess.run([sys.executable, run_script, config_file], check=True, cwd=project_root)
    except subprocess.CalledProcessError as e:
        print(f"Simulation failed: {e}")
        sys.exit(1)

    # 3. Verify Results
    print("Verifying results...")
    if not os.path.exists(results_file):
        print(f"Results file not found: {results_file}")
        sys.exit(1)

    try:
        with h5py.File(results_file, 'r') as f:
            # Check if any steps were recorded
            if 'data' not in f:
                print("No 'data' group found in HDF5 file.")
                sys.exit(1)
                
            data_grp = f['data']
            keys = list(data_grp.keys())
            # print(f"HDF5 Step Keys: {keys}")
            
            step_keys = [k for k in keys if k.startswith("step_")]
            
            if not step_keys:
                print("No 'step_XX' groups found in data group.")
                sys.exit(1)
            
            # Sort by step number
            step_keys.sort(key=lambda k: int(k.split('_')[1]))
            
            # Get the last step
            last_step = step_keys[-1]
            print(f"Checking last step: {last_step}")
            
            # Load data: group contains 'rho', 'u', 'v', 'p' (Cell Averages)
            grp = data_grp[last_step]
            
            rho = grp['rho'][:]
            u = grp['u'][:]
            v = grp['v'][:]
            p = grp['p'][:]
            
            # --- Uniform Flow Validation ---
            # Check against uniform flow (1, 1, 0, 1)
            # Numerical error accumulates but should be small
            
            err_rho = np.max(np.abs(rho - 1.0))
            err_u = np.max(np.abs(u - 1.0))
            err_v = np.max(np.abs(v - 0.0))
            err_p = np.max(np.abs(p - 1.0))
            
            print(f"  Max Error Rho: {err_rho:.4e}")
            print(f"  Max Error U:   {err_u:.4e}")
            print(f"  Max Error V:   {err_v:.4e}")
            print(f"  Max Error P:   {err_p:.4e}")
            
            # Threshold for passing (considering 3rd order and time integration)
            threshold = 1e-2 
            if err_rho > threshold or err_u > threshold or err_v > threshold or err_p > threshold:
                 print("FAIL: Solution deviated significantly from uniform flow.")
                 sys.exit(1)
            else:
                 print("Test PASSED.")

    except Exception as e:
        print(f"Verification failed with error: {e}")
        sys.exit(1)

@pytest.mark.slow
def test_verification():
    run_test()

if __name__ == "__main__":
    run_test()