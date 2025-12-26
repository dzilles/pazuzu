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
    
    mesh_script = os.path.join(test_dir, "generate_mesh.py")
    config_file = os.path.join(test_dir, "cylinder.yaml")
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

    # 2. Run Simulation (Short run for verification)
    print(f"Running simulation with {config_file} (limited to 50 steps)...")
    try:
        # Override max_steps temporarily for verification
        temp_config = os.path.join(test_dir, "temp_verify.yaml")
        import yaml
        with open(config_file, 'r') as f:
            cfg = yaml.safe_load(f)
        cfg['simulation']['max_steps'] = 50
        with open(temp_config, 'w') as f:
            yaml.dump(cfg, f)

        subprocess.run([sys.executable, run_script, temp_config], check=True, cwd=project_root)
        os.remove(temp_config)
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
            if 'data' not in f:
                print("No 'data' group found in HDF5 file.")
                sys.exit(1)
                
            data_grp = f['data']
            step_keys = [k for k in data_grp.keys() if k.startswith("step_")]
            if not step_keys:
                print("No 'step_XX' groups found.")
                sys.exit(1)
            
            step_keys.sort(key=lambda k: int(k.split('_')[1]))
            last_step = step_keys[-1]
            print(f"Checking step: {last_step}")
            
            grp = data_grp[last_step]
            rho = grp['rho'][:]
            p = grp['p'][:]
            
            if np.isnan(rho).any() or (rho <= 0).any():
                print("FAIL: Density invalid.")
                sys.exit(1)
            if np.isnan(p).any() or (p <= 0).any():
                print("FAIL: Pressure invalid.")
                sys.exit(1)

            print("Navier-Stokes Cylinder Test PASSED.")

    except Exception as e:
        print(f"Verification failed: {e}")
        sys.exit(1)

@pytest.mark.slow
def test_verification_ns():
    run_test()

if __name__ == "__main__":
    run_test()
