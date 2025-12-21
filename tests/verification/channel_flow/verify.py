import sys
import os
import numpy as np
import h5py
import subprocess

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../..")))

def run_simulation():
    config_path = os.path.join(os.path.dirname(__file__), "channel.yaml")
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
        # Load Final Step Data
        steps = sorted([k for k in f['data'].keys() if k.startswith('step_')])
        last_step = steps[-1]
        print(f"Verifying step: {last_step}")
        
        rho = f['data'][last_step]['rho'][:]
        u = f['data'][last_step]['u'][:]
        v = f['data'][last_step]['v'][:]
        p = f['data'][last_step]['p'][:]
        
        # Check against uniform flow (1, 1, 0, 1)
        # Numerical error accumulates but should be small
        
        err_rho = np.max(np.abs(rho - 1.0))
        err_u = np.max(np.abs(u - 1.0))
        err_v = np.max(np.abs(v - 0.0))
        err_p = np.max(np.abs(p - 1.0))
        
        print(f"Max Error Rho: {err_rho:.4e}")
        print(f"Max Error U:   {err_u:.4e}")
        print(f"Max Error V:   {err_v:.4e}")
        print(f"Max Error P:   {err_p:.4e}")
        
        threshold = 1e-3
        if err_rho > threshold or err_u > threshold or err_v > threshold or err_p > threshold:
             print("Test FAILED: Solution deviated from uniform flow.")
             sys.exit(1)
        else:
             print("Test PASSED.")

if __name__ == "__main__":
    run_simulation()
    verify()
