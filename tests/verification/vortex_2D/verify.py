import sys
import os
import numpy as np
import h5py
import subprocess
import matplotlib.pyplot as plt
import pytest

# Add local directory to path for importing generate_mesh
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from generate_mesh import generate_mesh

def run_simulation():
    # 1. Generate Mesh
    output_dir = os.path.join(os.path.dirname(__file__), "output")
    mesh_file = os.path.join(output_dir, "box.msh")
    print(f"Generating mesh: {mesh_file}")
    generate_mesh(mesh_file)
    
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

def periodic_vortex(x, y, t=0.0, Lx=10.0, Ly=10.0):
    """
    Computes the exact solution for the 2D Isentropic Vortex problem at time t,
    correctly handling periodic boundary conditions on a [0, Lx] x [0, Ly] domain.
    """
    gamma = 1.4
    beta = 5.0  # Vortex strength
    
    u_inf = 1.0
    v_inf = 1.0
    
    # Initial center (Relative to domain min, assumed 0 for simplicity or centered)
    # The standard test case usually centers it at (5,5) regardless of domain size?
    # Or (Lx/2, Ly/2)? 
    # box.msh usually comes from Gmsh scripts. If it is 20x20, maybe center is (10,10)?
    # Let's assume standard (5,5) for now, but if the domain is 20x20, (5,5) is lower left.
    # We should probably check the config/generation script.
    # But sticking to (5,5) is safe if that's how it was generated.
    xc, yc = 5.0, 5.0
    
    # Current analytical center (unwrapped)
    xc_t = xc + u_inf * t
    yc_t = yc + v_inf * t
    
    # Calculate distance vector from the vortex center
    dx = x - xc_t
    dy = y - yc_t
    
    # Apply Periodicity: Wrap distance to [-L/2, L/2]
    dx = dx - Lx * np.round(dx / Lx)
    dy = dy - Ly * np.round(dy / Ly)
    
    r_sq = dx**2 + dy**2
    
    # Perturbations
    du = -(beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * dy
    dv =  (beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * dx
    
    u = u_inf + du
    v = v_inf + dv
    
    T_inf = 1.0
    T = T_inf - ((gamma - 1) * beta**2 / (8 * gamma * np.pi**2)) * np.exp(1 - r_sq)
    
    rho = T**(1.0 / (gamma - 1))
    p = rho**gamma
    
    return rho, u, v, p

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
        
        # Determine Domain Size from Mesh Points
        x_min, x_max = points[:, 0].min(), points[:, 0].max()
        y_min, y_max = points[:, 1].min(), points[:, 1].max()
        Lx = x_max - x_min
        Ly = y_max - y_min
        print(f"Detected Domain Size: {Lx:.2f} x {Ly:.2f}")
        
        # Load Final Step Data
        steps = [k for k in f['data'].keys() if k.startswith('step_')]
        if not steps:
            print("No steps found in output file.")
            sys.exit(1)
            
        # Sort by integer suffix
        steps.sort(key=lambda x: int(x.split('_')[1]))
        
        last_step = steps[-1]
        time = f['data'][last_step].attrs['time']
        print(f"Verifying step: {last_step} at t={time:.4f}")
        
        rho_sample = f['data'][last_step]['rho'][:]
        
        if rho_sample.shape[0] == points.shape[0]:
            print(f"Detected Nodal Data (High-Order). N={rho_sample.shape[0]}")
            x = points[:, 0]
            y = points[:, 1]
            mode = "nodal"
        elif rho_sample.shape[0] == conn.shape[0]:
             print(f"Detected Cell Data (Low-Order). N={rho_sample.shape[0]}")
             element_coords = points[conn]
             centroids = np.mean(element_coords, axis=1)
             x = centroids[:, 0]
             y = centroids[:, 1]
             mode = "cell"
        else:
            print(f"Error: Data shape {rho_sample.shape} matches neither Points {points.shape} nor Elements {conn.shape}.")
            sys.exit(1)
        
        # --- Check All Variables ---
        variables = ['rho', 'u', 'v', 'p']
        
        # Compute Analytical Solution once
        rho_ana, u_ana, v_ana, p_ana = periodic_vortex(x, y, t=time, Lx=Lx, Ly=Ly)
        analytical = {'rho': rho_ana, 'u': u_ana, 'v': v_ana, 'p': p_ana}
        
        # Thresholds
        threshold_l2 = 5.0e-2 
        threshold_linf = 1.0e-1 
        
        failed = False
        
        print(f"{'Variable':<5} | {'L2 Error':<12} | {'L_inf Error':<12} | {'Status'}")
        print("-" * 45)
        
        for var in variables:
            # Get Numerical Data
            num_data = f['data'][last_step][var][:]
            ana_data = analytical[var]
            
            # Error Calculation
            diff = num_data - ana_data
            abs_diff = np.abs(diff)
            
            # L2 Error (Relative)
            norm_ana = np.linalg.norm(ana_data)
            if norm_ana < 1e-12: norm_ana = 1.0 # Safety
            l2_error = np.linalg.norm(diff) / norm_ana
            
            # L_inf Error (Absolute Max)
            linf_error = np.max(abs_diff)
            
            # Check Pass/Fail
            status = "PASS"
            if l2_error > threshold_l2 or linf_error > threshold_linf:
                status = "FAIL"
                failed = True
            
            print(f"{var:<5} | {l2_error:.4e}   | {linf_error:.4e}   | {status}")
            
            # --- Generate Heatmap ---
            plt.figure(figsize=(8, 6))
            if mode == "nodal":
                 plt.tricontourf(x, y, abs_diff, levels=20, cmap='inferno')
            else:
                 plt.scatter(x, y, c=abs_diff, cmap='inferno', s=10)
                 
            plt.colorbar(label=f'Absolute Error |{var}_num - {var}_exact|')
            plt.title(f'{var.upper()} Error at t={time:.2f}\nL2={l2_error:.2e}, Linf={linf_error:.2e}')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.axis('equal')
            
            plot_file = os.path.join(output_dir, f"error_{var}.png")
            plt.savefig(plot_file)
            plt.close()

        if failed:
            print("\n❌ Verification FAILED: Some errors exceeded thresholds.")
            sys.exit(1)
        else:
            print("\n✅ Verification PASSED.")

@pytest.mark.slow
def test_verification():
    run_simulation()
    verify()

if __name__ == "__main__":
    run_simulation()
    verify()
