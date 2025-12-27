import sys
import os
import numpy as np
import h5py
import subprocess
import matplotlib.pyplot as plt
import yaml
import argparse
from scipy import stats

# Add local directory and original vortex directory to path
test_dir = os.path.dirname(os.path.abspath(__file__))
vortex_dir = os.path.abspath(os.path.join(test_dir, "../vortex_2D"))
sys.path.append(test_dir)
sys.path.append(vortex_dir)

# Import analytical solution from verify.py in vortex_2D
try:
    from verify import periodic_vortex
except ImportError:
    def periodic_vortex(x, y, t=0.0, Lx=10.0, Ly=10.0):
        # ... (Same fallback as before)
        gamma = 1.4
        beta = 5.0  # Vortex strength
        u_inf = 1.0
        v_inf = 1.0
        xc, yc = 5.0, 5.0
        xc_t = xc + u_inf * t
        yc_t = yc + v_inf * t
        dx = x - xc_t
        dy = y - yc_t
        dx = dx - Lx * np.round(dx / Lx)
        dy = dy - Ly * np.round(dy / Ly)
        r_sq = dx**2 + dy**2
        du = -(beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * dy
        dv =  (beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * dx
        u = u_inf + du
        v = v_inf + dv
        T_inf = 1.0
        T = T_inf - ((gamma - 1) * beta**2 / (8 * gamma * np.pi**2)) * np.exp(1 - r_sq)
        rho = T**(1.0 / (gamma - 1))
        p = rho**gamma
        return rho, u, v, p

def run_study():
    resolutions = [40, 60, 80, 100]
    orders = [1, 2, 3]
    
    all_results = {}
    
    root_dir = os.path.abspath(os.path.join(test_dir, "../../.."))
    base_config_path = os.path.join(test_dir, "vortex.yaml")
    output_dir = os.path.join(test_dir, "output")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with open(base_config_path, 'r') as f:
        base_config = yaml.safe_load(f)

    # Use a shorter t_final for faster convergence study
    t_final = 1.0 
    base_config['simulation']['t_final'] = t_final
    base_config['io']['write_interval'] = 10000 # Only write final state
    
    for P in orders:
        print(f"\n========================================")
        print(f"  Studying Polynomial Order P = {P}")
        print(f"========================================")
        
        errors = []
        h_values = []
        
        for N in resolutions:
            print(f"\n--- Running P={P}, N={N} ---")
            
            # 1. Generate Mesh
            mesh_file = os.path.join(output_dir, f"box_{N}.msh")
            mesh_script = os.path.join(test_dir, "generate_mesh.py")
            subprocess.run([sys.executable, mesh_script, "--res", str(N), "--out", mesh_file], check=True)
            
            # 2. Update Config
            run_output_dir = os.path.join(output_dir, f"run_P{P}_N{N}")
            if not os.path.exists(run_output_dir):
                os.makedirs(run_output_dir)
                
            config = base_config.copy()
            config['numerics'] = base_config['numerics'].copy()
            config['io'] = base_config['io'].copy()
            config['simulation'] = base_config['simulation'].copy()
            
            config['mesh_file'] = mesh_file
            config['numerics']['polynomial_order'] = P
            config['numerics']['over_integration_order'] = 0
            config['io']['output_dir'] = run_output_dir
            
            temp_config_path = os.path.join(test_dir, f"vortex_P{P}_N{N}.yaml")
            with open(temp_config_path, 'w') as f:
                yaml.dump(config, f)
                
            # 3. Run Simulation
            cmd = [sys.executable, "run_simulation.py", temp_config_path]
            subprocess.run(cmd, cwd=root_dir, check=True)
            
            # 4. Analyze Results
            h5_file = os.path.join(run_output_dir, "results.h5")
            with h5py.File(h5_file, 'r') as f:
                points = f['mesh/points'][:]
                steps = sorted([k for k in f['data'].keys() if k.startswith('step_')], 
                               key=lambda x: int(x.split('_')[1]))
                last_step = steps[-1]
                time = f['data'][last_step].attrs['time']
                
                rho_num = f['data'][last_step]['rho'][:]
                x = points[:, 0]
                y = points[:, 1]
                
                # Domain size for periodic wrapping
                x_min, x_max = points[:, 0].min(), points[:, 0].max()
                y_min, y_max = points[:, 1].min(), points[:, 1].max()
                Lx, Ly = x_max - x_min, y_max - y_min
                
                rho_ana, _, _, _ = periodic_vortex(x, y, t=time, Lx=Lx, Ly=Ly)
                
                diff = rho_num - rho_ana
                # Discrete L2 norm
                l2_error = np.sqrt(np.mean(diff**2))
                
                errors.append(l2_error)
                h_val = Lx / N
                h_values.append(h_val)

                # --- Print Error Comparison ---
                expected_l2 = 0.0
                if len(errors) == 1:
                    print(f"L2 Error: {l2_error:.6e} (Baseline)")
                else:
                    # Expected = Baseline_Error * (h / Baseline_h)^(P+1)
                    expected_l2 = errors[0] * (h_val / h_values[0])**(P+1)
                    ratio = l2_error / expected_l2
                    print(f"L2 Error: {l2_error:.6e} | Expected: {expected_l2:.6e} | Actual/Expected: {ratio:.2f}")
        
        all_results[P] = {
            'h': np.array(h_values),
            'errors': np.array(errors)
        }

        # --- Create Plot for this Order P ---
        h = all_results[P]['h']
        err = all_results[P]['errors']
        
        log_h = np.log(h)
        log_err = np.log(err)
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_h, log_err)
        
        print(f"P = {P}: Order of Convergence (Slope) = {slope:.3f}")

        plt.figure(figsize=(8, 6))
        plt.loglog(h, err, 'o-', label=f'Numerical P={P} (Slope={slope:.2f})')
        
        # Reference line O(h^{P+1})
        expected_slope = P + 1
        h_ref = np.linspace(min(h), max(h), 100)
        err_ref = err[0] * (h_ref / h[0])**expected_slope
        plt.loglog(h_ref, err_ref, 'k--', alpha=0.5, label=f'Reference O(h^{expected_slope})')
        
        plt.xlabel('Mesh size h')
        plt.ylabel('L2 Error (Density)')
        plt.title(f'Grid Convergence Study: P={P}')
        plt.legend()
        plt.grid(True, which="both", ls="-", alpha=0.5)
        
        plot_path = os.path.join(output_dir, f"convergence_P{P}.png")
        plt.savefig(plot_path)
        plt.close()
        print(f"Plot saved to {plot_path}")

    # Final Summary Plot (All in one)
    plt.figure(figsize=(10, 8))
    for P in orders:
        h = all_results[P]['h']
        err = all_results[P]['errors']
        log_h = np.log(h)
        log_err = np.log(err)
        slope, _, _, _, _ = stats.linregress(log_h, log_err)
        plt.loglog(h, err, 'o-', label=f'P={P} (Slope={slope:.2f})')

    plt.xlabel('Mesh size h')
    plt.ylabel('L2 Error (Density)')
    plt.title('P-Refinement Convergence Study: Isentropic Vortex')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.5)
    
    summary_plot_path = os.path.join(output_dir, "convergence_p_refinement.png")
    plt.savefig(summary_plot_path)
    plt.close()
    print(f"\nStudy complete. Summary plot saved to {summary_plot_path}")

if __name__ == "__main__":
    run_study()
