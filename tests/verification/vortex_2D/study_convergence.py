import sys
import os
import numpy as np
import h5py
import subprocess
import matplotlib.pyplot as plt
import yaml
from scipy import stats
from scipy.special import legendre

# Add root directory to path to find solver if needed
test_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(test_dir, "../../.."))
sys.path.append(root_dir)

def get_gll_nodes_weights(N):
    if N == 0: return np.array([0.0]), np.array([2.0])
    if N == 1: roots = np.array([])
    else: roots = np.roots(legendre(N).deriv(1))
    nodes = np.concatenate(([-1.0], np.sort(roots), [1.0]))
    weights = 2 / (N * (N + 1) * legendre(N)(nodes)**2)
    return nodes, weights

def periodic_vortex(x, y, t=0.0, Lx=10.0, Ly=10.0):
    gamma = 1.4
    beta = 5.0  # Vortex strength
    u_inf = 1.0
    v_inf = 1.0
    xc, yc = 0.0, 0.0
    xc_t = xc + u_inf * t
    yc_t = yc + v_inf * t
    dx = x - xc_t
    dy = y - yc_t
    # Periodic wrapping
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
    # depths = [5, 6, 7] -> Res = 32, 64, 128 blocks per edge
    depths = [5, 6, 7] 
    orders = [1, 2, 3]
    
    all_results = {}
    
    base_config_path = os.path.join(test_dir, "vortex.yaml")
    output_dir = os.path.join(test_dir, "output")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with open(base_config_path, 'r') as f:
        base_config = yaml.safe_load(f)

    # Ensure precision is double for convergence study
    base_config['numerics']['precision'] = "double"
    base_config['simulation']['device'] = "cuda" 
    base_config['amr']['max_blocks'] = 70000 
    
    for P in orders:
        print(f"\n========================================")
        print(f"  Studying Polynomial Order P = {P}")
        print(f"========================================")
        
        errors = []
        h_values = []
        
        _, w1d = get_gll_nodes_weights(P)
        weights_2d = np.kron(w1d, w1d)
        
        for depth in depths:
            N_blocks_edge = 2**depth
            print(f"\n--- Running P={P}, Depth={depth} ({N_blocks_edge}x{N_blocks_edge} blocks) ---")
            
            run_output_dir = os.path.join(output_dir, f"run_P{P}_D{depth}")
            if not os.path.exists(run_output_dir):
                os.makedirs(run_output_dir)
                
            config = base_config.copy()
            config['numerics'] = base_config['numerics'].copy()
            config['io'] = base_config['io'].copy()
            config['amr'] = base_config['amr'].copy()
            
            config['numerics']['polynomial_order'] = P
            config['numerics']['cfl'] = 0.1
            config['simulation']['t_final'] = 1.0
            config['numerics']['flux'] = "hllc"
            config['numerics']['time_integrator'] = "rk4"
            config['amr']['initial_depth'] = depth
            config['io']['output_dir'] = run_output_dir
            
            temp_config_path = os.path.join(test_dir, f"vortex_P{P}_D{depth}.yaml")
            with open(temp_config_path, 'w') as f:
                yaml.dump(config, f)
                
            # 3. Run Simulation
            # Use the .venv python if it exists
            if sys.platform == "win32":
                python_exe = os.path.join(root_dir, ".venv", "Scripts", "python.exe")
            else:
                python_exe = os.path.join(root_dir, ".venv", "bin", "python")
            
            if not os.path.exists(python_exe):
                python_exe = sys.executable

            cmd = [python_exe, "solver.py", temp_config_path]
            try:
                subprocess.run(cmd, cwd=root_dir, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Simulation failed for P={P}, Depth={depth}: {e}")
                continue
            
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
                
                Lx = base_config['mesh']['x_max'] - base_config['mesh']['x_min']
                Ly = base_config['mesh']['y_max'] - base_config['mesh']['y_min']
                
                rho_ana, _, _, _ = periodic_vortex(x, y, t=time, Lx=Lx, Ly=Ly)
                
                diff = rho_num - rho_ana
                
                # Discrete L2 norm (RMS) as in main branch
                l2_error = np.sqrt(np.mean(diff**2))

                dx = Lx / N_blocks_edge
                
                errors.append(l2_error)
                h_values.append(dx)

                if len(errors) == 1:
                    print(f"L2 Error: {l2_error:.6e} (Baseline)")
                else:
                    actual_ratio = errors[-2] / errors[-1]
                    h_ratio = h_values[-2] / h_values[-1]
                    actual_slope = np.log(actual_ratio) / np.log(h_ratio)
                    print(f"L2 Error: {l2_error:.6e} | Observed Order: {actual_slope:.2f}")
        
        if not errors:
            continue

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
        
        print(f"P = {P}: Global Order of Convergence (Slope) = {slope:.3f}")

        plt.figure(figsize=(8, 6))
        plt.loglog(h, err, 'o-', label=f'Numerical P={P} (Slope={slope:.2f})')
        
        expected_slope = P + 1
        h_ref = np.linspace(min(h), max(h), 100)
        err_ref = err[0] * (h_ref / h[0])**expected_slope
        plt.loglog(h_ref, err_ref, 'k--', alpha=0.5, label=f'Reference O(h^{expected_slope})')
        
        plt.xlabel('Block size h')
        plt.ylabel('L2 Error (Density)')
        plt.title(f'Grid Convergence Study: P={P}')
        plt.legend()
        plt.grid(True, which="both", ls="-", alpha=0.5)
        
        plot_path = os.path.join(output_dir, f"convergence_P{P}.png")
        plt.savefig(plot_path)
        plt.close()

    # Final Summary Plot
    if all_results:
        plt.figure(figsize=(10, 8))
        for P in all_results:
            h = all_results[P]['h']
            err = all_results[P]['errors']
            log_h = np.log(h)
            log_err = np.log(err)
            slope, _, _, _, _ = stats.linregress(log_h, log_err)
            plt.loglog(h, err, 'o-', label=f'P={P} (Slope={slope:.2f})')

        plt.xlabel('Block size h')
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
