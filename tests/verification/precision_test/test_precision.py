import sys
import os
import numpy as np
import warp as wp
import pytest

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))

from solver import PazuzuSolver
from src.core.config import PazuzuConfig, SolverType, MeshConfig, AmrConfig, NumericsConfig, PhysicsConfig, IOConfig, SimulationConfig, InitialConditionConfig
from src.kernels.initial_conditions import init_isentropic_vortex

def compute_errors(solver):
    """
    Computes L2 and Linf error for the density field without hardcoded parameters.
    Reused from test_vortex.py
    """
    # 1. Get Numerical Solution
    q_num = solver.state.q.numpy() 
    
    # Calculate domain size
    Lx = solver.config.mesh.x_max - solver.config.mesh.x_min
    Ly = solver.config.mesh.y_max - solver.config.mesh.y_min

    # Calculate effective time for periodic boundaries
    t_eff = solver.state.t
    if solver.config.mesh.periodic_x and abs(solver.params.u_inf) > 1e-9:
        period = Lx / abs(solver.params.u_inf)
        t_eff = solver.state.t % period
        if abs(t_eff - period) < 1e-9:
            t_eff = 0.0

    # 2. Compute Exact Solution at t_final
    # Initialize a temporary state for exact solution with correct dtype
    scalar_dtype = solver.scalar_dtype
    vec4_type = wp.vec4d if scalar_dtype == wp.float64 else wp.vec4
    
    q_exact_wp = wp.zeros((solver.state.max_blocks, solver.basis.Np), dtype=vec4_type, device=solver.device)
    
    # Cast params to scalar_dtype
    beta = solver.config.initial_condition.params.get('beta', 5.0)
    radius = solver.config.initial_condition.params.get('radius', 1.0)
    center_x = solver.config.initial_condition.params.get('center_x', 0.0)
    center_y = solver.config.initial_condition.params.get('center_y', 0.0)
    
    wp.launch(
        kernel=init_isentropic_vortex,
        dim=(solver.quadtree.num_blocks, solver.basis.Np),
        inputs=[
            solver.state.x,
            solver.state.y,
            q_exact_wp,
            solver.state.active_block_indices,
            solver.quadtree.num_blocks,
            solver.params,
            scalar_dtype(t_eff),
            scalar_dtype(beta),
            scalar_dtype(radius),
            scalar_dtype(center_x),
            scalar_dtype(center_y)
        ],
        device=solver.device
    )
    
    q_exact = q_exact_wp.numpy()
    
    # 3. Compute Error Norms
    w = solver.basis.weights_1d.numpy()
    N1 = solver.basis.N1
    w_2d = np.zeros(solver.basis.Np)
    for j in range(N1):
        for i in range(N1):
            w_2d[j*N1 + i] = w[i] * w[j]
            
    # Jacobian Calculation (Accounts for Domain Size and Depth)
    level = solver.config.amr.initial_depth
    grid_dim = 1 << level
    
    hx = Lx / grid_dim
    hy = Ly / grid_dim
    detJ = (hx / 2.0) * (hy / 2.0)
    
    # Slice to active blocks
    num_active = solver.quadtree.num_blocks
    diff_rho = q_num[:num_active, :, 0] - q_exact[:num_active, :, 0]
    
    error_sq = 0.0
    max_err = 0.0
    
    # Sum over active blocks and nodes
    # Note: numpy loops can be slow, but for verification this is acceptable
    for b in range(num_active):
        for n in range(solver.basis.Np):
            val = diff_rho[b, n]
            error_sq += val**2 * w_2d[n] * detJ
            max_err = max(max_err, np.abs(val))
            
    L2_error = np.sqrt(error_sq)
    return L2_error, max_err

def run_vortex_simulation(precision_mode):
    print(f"\n--- Running Simulation with precision={precision_mode} ---")
    
    # Ensure output is in the local output folder
    output_base = os.path.join(os.path.dirname(__file__), "output")
    
    # "Nano-Vortex" Configuration
    # We center the domain at 1.0.
    # In float32, machine epsilon at 1.0 is ~1.19e-7.
    # dx = 1e-3 / (16 * 2) = 3.125e-5.
    # At 1e6, float32 epsilon is ~0.0625. dx << epsilon.
    
    domain_center = 1.0e6
    domain_width = 1.0e-3
    
    initial_depth = 4

    wp.init()
    device = "cuda" if wp.is_cuda_available() else "cpu"
    
    config = PazuzuConfig(
        case_name=f"Vortex_{precision_mode}",
        solver_type=SolverType.EULER_2D,
        mesh=MeshConfig(
            x_min=domain_center - domain_width/2, 
            x_max=domain_center + domain_width/2,
            y_min=domain_center - domain_width/2, 
            y_max=domain_center + domain_width/2,
            periodic_x=True, periodic_y=True
        ),
        amr=AmrConfig(
            max_blocks=500, 
            initial_depth=initial_depth
        ),
        initial_condition=InitialConditionConfig(
            name="vortex",
            params={
                "beta": 5.0, 
                "radius": 0.2 * domain_width, 
                "center_x": domain_center,
                "center_y": domain_center
            }
        ),
        physics=PhysicsConfig(
            gamma=1.4,
            u_inf=1.0, v_inf=1.0, p_inf=1.0,
            rho_inf=1.0
        ),
        numerics=NumericsConfig(
            polynomial_order=2, 
            cfl=0.1,
            precision=precision_mode,
            dt_static=None # Use dynamic DT
        ),
        io=IOConfig(
            output_dir=os.path.join(output_base, f"test_precision_{precision_mode}"),
            write_interval=100
        ),
        simulation=SimulationConfig(
            t_final=domain_width / 1.0, # Time to cross domain
            device=device
        )
    )
    
    try:
        solver = PazuzuSolver(config)
        solver.run()
        l2, linf = compute_errors(solver)
        print(f"Precision: {precision_mode}, L2 Error: {l2:.6e}")
        return l2
    except Exception as e:
        print(f"Simulation FAILED for {precision_mode} as expected or unexpected: {e}")
        return None

def test_nano_vortex_precision():
    """
    Verifies that Single Precision fails for a 'Nano-Vortex' at coordinates ~1.0,
    due to machine epsilon limits, while Double Precision succeeds.
    """
    print("\nStarting Nano-Vortex Precision Stress Test...")
    
    # 1. Run Single Precision
    # We expect this to likely crash (Jacobian=0) or produce garbage.
    # If it produces 0.0 error, it likely means the grid collapsed to a constant field.
    err_single = run_vortex_simulation("single")
    
    # 2. Run Double Precision
    # We expect this to run correctly and produce a valid (non-zero) discretization error.
    err_double = run_vortex_simulation("double")
    
    print(f"\nResults:")
    print(f"Single Precision Result: {err_single}")
    print(f"Double Precision Result: {err_double}")
    
    # Assertions
    
    # Check if Single failed appropriately
    single_failed = False
    if err_single is None:
        print("SUCCESS: Single precision crashed (likely grid collapse).")
        single_failed = True
    elif err_single < 1e-20:
         print("SUCCESS: Single precision yielded 0 error (likely grid collapsed to constant field).")
         single_failed = True
    else:
         print(f"WARNING: Single precision produced error {err_single}. Checking if Double is better...")
    
    # Check Double
    assert err_double is not None, "Double precision crashed!"
    assert err_double > 1e-20, f"Double precision error {err_double} is too small (vortex missed?)"
    
    if not single_failed:
        assert err_double < err_single, "Double should be better than single"
        ratio = err_single / err_double
        print(f"Error Ratio: {ratio:.2f}")
        assert ratio > 10.0, "Single precision should be significantly worse."

if __name__ == "__main__":
    test_nano_vortex_precision()
