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
    Reused from analyze_vortex.py
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
            scalar_dtype(radius)
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
    
    # Construct Configuration Programmatically
    config = PazuzuConfig(
        case_name=f"Vortex_{precision_mode}",
        solver_type=SolverType.EULER_2D,
        mesh=MeshConfig(
            x_min=-5.0, x_max=5.0,
            y_min=-5.0, y_max=5.0,
            periodic_x=True, periodic_y=True
        ),
        amr=AmrConfig(
            max_blocks=5000, # Sufficient for depth 4/5
            initial_depth=5
        ),
        initial_condition=InitialConditionConfig(
            name="vortex",
            params={"beta": 1.0, "radius": 1.0}
        ),
        physics=PhysicsConfig(
            gamma=1.4,
            u_inf=1.0, v_inf=1.0, p_inf=1.0,
            rho_inf=1.0
        ),
        numerics=NumericsConfig(
            polynomial_order=2, # N=2 (3rd order) where we saw difference
            cfl=0.1,
            precision=precision_mode
        ),
        io=IOConfig(
            output_dir=f"output/test_precision_{precision_mode}",
            write_interval=10000 # Don't write to disk to save time, unless needed
        ),
        simulation=SimulationConfig(
            t_final=10.0, # Match the manual run
            device="cuda"
        )
    )
    
    solver = PazuzuSolver(config)
    solver.run()
    
    l2, linf = compute_errors(solver)
    print(f"Precision: {precision_mode}, L2 Error: {l2:.6e}")
    
    # Cleanup Warp to free memory for next run?
    # Warp doesn't support full shutdown/reinit easily in one process.
    # But PazuzuSolver creates new arrays, so old ones should be GC'd if solver is deleted.
    
    return l2

def test_precision_improvement():
    """
    Verifies that Double Precision yields lower error than Single Precision
    for the Isentropic Vortex case.
    """
    # 1. Run Single Precision
    err_single = run_vortex_simulation("single")
    
    # 2. Run Double Precision
    err_double = run_vortex_simulation("double")
    
    print(f"\nComparison:")
    print(f"Single Precision Error: {err_single:.6e}")
    print(f"Double Precision Error: {err_double:.6e}")
    
    ratio = err_single / err_double
    print(f"Improvement Ratio (Single/Double): {ratio:.2f}x")
    
    # Assertions
    # 1. Double should be better
    assert err_double < err_single, "Double precision error should be lower than single precision."
    
    # 2. Significant improvement check (heuristic)
    # We saw ~7x improvement in manual tests. Let's be conservative and ask for >1.5x
    assert ratio > 1.5, f"Expected significant improvement (>1.5x), got {ratio:.2f}x"

if __name__ == "__main__":
    test_precision_improvement()
