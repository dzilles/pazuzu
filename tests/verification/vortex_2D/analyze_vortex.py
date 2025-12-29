import sys
import os
import numpy as np
import warp as wp
import matplotlib.pyplot as plt

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))

from solver import PazuzuSolver
from src.kernels.initial_conditions import init_isentropic_vortex

def compute_errors(solver):
    """
    Computes L2 and Linf error for the density field without hardcoded parameters.
    """
    # 1. Get Numerical Solution
    q_num = solver.state.q.numpy()
    
    # Calculate domain size
    Lx = solver.config.mesh.x_max - solver.config.mesh.x_min
    Ly = solver.config.mesh.y_max - solver.config.mesh.y_min

    # Calculate effective time for periodic boundaries
    # The exact solution kernel moves the vortex indefinitely.
    # For periodic domains, we wrap the time so the vortex center matches the wrapped domain.
    t_eff = solver.state.t
    if solver.config.mesh.periodic_x and abs(solver.params.u_inf) > 1e-9:
        period = Lx / abs(solver.params.u_inf)
        t_eff = solver.state.t % period
        # Handle case where modulo result is very close to period (float precision)
        if abs(t_eff - period) < 1e-9:
            t_eff = 0.0

    # 2. Compute Exact Solution at t_final
    q_exact_wp = wp.zeros_like(solver.state.q)
    
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
            t_eff,
            float(solver.config.initial_condition.params.get('beta', 5.0)),
            float(solver.config.initial_condition.params.get('radius', 1.0))
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
    
    for b in range(num_active):
        for n in range(solver.basis.Np):
            val = diff_rho[b, n]
            error_sq += val**2 * w_2d[n] * detJ
            max_err = max(max_err, np.abs(val))
            
    L2_error = np.sqrt(error_sq)
    return L2_error, max_err, diff_rho

def get_expected_error(N, level, Lx, beta=5.0):
    """
    Returns an order-of-magnitude expected L2 error for the vortex case.
    h ~ Lx/2^level. Error ~ C * h^(N+1).
    """
    h = Lx / (2**level) # Reference length
    # Heuristic C constant for isentropic vortex
    C = 0.05 * beta 
    return C * (h**(N + 1))

def main():
    config_path = os.path.join(os.path.dirname(__file__), "vortex.yaml")
    solver = PazuzuSolver(config_path)
    
    print(f"Running simulation: {solver.config.case_name}")
    print(f"Order: {solver.basis.N+1} (N={solver.basis.N}), Depth: {solver.config.amr.initial_depth}")
    
    solver.run()
    
    l2, linf, diff_rho = compute_errors(solver)
    beta = float(solver.config.initial_condition.params.get('beta', 5.0))
    Lx = solver.config.mesh.x_max - solver.config.mesh.x_min
    expected = get_expected_error(solver.basis.N, solver.config.amr.initial_depth, Lx, beta)
    
    print(f"\nFinal L2 Error:   {l2:.6e}")
    print(f"Final Linf Error: {linf:.6e}")
    print(f"Expected Error:   {expected:.6e}")
    
    # Check threshold (leeway on expected error)
    threshold = expected * 3.0
    
    if l2 > threshold:
         print(f"FAIL: L2 Error {l2:.2e} exceeds threshold {threshold:.2e}")
         sys.exit(1)
         
    print(f"PASS: L2 Error is within expected bounds (Threshold: {threshold:.2e})")
    
if __name__ == "__main__":
    main()