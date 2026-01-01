import sys
import os
import numpy as np
import warp as wp
import matplotlib.pyplot as plt

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))

from solver import PazuzuSolver
from src.kernels.initial_conditions import init_isentropic_vortex

import pytest

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
            solver.scalar_dtype(t_eff),
            solver.scalar_dtype(float(solver.config.initial_condition.params.get('beta', 5.0))),
            solver.scalar_dtype(float(solver.config.initial_condition.params.get('radius', 1.0))),
            solver.scalar_dtype(float(solver.config.initial_condition.params.get('center_x', 0.0))),
            solver.scalar_dtype(float(solver.config.initial_condition.params.get('center_y', 0.0)))
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
            
    # Slice to active blocks
    num_active = solver.quadtree.num_blocks
    active = solver.state.active_block_indices.numpy()[:num_active]
    
    # We need to compute diff only for active blocks
    # diff_rho shape will be (num_active, Np)
    diff_rho = np.zeros((num_active, solver.basis.Np))
    for b in range(num_active):
        pool_idx = active[b]
        diff_rho[b, :] = q_num[pool_idx, :, 0] - q_exact[pool_idx, :, 0]
    
    error_sq = 0.0
    max_err = 0.0
    
    levels = solver.quadtree.block_levels.numpy()
    
    for b in range(num_active):
        # Jacobian Calculation for this specific block level
        pool_idx = active[b]
        level = levels[pool_idx]
        grid_dim = 1 << level
        hx = Lx / grid_dim
        hy = Ly / grid_dim
        detJ = (hx / 2.0) * (hy / 2.0)
        
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

@pytest.mark.slow
def test_vortex_uniform_L3_error():
    config_path = os.path.join(os.path.dirname(__file__), "vortex_uniform_L3.yaml")
    solver = PazuzuSolver(config_path)
    print(f"Running simulation: {solver.config.case_name}")
    solver.run()
    l2, linf, diff_rho = compute_errors(solver)
    print(f"\nFinal L2 Error:   {l2:.6e}")
    print(f"Final Linf Error: {linf:.6e}")

@pytest.mark.slow
def test_vortex_uniform_L2_error():
    config_path = os.path.join(os.path.dirname(__file__), "vortex_uniform_L2.yaml")
    solver = PazuzuSolver(config_path)
    print(f"Running simulation: {solver.config.case_name}")
    solver.run()
    l2, linf, diff_rho = compute_errors(solver)
    print(f"\nFinal L2 Error:   {l2:.6e}")
    print(f"Final Linf Error: {linf:.6e}")

@pytest.mark.slow
def test_vortex_amr_error():
    config_path = os.path.join(os.path.dirname(__file__), "vortex_amr.yaml")
    solver = PazuzuSolver(config_path)
    
    print(f"Running simulation: {solver.config.case_name}")
    print(f"Initial Blocks: {solver.quadtree.num_blocks}")
    
    # Debug: Check coordinates
    x_np = solver.state.x.numpy()[:solver.quadtree.num_blocks]
    print(f"X range: {np.min(x_np):.2f} to {np.max(x_np):.2f}")
    
    # Check Error at t=0
    l2_0, linf_0, _ = compute_errors(solver)
    print(f"Initial L2 Error: {l2_0:.6e}")
    
    solver.run()
    
    l2, linf, diff_rho = compute_errors(solver)
    print(f"\nFinal L2 Error:   {l2:.6e}")
    print(f"Final Linf Error: {linf:.6e}")
    
    # Assert reasonable error for AMR
    # Since it's a mix of levels, we just check if it's small (stable)
    assert l2 < 1.0e-2
    print("PASS: AMR Vortex test successful.")

@pytest.mark.slow
def test_vortex_l2_error():
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
    
    assert l2 <= threshold, f"L2 Error {l2:.2e} exceeds threshold {threshold:.2e}"
    print(f"PASS: L2 Error is within expected bounds (Threshold: {threshold:.2e})")