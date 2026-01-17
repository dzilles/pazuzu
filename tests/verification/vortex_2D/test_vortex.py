import sys
import os
import numpy as np
import warp as wp
import matplotlib.pyplot as plt

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))

from solver import PazuzuSolver
from src.kernels.initial_conditions import init_isentropic_vortex
from src.kernels import boundary_conditions as bc
from typing import Any

import pytest

@wp.kernel
def init_isentropic_vortex_periodic(
    x: Any,
    y: Any,
    q: Any,
    active_indices: Any,
    num_active: int,
    params: Any,
    t: Any,
    beta: Any,
    radius: Any,
    center_x: Any,
    center_y: Any,
    x_min: Any,
    x_max: Any,
    y_min: Any,
    y_max: Any
):
    # Launch dimensions: (num_active, Np)
    block_idx, node_idx = wp.tid()
    
    pool_idx = active_indices[block_idx]
    
    xx = x[pool_idx, node_idx]
    yy = y[pool_idx, node_idx]
    
    # Vortex Parameters
    # Advecting with u_inf, v_inf
    x0 = center_x + params.u_inf * t 
    y0 = center_y + params.v_inf * t
    gamma = params.gamma
    
    # Periodic Distance Logic
    Lx = x_max - x_min
    Ly = y_max - y_min
    
    dx = xx - x0
    dy = yy - y0
    
    # Wrap to nearest image
    # dx = dx - Lx * round(dx / Lx)
    # Using generic Warp round
    
    # Note: For strict periodicity, we check if periodic domain is set (Lx > 0).
    # Assuming Lx > 0 for this test kernel.
    
    dx = dx - Lx * wp.round(dx / Lx)
    dy = dy - Ly * wp.round(dy / Ly)
    
    r2 = dx*dx + dy*dy
    r2_scaled = r2 / (radius * radius)
    
    template = beta
    one = bc.get_one_generic(template)
    half = bc.get_half_generic(template)
    two = one + one
    pi = bc.get_any_generic(template, 3.141592653589793)

    S_2pi = beta / (two * pi)
    exp_term = wp.exp(half * (one - r2_scaled))
    
    # Scale perturbation by radius to keep beta as peak velocity
    du = -S_2pi * (dy / radius) * exp_term
    dv =  S_2pi * (dx / radius) * exp_term
    
    u = params.u_inf + du
    v = params.v_inf + dv
    
    # Correct Isentropic Relation
    T_sub = (gamma - one) / gamma * half * (S_2pi * S_2pi) * wp.exp(one - r2_scaled)
    T = one - T_sub
    
    rho = wp.pow(T, one / (gamma - one))
    p = wp.pow(rho, gamma)
    
    E = p / (gamma - one) + half * rho * (u*u + v*v)
    
    q[pool_idx, node_idx] = bc.make_vec4_generic(rho, rho*u, rho*v, E)

def compute_errors(solver):
    """
    Computes L2 and Linf error for the density field without hardcoded parameters.
    """
    # 1. Get Numerical Solution
    q_num = solver.state.q.numpy()
    
    # Calculate domain size
    Lx = solver.config.mesh.x_max - solver.config.mesh.x_min
    Ly = solver.config.mesh.y_max - solver.config.mesh.y_min
    
    # Original Center
    x0 = float(solver.config.initial_condition.params.get('center_x', 0.0))
    y0 = float(solver.config.initial_condition.params.get('center_y', 0.0))
    
    # Note: We do NOT need to manually wrap center for the kernel,
    # because the kernel does the wrapping of (xx - x0).
    # We just need to pass the time and the original center.
    # But wait, init_isentropic_vortex_periodic computes x0 = center + u*t.
    # So we pass t=solver.state.t and center=original.
    
    # 2. Compute Exact Solution at t_final
    q_exact_wp = wp.zeros_like(solver.state.q)
    
    wp.launch(
        kernel=init_isentropic_vortex_periodic,
        dim=(solver.quadtree.num_blocks, solver.basis.Np),
        inputs=[
            solver.state.x,
            solver.state.y,
            q_exact_wp,
            solver.state.active_block_indices,
            solver.quadtree.num_blocks,
            solver.params,
            solver.scalar_dtype(solver.state.t), # Use actual simulation time
            solver.scalar_dtype(float(solver.config.initial_condition.params.get('beta', 5.0))),
            solver.scalar_dtype(float(solver.config.initial_condition.params.get('radius', 1.0))),
            solver.scalar_dtype(x0),
            solver.scalar_dtype(y0),
            solver.scalar_dtype(solver.config.mesh.x_min),
            solver.scalar_dtype(solver.config.mesh.x_max),
            solver.scalar_dtype(solver.config.mesh.y_min),
            solver.scalar_dtype(solver.config.mesh.y_max)
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
    # Error is around 2.2% with threshold 3%
    assert l2 < 3.0e-2
    print("PASS: AMR Vortex test successful.")

@pytest.mark.slow
def test_vortex_l2_error():
    config_path = os.path.join(os.path.dirname(__file__), "vortex.yaml")
    solver = PazuzuSolver(config_path)
    
    print(f"Running simulation: {solver.config.case_name}")
    print(f"Order: {solver.basis.N+1} (N={solver.basis.N}), Depth: {solver.config.mesh.initial_depth}")
    
    solver.run()
    
    l2, linf, diff_rho = compute_errors(solver)
    beta = float(solver.config.initial_condition.params.get('beta', 5.0))
    Lx = solver.config.mesh.x_max - solver.config.mesh.x_min
    expected = get_expected_error(solver.basis.N, solver.config.mesh.initial_depth, Lx, beta)
    
    print(f"\nFinal L2 Error:   {l2:.6e}")
    print(f"Final Linf Error: {linf:.6e}")
    print(f"Expected Error:   {expected:.6e}")
    
    # Check threshold (leeway on expected error)
    threshold = expected * 3.0
    
    assert l2 <= threshold, f"L2 Error {l2:.2e} exceeds threshold {threshold:.2e}"
    print(f"PASS: L2 Error is within expected bounds (Threshold: {threshold:.2e})")