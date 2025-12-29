import sys
import os
import numpy as np
import warp as wp
import matplotlib.pyplot as plt

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../")))

from solver import PazuzuSolver
from src.kernels.initial_conditions import init_isentropic_vortex

def compute_l2_error(solver):
    # 1. Get Numerical Solution
    q_num = solver.state.q.numpy()
    
    # 2. Compute Exact Solution at t_final
    # We can reuse the IC kernel to generate the exact solution at time t
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
            solver.state.t, # t_final
            float(solver.config.initial_condition.params.get('beta', 5.0)),
            float(solver.config.initial_condition.params.get('radius', 1.0))
        ],
        device=solver.device
    )
    
    q_exact = q_exact_wp.numpy()
    
    # 3. Compute Error Norm
    # L2^2 = sum (q_num - q_exact)^2 * w_i * w_j * J
    
    # Weights (1D)
    w = solver.basis.weights_1d.numpy()
    # 2D weights (tensor product)
    # w_2d[node_idx] = w[i] * w[j]
    N1 = solver.basis.N1
    w_2d = np.zeros(solver.basis.Np)
    for j in range(N1):
        for i in range(N1):
            w_2d[j*N1 + i] = w[i] * w[j]
            
    # Jacobian (Uniform Grid)
    # Domain [-1, 1]x[-1, 1] -> 2x2 size.
    # Level L -> Grid Dim 2^L.
    # Block size h = 2 / 2^L.
    # Reference element [-1, 1]x[-1, 1] (size 2x2).
    # Mapping [-1, 1] -> [x0, x0+h]. Scaling factor h/2.
    # J = (h/2) * (h/2) = h^2 / 4.
    
    level = 3 # From solver defaults (should extract from quadtree if variable)
    grid_dim = 1 << level
    h = 2.0 / grid_dim
    detJ = (h / 2.0) * (h / 2.0)
    
    # Compute Density Error (component 0)
    diff_rho = q_num[:, :, 0] - q_exact[:, :, 0]
    
    # Weighted Sum
    # Sum over blocks, Sum over nodes
    error_sq = 0.0
    
    for b in range(solver.quadtree.num_blocks):
        for n in range(solver.basis.Np):
            error_sq += diff_rho[b, n]**2 * w_2d[n] * detJ
            
    L2_error = np.sqrt(error_sq)
    return L2_error

def main():
    config_path = os.path.join(os.path.dirname(__file__), "vortex.yaml")
    solver = PazuzuSolver(config_path)
    
    print(f"Running simulation for {solver.config.simulation.t_final}s...")
    solver.run()
    
    error = compute_l2_error(solver)
    print(f"\nFinal L2 Error (Density): {error:.6e}")
    
    # N=3 (4th order). Error should be small.
    # Expected convergence is O(h^(N+1)).
    
if __name__ == "__main__":
    main()
