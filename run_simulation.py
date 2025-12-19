from src.dg_solver.mesh import Mesh
from src.dg_solver.basis import Basis
from src.dg_solver.solver import DGSolver

import numpy as np
import os
import argparse

def initial_conditions_vortex(x, y):
    """
    2D Isentropic Vortex problem.
    This is a common test case for Euler solvers.
    The solution is a vortex that advects diagonally across the domain.
    """
    gamma = 1.4
    beta = 5.0  # Vortex strength
    r_sq = (x - 5.0)**2 + (y - 5.0)**2
    
    u_inf = 1.0
    v_inf = 1.0
    
    du = -(beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * (y - 5.0)
    dv = (beta / (2 * np.pi)) * np.exp(0.5 * (1 - r_sq)) * (x - 5.0)
    
    u = u_inf + du
    v = v_inf + dv
    
    T_inf = 1.0
    T = T_inf - ((gamma - 1) * beta**2 / (8 * gamma * np.pi**2)) * np.exp(1 - r_sq)
    
    rho = T**(1.0 / (gamma - 1))
    p = rho**gamma
    
    return rho, u, v, p

def initial_conditions_uniform(x, y):
    """
    Uniform flow.
    """
    rho = 1.0
    u = 1.0
    v = 0.0
    p = 1.0
    return rho, u, v, p


def main(args):
    # Simulation parameters
    polynomial_degree = 4
    nx = args.nx
    ny = args.ny
    t_final =20.0
    device = args.device

    # --- Domain ---
    x_min, x_max = 0.0, 10.0
    y_min, y_max = 0.0, 10.0
    
    # --- CFL condition ---
    # dt = CFL * dx / ( |u| + c )
    # For DG, CFL is restricted by ~1/(2N+1), where N is the polynomial degree.
    # For N=3, CFL < 0.14. We choose a safe value.
    CFL = 0.05

    # Create mesh and basis
    print(f"Creating mesh ({nx}x{ny}) and basis (N={polynomial_degree})...")
    mesh = Mesh(nx=nx, ny=ny, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, device=device)
    basis = Basis(polynomial_degree=polynomial_degree, device=device)

    # Create the solver
    print("Initializing solver...")
    solver = DGSolver(mesh, basis, initial_conditions_func=initial_conditions_vortex)


    output_dir = "data"
    np.savez(
        os.path.join(output_dir, "init.npz"), 
        q=solver.Q.numpy(), 
        x=solver.x.numpy(), 
        y=solver.y.numpy(),
        vertices=mesh.vertices_host
    )

    # Run the simulation
    print(f"Starting simulation on device '{device}'...")
    print(f"t_final = {t_final}, CFL = {CFL}")
    solver.solve(t_final=t_final, CFL=CFL, log_frequency=20)
    print("Simulation finished.")

    # Save the solution
    print("Saving solution to data/output.npz...")
    output_dir = "data"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    np.savez(
        os.path.join(output_dir, "output.npz"), 
        q=solver.Q.numpy(), 
        x=solver.x.numpy(), 
        y=solver.y.numpy(),
        vertices=mesh.vertices_host
    )

    print("Solution saved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="2D Discontinuous Galerkin Euler Solver with Nvidia Warp")
    parser.add_argument("--poly_degree", type=int, default=3, help="Polynomial degree for the basis functions.")
    parser.add_argument("--nx", type=int, default=32, help="Number of elements in the x-direction.")
    parser.add_argument("--ny", type=int, default=32, help="Number of elements in the y-direction.")
    parser.add_argument("--t_final", type=float, default=20.0, help="Final simulation time.")
    parser.add_argument("--device", type=str, default="cuda", choices=["cuda", "cpu"], help="Device to run the simulation on.")
    
    args = parser.parse_args()
    main(args)
