from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.core.solver_builder import SolverBuilder
from src.numerics.time_steppers import RK4Stepper
from src.core.driver import TimeIntegrator
from src.io.data_writer import HDF5Writer
from src.physics.initial_conditions import get_ic_function
from src.core.boundary_condition_manager import BoundaryConditionManager

import numpy as np
import warp as wp
import os
import argparse
import sys

def main(config_path):
    from src.core.config import PazuzuConfig
    
    print(f"Loading configuration from {config_path}...")
    try:
        cfg = PazuzuConfig.from_yaml(config_path)
    except Exception as e:
        print(f"[Error] Configuration Error:\n{e}")
        sys.exit(1)

    print(f"[Success] Configuration '{cfg.case_name}' loaded successfully.")
    
    # Get the directory of the config file to resolve relative paths
    config_dir = os.path.dirname(os.path.abspath(config_path))
    
    # --- Mesh Setup ---
    mesh_path = cfg.mesh_file
    if not os.path.isabs(mesh_path):
        mesh_path = os.path.join(config_dir, mesh_path)
        
    device = cfg.simulation.device
    print(f"Initializing Mesh from {mesh_path} on {device}...")
    if not os.path.exists(mesh_path):
         raise FileNotFoundError(f"Mesh file not found: {mesh_path}")
    
    mesh = Mesh(filename=mesh_path, device=device)

    # --- Periodic BC Setup ---
    print("Applying Periodic Boundary Conditions (if any)...")
    BoundaryConditionManager.apply_periodic_conditions(mesh, cfg)

    # --- Basis Setup ---
    dtype_warp = wp.float64 if cfg.numerics.precision == "double" else wp.float32
    print(f"Initializing Basis (Poly Degree N={cfg.numerics.polynomial_order}, Precision={cfg.numerics.precision})...")
    basis = Basis(
        cfg.numerics.polynomial_order, 
        device=cfg.simulation.device, 
        dtype=dtype_warp,
        over_integration_order=cfg.numerics.over_integration_order
    )

    # --- Solver Setup ---
    ic_cfg = cfg.initial_condition
    ic_name = ic_cfg if isinstance(ic_cfg, str) else ic_cfg.name
    
    print(f"Initializing Solver with IC: {ic_name}...")
    ic_func = get_ic_function(ic_name)
    
    # Use SolverBuilder to create the solver instance
    builder = SolverBuilder(mesh, basis, config=cfg)
    solver = builder.build()
    solver.initialize(ic_func)

    # --- Output Directory ---
    output_dir = cfg.io.output_dir
    if not os.path.isabs(output_dir):
        output_dir = os.path.join(config_dir, output_dir)
        
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    writer = HDF5Writer(
        os.path.join(output_dir, "results.h5"), 
        mesh,
        basis=solver.basis,
        node_coords=(solver.mesh.x_host, solver.mesh.y_host)
    )

    # Save Initial State
    print("Saving initial state...")
    writer.write_step(0, 0.0, solver.state.numpy())

    # --- Run Simulation ---
    print(f"Starting simulation on device '{device}'...")
    
    # Instantiate Stepper and Driver
    stepper = RK4Stepper(device=device)
    driver = TimeIntegrator(solver, stepper)
    # driver.solve now retrieves t_final, CFL, write_interval from solver.config
    driver.solve(writer=writer, max_steps=cfg.simulation.max_steps)
    
    print("Simulation finished.")
    print(f"Results saved to {output_dir}/results.h5 and .xmf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="2D Discontinuous Galerkin Euler Solver with Nvidia Warp")
    parser.add_argument("config", type=str, nargs='?', default="config/default.yaml", help="Path to the YAML configuration file.")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.config):
        print(f"Error: Config file '{args.config}' not found.")
        sys.exit(1)
        
    main(args.config)
