from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.core.solver_builder import SolverBuilder
from src.numerics.time_steppers import RK4Stepper
from src.core.time_integrator import TimeIntegrator
from src.io.data_writer import HDF5Writer, HDF5Reader
from src.physics.initial_conditions import get_ic_function
from src.core.boundary_condition_manager import BoundaryConditionManager
import src.physics.laws.common as common_laws

import numpy as np
import warp as wp
import os
import argparse
import sys
import re

def get_unique_filename(path):
    """
    Returns a unique filename by appending _1, _2, etc. if the file exists.
    If the filename already ends in _N, it increments N.
    """
    if not os.path.exists(path):
        return path
        
    directory = os.path.dirname(path)
    filename = os.path.basename(path)
    base, ext = os.path.splitext(filename)
    
    # Pattern to match trailing _<number>
    match = re.search(r'_(\d+)$', base)
    if match:
        prefix = base[:match.start()]
        counter = int(match.group(1))
    else:
        prefix = base
        counter = 0
        
    while True:
        counter += 1
        new_filename = f"{prefix}_{counter}{ext}"
        new_path = os.path.join(directory, new_filename)
        if not os.path.exists(new_path):
            return new_path

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
    
    # Use SolverBuilder to create the solver instance
    builder = SolverBuilder(mesh, basis, config=cfg)
    solver = builder.build()
    
    # Initialize Solver State (Allocates memory)
    # We always run default initialization first to setup the structure
    print(f"Initializing Solver with IC: {ic_name}...")
    ic_func = get_ic_function(ic_name)
    solver.initialize(ic_func)

    # --- Restart / Continue Logic ---
    restart_file = cfg.simulation.restart_from
    if restart_file:
        if not os.path.isabs(restart_file):
            restart_file = os.path.join(config_dir, restart_file)
            
        print(f"Restarting from checkpoint: {restart_file}")
        
        # Load Checkpoint
        chk = HDF5Reader.load_checkpoint(restart_file)
        time = chk["time"]
        step = chk["step"]
        q_prim_flat = chk["q_prim_flat"] # (4, TotalPoints)
        
        # Validate Shape
        expected_points = mesh.num_elements * basis.Np
        if q_prim_flat.shape[1] != expected_points:
            raise ValueError(f"Restart file geometry mismatch. Expected {expected_points} points, got {q_prim_flat.shape[1]}. Ensure polynomial order matches.")
            
        # Convert to Conservative
        # Update global gamma for common laws to ensure correct conversion
        common_laws.gamma = cfg.physics.gamma
        
        q_cons_flat = common_laws.primitive_to_conservative(q_prim_flat)
        
        # Reshape: (4, NumElements * Np) -> (4, NumElements, Np) -> (NumElements, Np, 4)
        q_cons = q_cons_flat.reshape(4, mesh.num_elements, basis.Np).transpose(1, 2, 0)
        
        # Override State
        solver.state.q = wp.array(q_cons, dtype=solver.state.dtype, device=solver.device)
        solver.state.t = time
        solver.state.step = step
        
        print(f"Resumed state at t={time:.4f}, step={step}")

    # --- Output Directory ---
    output_dir = cfg.io.output_dir
    if not os.path.isabs(output_dir):
        output_dir = os.path.join(config_dir, output_dir)
        
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_file = os.path.join(output_dir, "results.h5")
    
    # Ensure unique output filename to avoid overwriting or conflict with restart file
    output_file = get_unique_filename(output_file)
    
    writer = HDF5Writer(
        output_file, 
        mesh,
        basis=solver.basis,
        node_coords=(solver.mesh.x_host, solver.mesh.y_host)
    )

    # Save Initial/Restart State
    print(f"Saving start state to {os.path.basename(output_file)} (Step {solver.state.step})...")
    writer.write_step(solver.state.step, solver.state.t, solver.state.numpy())

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
