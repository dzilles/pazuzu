from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.core.solver_builder import SolverBuilder
from src.numerics.time_steppers import RK4Stepper
from src.core.driver import TimeIntegrator
from src.io.data_writer import HDF5Writer
from src.physics import initial_conditions as ic

import numpy as np
import warp as wp
import os
import argparse
import yaml
import sys

def get_initial_condition_func(name):
    if name == "vortex":
        return ic.vortex
    elif name == "uniform":
        return ic.uniform
    elif name == "rest":
        return ic.rest
    elif name == "sod_shock_tube":
        return ic.sod_shock_tube
    else:
        raise ValueError(f"Unknown initial condition: {name}")

def main(config_path):
    from src.core.config import PazuzuConfig
    
    print(f"Loading configuration from {config_path}...")
    try:
        cfg = PazuzuConfig.from_yaml(config_path)
    except Exception as e:
        print(f"❌ Configuration Error:\n{e}")
        sys.exit(1)

    print(f"✅ Configuration '{cfg.case_name}' loaded successfully.")
    
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
    if cfg.periodic_pairs:
        print("Applying Periodic Boundary Conditions...")
        for pair in cfg.periodic_pairs:
            # Expecting [Tag1, Tag2, Axis]
            if len(pair) != 3:
                print(f"  Warning: Invalid periodic pair format: {pair}. Expected [Tag1, Tag2, Axis].")
                continue
                
            t1_raw, t2_raw, axis = pair
            
            # Helper to resolve tag name to ID
            def resolve_tag(raw):
                if isinstance(raw, int): return raw
                if hasattr(mesh, 'physical_groups') and raw in mesh.physical_groups:
                    return mesh.physical_groups[raw]
                try: return int(raw)
                except: return -1
            
            t1 = resolve_tag(t1_raw)
            t2 = resolve_tag(t2_raw)
            
            if t1 == -1 or t2 == -1:
                print(f"  Warning: Could not resolve tags for periodic pair: {pair}")
            else:
                mesh.apply_periodic_condition(t1, t2, axis)

    # --- Basis Setup ---
    dtype_warp = wp.float64 if cfg.numerics.precision == "double" else wp.float32
    print(f"Initializing Basis (Poly Degree N={cfg.numerics.polynomial_order}, Precision={cfg.numerics.precision})...")
    basis = Basis(polynomial_degree=cfg.numerics.polynomial_order, device=device, dtype=dtype_warp)

    # --- Solver Setup ---
    ic_name = cfg.initial_condition
    print(f"Initializing Solver with IC: {ic_name}...")
    ic_func = get_initial_condition_func(ic_name)
    
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
