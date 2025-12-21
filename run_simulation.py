from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.physics.euler_2d import Euler2DSolver
from src.core.driver import TimeIntegrator
from src.io.data_writer import HDF5Writer
from src.physics import initial_conditions as ic

import numpy as np
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
    else:
        raise ValueError(f"Unknown initial condition: {name}")

def main(config_path):
    
    print(f"Loading configuration from {config_path}...")
    # Get the directory of the config file to resolve relative paths
    config_dir = os.path.dirname(os.path.abspath(config_path))
    
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)

    # --- Extract Settings ---
    sim_cfg = cfg['simulation']
    mesh_cfg = cfg['mesh']
    num_cfg = cfg['numerical']
    ic_cfg = cfg['initial_condition']

    # Simulation parameters
    polynomial_degree = num_cfg['polynomial_degree']
    t_final = sim_cfg['t_final']
    device = sim_cfg.get('device', 'cpu') # Default to CPU for safety
    CFL = sim_cfg['cfl']
    log_freq = sim_cfg.get('log_frequency', 10)
    
    # Resolve output directory
    output_dir = sim_cfg['output_dir']
    if not os.path.isabs(output_dir):
        output_dir = os.path.join(config_dir, output_dir)

    # --- Mesh Setup ---
    print("Initializing Mesh...")
    if mesh_cfg['type'] == 'cartesian':
        nx = mesh_cfg['nx']
        ny = mesh_cfg['ny']
        x_min = mesh_cfg['x_min']
        x_max = mesh_cfg['x_max']
        y_min = mesh_cfg['y_min']
        y_max = mesh_cfg['y_max']
        print(f"  Type: Cartesian ({nx}x{ny})")
        mesh = Mesh(nx=nx, ny=ny, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, device=device)
    
    elif mesh_cfg['type'] == 'unstructured':
        filename = mesh_cfg['filename']
        # Resolve mesh filename
        if not os.path.isabs(filename):
            filename = os.path.join(config_dir, filename)
            
        print(f"  Type: Unstructured (loading from {filename})")
        if not os.path.exists(filename):
             raise FileNotFoundError(f"Mesh file not found: {filename}")
        mesh = Mesh(filename=filename, device=device)
    
    else:
        raise ValueError(f"Unknown mesh type: {mesh_cfg['type']}")

    # --- Basis Setup ---
    print(f"Initializing Basis (Poly Degree N={polynomial_degree})...")
    basis = Basis(polynomial_degree=polynomial_degree, device=device)

    # --- Solver Setup ---
    ic_name = ic_cfg['name']
    print(f"Initializing Solver with IC: {ic_name}...")
    ic_func = get_initial_condition_func(ic_name)
    
    # Instantiate Euler2DSolver
    solver = Euler2DSolver(mesh, basis, config=cfg)
    solver.initialize(ic_func)

    # --- Output Directory ---
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    writer = HDF5Writer(os.path.join(output_dir, "results.h5"), mesh)

    # Save Initial State
    print("Saving initial state...")
    writer.write_step(0, 0.0, solver.Q.numpy())

    # --- Run Simulation ---
    print(f"Starting simulation on device '{device}'...")
    print(f"t_final = {t_final}, CFL = {CFL}")
    
    # Instantiate Driver
    driver = TimeIntegrator(solver)
    max_steps = sim_cfg.get('max_steps', None)
    driver.solve(t_final=t_final, CFL=CFL, log_frequency=log_freq, writer=writer, max_steps=max_steps)
    
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
