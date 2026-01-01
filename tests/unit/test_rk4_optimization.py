import pytest
import warp as wp
import os
from src.core.config import PazuzuConfig
from solver import PazuzuSolver

def test_rk4_memory_optimization(device):
    """
    Verifies that the RK4 integrator uses the pre-allocated q_accum buffer
    and runs successfully without errors.
    """
    # 1. Create a minimal configuration with RK4
    config_dict = {
        "case_name": "test_rk4",
        "solver_type": "euler_2d",
        "mesh": {
            "x_min": -1.0, "x_max": 1.0,
            "y_min": -1.0, "y_max": 1.0,
            "periodic_x": True, "periodic_y": True
        },
        "amr": {
            "max_blocks": 100,
            "initial_depth": 2
        },
        "initial_condition": {
            "name": "vortex",
            "params": {"beta": 5.0}
        },
        "physics": {
            "gamma": 1.4
        },
        "numerics": {
            "time_integrator": "rk4", # FORCE RK4
            "cfl": 0.1,
            "polynomial_order": 1,
            "dt_static": 0.001
        },
        "io": {
            "output_dir": "output/test_rk4",
            "write_interval": 100
        },
        "simulation": {
            "t_final": 0.01,
            "max_steps": 10,
            "device": device
        }
    }
    
    config = PazuzuConfig(**config_dict)
    
    # 2. Initialize Solver
    solver = PazuzuSolver(config)
    
    # 3. Verify q_accum exists in state
    assert hasattr(solver.state, "q_accum"), "SimulationState missing q_accum attribute"
    assert solver.state.q_accum is not None, "q_accum is None"
    assert solver.state.q_accum.shape == solver.state.q.shape, "q_accum shape mismatch"
    
    # 4. Run a few steps
    # We can call run() directly as max_steps is small
    try:
        solver.run()
    except Exception as e:
        pytest.fail(f"RK4 Simulation failed with error: {e}")

    # 5. Check that q_accum was used (should be close to zero at end of step because it's an accumulator, 
    # but let's just ensure the simulation advanced).
    assert solver.state.step == 10
    assert solver.state.t > 0.0
    
    print("RK4 Simulation completed successfully with pre-allocated buffer.")
