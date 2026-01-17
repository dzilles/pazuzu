import pytest
import numpy as np
import warp as wp
from src.core.config import PazuzuConfig
from solver import PazuzuSolver

# Ensure Warp is initialized
wp.init()

class TestIBMAMRPhysics:
    def test_geometry_induced_refinement(self):
        """
        Verifies that enabling IBM forces refinement around the interface,
        even when gradient-based refinement is effectively disabled.
        """
        # 1. Config
        # Center (0,0), Radius 0.5. Domain [-2, 2].
        # Max Depth 5.
        # High refine_threshold to prevent gradient AMR.
        config_dict = {
            "case_name": "test_ibm_amr",
            "io": {"write_interval": 100},
            "mesh": {"x_min": -2.0, "x_max": 2.0, "y_min": -2.0, "y_max": 2.0, "initial_depth": 2},
            "amr": {
                "enabled": True,
                "max_depth": 5,
                "max_blocks": 5000,
                "refinement_threshold": 1000.0, # Very high -> No flow refinement
                "coarsening_threshold": 500.0,
                "refine_interval": 1 # Enable dynamic AMR
            },
            "ibm": {
                "enabled": True,
                "mode": "analytical",
                "geometric_params": {"center_x": 0.0, "center_y": 0.0, "radius": 0.5}
            },
            "initial_condition": {"name": "vortex", "params": {"beta": 0.0}} # Stagnant flow
        }
        config = PazuzuConfig(**config_dict)
        
        # 2. Init Solver
        solver = PazuzuSolver(config)
        
        # Initial state: Should be level 2 (uniform)
        # Check num blocks
        initial_blocks = int(solver.quadtree.num_blocks) # Assuming int/scalar
        if hasattr(solver.quadtree.num_blocks, "numpy"):
             initial_blocks = int(solver.quadtree.num_blocks.numpy()[0])
        
        print(f"Initial Blocks (Level 2): {initial_blocks}")
        
        # 3. Action: Force Adaptation
        # The solver does initial adaptation in __init__, but thresholds might have prevented it?
        # Wait, in __init__, it calls adapt_mesh.
        # If refine_threshold is 1000.0, gradient check finds nothing.
        # BUT ibm_enabled=True should trigger geometry refinement.
        # So adaptation should have ALREADY happened in __init__.
        
        # Let's verify the current state.
        
        # Access data
        x = solver.state.x.numpy()
        y = solver.state.y.numpy()
        phi = solver.state.phi.numpy()
        levels = solver.quadtree.block_levels.numpy()
        active = solver.state.active_block_indices.numpy()[:initial_blocks]
        
        # Check blocks crossing interface
        refined_count = 0
        interface_blocks = 0
        
        for idx in active:
            # Check if block crosses interface
            # min_phi * max_phi <= 0
            p_vals = phi[idx, :]
            min_p = np.min(p_vals)
            max_p = np.max(p_vals)
            
            if min_p * max_p <= 0.0:
                interface_blocks += 1
                # Should be max_depth (5)
                # Note: It might take multiple passes to reach max_depth if adapt_mesh only refines once per call?
                # adapt_mesh calls refine_blocks once.
                # So if we start at level 2, one call -> level 3.
                # In __init__, we only call adapt_mesh once.
                # So we expect level 3, not necessarily 5, unless we loop.
                # The Requirement says: "Verify that blocks... have level == 5".
                # This implies I should call adapt_mesh multiple times until convergence.
                pass
        
        # Let's run adapt_mesh loop until max level reached or no change
        print("Running adaptation loop...")
        for i in range(10):
            prev_blocks = initial_blocks
            solver.quadtree.adapt_mesh(
                solver.state, 
                solver.basis, 
                config.amr.refinement_threshold, 
                config.amr.coarsening_threshold, 
                ibm_enabled=True
            )
            # Re-init SDF because refinement creates new blocks with zero/garbage phi
            solver._init_ibm_sdf()
            
            curr_blocks = solver.quadtree.num_blocks
            if hasattr(curr_blocks, "numpy"):
                curr_blocks = int(curr_blocks.numpy()[0])
            
            print(f"Adaptation Iteration {i+1}: {curr_blocks} blocks")
            if curr_blocks == prev_blocks:
                break
            initial_blocks = curr_blocks

        # 4. Final Verification
        active = solver.state.active_block_indices.numpy()[:initial_blocks]
        phi = solver.state.phi.numpy()
        levels = solver.quadtree.block_levels.numpy()
        
        max_level_found = 0
        interface_blocks_checked = 0
        
        for idx in active:
            p_vals = phi[idx, :]
            min_p = np.min(p_vals)
            max_p = np.max(p_vals)
            
            lvl = levels[idx]
            if lvl > max_level_found:
                max_level_found = lvl
            
            if min_p * max_p <= 0.0:
                # Interface Block
                interface_blocks_checked += 1
                assert lvl == 5, f"Interface block {idx} at level {lvl}, expected 5"
            else:
                # Far field. Might be refined due to balancing (2:1 constraint).
                # But generally lower level.
                pass

        print(f"Max Level Found: {max_level_found}")
        print(f"Interface Blocks Verified: {interface_blocks_checked}")
        assert max_level_found == 5
        assert interface_blocks_checked > 0

    def test_full_simulation_step(self):
        """
        Run a few time steps to ensure robustness.
        """
        config_dict = {
            "case_name": "test_ibm_run",
            "io": {"write_interval": 10},
            "simulation": {"max_steps": 5, "t_final": 0.05},
            "mesh": {"x_min": -2.0, "x_max": 2.0, "y_min": -2.0, "y_max": 2.0, "initial_depth": 2},
            "numerics": {"dt_static": 0.01}, # Force small steps
            "ibm": {
                "enabled": True,
                "mode": "analytical",
                "geometric_params": {"center_x": 0.0, "center_y": 0.0, "radius": 0.5}
            },
            # Uniform flow
            "physics": {"u_inf": 1.0, "rho_inf": 1.0, "p_inf": 1.0},
            "initial_condition": {"name": "vortex", "params": {"beta": 0.0}}
        }
        config = PazuzuConfig(**config_dict)
        solver = PazuzuSolver(config)
        
        # Run
        solver.run()
        
        # Check
        assert solver.state.step == 5
        
        # Check for NaNs
        q = solver.state.q.numpy()
        
        num_blocks = solver.quadtree.num_blocks
        if hasattr(num_blocks, "numpy"):
            num_blocks = int(num_blocks.numpy()[0])
            
        active_indices = solver.state.active_block_indices.numpy()[:num_blocks]
        
        for idx in active_indices:
            block_q = q[idx]
            assert not np.any(np.isnan(block_q)), f"NaN detected in block {idx}"

