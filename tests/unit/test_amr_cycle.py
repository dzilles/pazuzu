import pytest
import warp as wp
import numpy as np
from src.geometry.quadtree import Quadtree
from src.core.simulation_state import SimulationState
from src.core.basis import Basis

@pytest.fixture(scope="module")
def device():
    wp.init()
    return "cpu"

def test_amr_refine_coarsen_cycle(device):
    """
    Test full AMR cycle: Refine -> Check -> Coarsen -> Check.
    """
    # Setup
    poly_degree = 1 # P=1 -> 2x2 nodes = 4 nodes per element
    basis = Basis(polynomial_degree=poly_degree, device=device)
    Np = basis.Np # Should be 4
    
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    # 1. Uniform L1 (4 blocks)
    qt.uniform_refine(level=1, state=state, basis=basis)
    assert qt.num_blocks == 4
    
    # Initialize field with gradient
    # q = x coordinate (density)
    # Block 0 (BL): x < 0. Block 1 (BR): x > 0.
    # Gradient in Block 0: ~1. Gradient in Block 1: ~1.
    
    x_np = state.x.numpy()
    q_np = np.zeros((100, Np, 4), dtype=np.float32)
    
    # Set only Block 0 to have a gradient
    # Block 0 indices are active_indices[0] -> 0
    # Block 0 nodes
    q_np[0, :, 0] = x_np[0, :] # Linear gradient
    
    # Other blocks constant
    q_np[1, :, 0] = 5.0 
    q_np[2, :, 0] = 5.0
    q_np[3, :, 0] = 5.0
    
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    # 2. Refine (Mark Gradient > 0.5)
    # Block 0 has diff ~ 1.0 (from -1 to 0).
    # Others have diff 0.
    qt.refine_marked_blocks(state, basis, threshold=0.5)
    
    # Check
    # Block 0 should be gone. Replaced by 4 children.
    # Total blocks: 4 - 1 + 4 = 7.
    assert qt.num_blocks == 7, f"Expected 7 blocks, got {qt.num_blocks}"
    
    # Check levels
    levels = qt.block_levels.numpy()
    active = state.active_block_indices.numpy()[:7]
    
    l1_count = 0
    l2_count = 0
    for idx in active:
        if levels[idx] == 1: l1_count += 1
        elif levels[idx] == 2: l2_count += 1
        
    assert l1_count == 3
    assert l2_count == 4
    
    # 3. Coarsen
    # Calling coarsen_marked_blocks should collapse the 4 children back to 1 parent.
    # (Since we implemented 'coarsen all eligible' logic).
    
    qt.coarsen_marked_blocks(state, basis)
    
    # Check
    assert qt.num_blocks == 4, f"Expected 4 blocks after coarsening, got {qt.num_blocks}"
    
    # Verify we have 4 blocks at level 1
    levels = qt.block_levels.numpy()
    active = state.active_block_indices.numpy()[:4]
    
    for idx in active:
        assert levels[idx] == 1

    print("\nAMR Cycle Test Passed!")
