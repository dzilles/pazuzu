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

def test_quadtree_2to1_balancing(device):
    """
    Verify that the 2:1 balancing constraint works.
    If we have a uniform grid at level 1 (2x2), and we refine one block to level 2,
    and then we mark a level 2 block for level 3, its level 1 neighbors must
    first be refined to level 2.
    """
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=1000)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=1000, basis=basis)
    
    # 1. Start with Level 1 (4 blocks)
    qt.uniform_refine(level=1, state=state, basis=basis)
    assert qt.num_blocks == 4
    
    # Identify block (0,0) at Level 1 (Morton code 4)
    codes = qt.block_morton_codes.numpy()
    target_idx = -1
    for i in range(qt.num_blocks):
        p_idx = state.active_block_indices.numpy()[i]
        if codes[p_idx] == 4:
            target_idx = p_idx
            break
    assert target_idx != -1
    
    # 2. Refine block (0,0) to Level 2
    qt.refine_blocks([target_idx], state, basis)
    assert qt.num_blocks == 7 # 4 - 1 + 4 = 7
    
    # Now we have:
    # 4 children at Level 2 (from block 0,0)
    # 3 blocks at Level 1
    
    # Identify a Level 2 child, e.g., (0,0) at L2 (Morton code 16)
    codes = qt.block_morton_codes.numpy()
    levels = qt.block_levels.numpy()
    l2_target_idx = -1
    for i in range(qt.num_blocks):
        p_idx = state.active_block_indices.numpy()[i]
        if levels[p_idx] == 2 and codes[p_idx] == 16:
            l2_target_idx = p_idx
            break
    assert l2_target_idx != -1
    
    # 3. Mark the L2 child for refinement to L3
    # We'll use a threshold that only marks this specific block.
    # Instead of calling refine_marked_blocks directly with a threshold (which depends on q),
    # we'll mock the gradient marking.
    
    # Force q values to trigger refinement if we used mark_blocks_gradient, 
    # but it's easier to just call the balancing kernel manually or 
    # slightly modify the test to use a dummy gradient.
    
    # Let's mock the refine_flags
    refine_flags = wp.zeros(qt.max_blocks, dtype=wp.int32, device=device)
    h_flags = np.zeros(qt.max_blocks, dtype=np.int32)
    h_flags[l2_target_idx] = 1
    wp.copy(refine_flags, wp.array(h_flags, dtype=wp.int32, device=device))
    
    from src.kernels.amr_kernels import balance_refine_flags
    changed = wp.zeros(1, dtype=wp.int32, device=device)
    h_changed = np.array([1], dtype=np.int32)
    wp_h_changed = wp.from_numpy(h_changed, dtype=wp.int32, device=device)
    
    while h_changed[0] > 0:
        changed.zero_()
        wp.launch(
            kernel=balance_refine_flags,
            dim=qt.num_blocks,
            inputs=[
                state.active_block_indices,
                qt.num_blocks,
                qt.block_morton_codes,
                qt.block_levels,
                qt.map_keys,
                qt.map_values,
                qt.map_capacity,
                0, 0, # periodic_x, periodic_y
                refine_flags,
                changed
            ],
            device=device
        )
        wp.copy(wp_h_changed, changed)
        h_changed = wp_h_changed.numpy()
        
    # 5. Verify that level 1 neighbors are now marked
    res_flags = refine_flags.numpy()
    
    # Neighbors of (0,0) L2 are at (1,0) L2, (0,1) L2, and importantly Level 1 neighbors.
    # The Level 1 blocks were (1,0), (0,1), (1,1) at L1.
    # Coords of neighbors at Level 1 are:
    # Right of (0,0) L2 is (1,0) L2. Its parent is (0,0) L1 (already refined).
    # Wait, if we are at (0,0) L2, our neighbors at SAME level are (1,0) L2 and (0,1) L2.
    # If we are at (1,1) L2 (code 19), our neighbors are (2,1) L2 and (1,2) L2.
    # (2,1) L2 is inside block (1,0) L1.
    # So if we refine (1,1) L2, then block (1,0) L1 MUST be refined first.
    
    l2_tr_child = -1
    for i in range(qt.num_blocks):
        p_idx = state.active_block_indices.numpy()[i]
        if levels[p_idx] == 2 and codes[p_idx] == 19: # (1,1) L2
            l2_tr_child = p_idx
            break
    assert l2_tr_child != -1
    
    # Reset flags and mark (1,1) L2
    h_flags.fill(0)
    h_flags[l2_tr_child] = 1
    wp.copy(refine_flags, wp.array(h_flags, dtype=wp.int32, device=device))
    
    h_changed[0] = 1
    while h_changed[0] > 0:
        changed.zero_()
        wp.launch(balance_refine_flags, dim=qt.num_blocks, 
                  inputs=[state.active_block_indices, qt.num_blocks, qt.block_morton_codes, qt.block_levels, qt.map_keys, qt.map_values, qt.map_capacity, 0, 0, refine_flags, changed],
                  device=device)
        wp.copy(wp_h_changed, changed)
        h_changed = wp_h_changed.numpy()
        
    res_flags = refine_flags.numpy()
    
    # Find Level 1 neighbors
    l1_right_idx = -1
    l1_top_idx = -1
    for i in range(qt.num_blocks):
        p_idx = state.active_block_indices.numpy()[i]
        if levels[p_idx] == 1:
            if codes[p_idx] == 5: l1_right_idx = p_idx # (1,0) L1
            if codes[p_idx] == 6: l1_top_idx = p_idx   # (0,1) L1
            
    assert l1_right_idx != -1
    assert l1_top_idx != -1
    
    assert res_flags[l1_right_idx] == 1, "Level 1 Right neighbor should be marked for balancing"
    assert res_flags[l1_top_idx] == 1, "Level 1 Top neighbor should be marked for balancing"
    print("SUCCESS: 2:1 balancing correctly marked level 1 neighbors.")
