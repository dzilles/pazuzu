import pytest
import warp as wp
import numpy as np
from src.geometry.quadtree import Quadtree
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.kernels.grid_kernels import morton_decode, morton_encode

@pytest.fixture(scope="module")
def device():
    wp.init()
    return "cpu"

def test_hierarchy_uniqueness(device):
    """
    Test 1: Hierarchy Uniqueness
    Ensure hash map handles blocks at different levels without collision.
    """
    # Setup
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)

    # Action: Create a uniform Level 1 grid.
    qt.uniform_refine(level=1, state=state, basis=basis)
    # Uniform L1 (2x2) -> 4 blocks. Codes 4, 5, 6, 7. Active indices 0,1,2,3.
    
    # Manually force one block to Level 2
    # Let's pick Block 0 (Code 4, Level 1).
    # We change it to Code 16, Level 2.
    # Code 16 corresponds to (0,0) at L2.
    # Geometrically these overlap (top-left corner of L1 block is top-left corner of L2 block).
    
    h_codes = qt.block_morton_codes.numpy()
    h_levels = qt.block_levels.numpy()
    
    # Modify Block 0
    h_levels[0] = 2
    h_codes[0] = 16 # (1 << 4) | interleaved(0,0)
    
    # Upload back
    qt.block_levels = wp.array(h_levels, dtype=wp.int32, device=device)
    qt.block_morton_codes = wp.array(h_codes, dtype=wp.int32, device=device)
    
    # Check: Rebuild connectivity
    try:
        qt.build_connectivity(state)
    except Exception as e:
        pytest.fail(f"build_connectivity failed with mixed levels: {e}")

    # Verify we can find the block
    # Look up Block 0 (Code 16, Level 2)
    # This requires manual kernel launch or inspecting the map.
    # Let's inspect the map on host.
    
    map_keys = qt.map_keys.numpy()
    map_vals = qt.map_values.numpy()
    
    target_key = 16
    found = False
    for k, v in zip(map_keys, map_vals):
        if k == target_key and v == 0:
            found = True
            break
    
    assert found, "Did not find Block 0 with code 16 in hash map"
    print("\nTest 1: Hierarchy Uniqueness Passed")


def test_hanging_node_connectivity(device):
    """
    Test 2: Hanging Node Connectivity
    Verify Fine -> Coarse and Coarse -> Fine neighbor lookups.
    """
    # Setup
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    # 1. Start with Root -> Refine to Level 1
    qt.uniform_refine(level=1, state=state, basis=basis)
    # Blocks: 0(BL), 1(BR), 2(TL), 3(TR) at L1.
    
    # 2. Refine Child 0 (Bottom-Left) to Level 2
    # Block 0 is the one with Code 4.
    codes = qt.block_morton_codes.numpy()
    active = state.active_block_indices.numpy()[:qt.num_blocks]
    
    idx_to_refine = -1
    code_l1_bl = morton_encode(0, 0, 1) # code 4
    for idx in active:
        if codes[idx] == code_l1_bl:
            idx_to_refine = idx
            break
            
    assert idx_to_refine != -1
    
    qt.refine_blocks([idx_to_refine], state, basis)
    
    # Update host views
    active = state.active_block_indices.numpy()[:qt.num_blocks]
    codes = qt.block_morton_codes.numpy()
    levels = qt.block_levels.numpy()
    neighbors = state.neighbors.numpy()
    
    # Identify Blocks
    # L1 Blocks: 1, 2, 3 should still be there.
    # L2 Blocks: 4 children of 0.
    
    # Find L1 Block (Bottom-Right). Code: (1,0) at L1.
    # Encode(1,0, 1) -> 5.
    code_l1_br = morton_encode(1, 0, 1)
    l1_br_idx = -1
    for idx in active:
        if levels[idx] == 1 and codes[idx] == code_l1_br:
            l1_br_idx = idx
            break
    assert l1_br_idx != -1, "Could not find L1 Bottom-Right block"
    
    # Find L2 Block (Top-Right grandchild of Bottom-Left).
    # Parent (0,0) L1.
    # L2 Coords: x=1, y=1.
    # Encode(1,1, 2) -> (1 << 4) | interleaved(1,1) = 16 | 3 = 19.
    code_l2_tr = morton_encode(1, 1, 2)
    l2_tr_idx = -1
    for idx in active:
        if levels[idx] == 2 and codes[idx] == code_l2_tr:
            l2_tr_idx = idx
            break
    assert l2_tr_idx != -1, "Could not find L2 Top-Right block"
    
    # --- Verification 1: Fine -> Coarse ---
    # Pick L2 TR block. Look Right (Face 1).
    # Should find MORTAR_FLAG (-2) because it's a hierarchical interface.
    
    neighbor_right = neighbors[l2_tr_idx, 1]
    assert neighbor_right == -2, \
        f"Fine->Coarse Failed: L2 Block {l2_tr_idx} Right neighbor is {neighbor_right}, expected -2 (MORTAR_FLAG)"
        
    # --- Verification 2: Coarse -> Fine ---
    # Pick L1 BR block. Look Left (Face 0).
    # Should be -2 (MORTAR_FLAG) now that we implemented mark_mortar_neighbors!
    
    neighbor_left = neighbors[l1_br_idx, 0]
    assert neighbor_left == -2, \
        f"Coarse->Fine Failed: L1 Block {l1_br_idx} Left neighbor is {neighbor_left}, expected -2 (MORTAR_FLAG)"
        
    # --- Verification 3: Same Level ---
    # L2 TR (x=1, y=1). Look Bottom (Face 2) -> (x=1, y=0).
    # Should be L2 BR block (Code 17 at L2).
    
    code_l2_br = morton_encode(1, 0, 2) # 16 | 1 = 17
    l2_br_idx = -1
    for idx in active:
        if levels[idx] == 2 and codes[idx] == code_l2_br:
            l2_br_idx = idx
            break
            
    neighbor_bottom = neighbors[l2_tr_idx, 2]
    assert neighbor_bottom == l2_br_idx, \
        f"Same Level Failed: L2 Block {l2_tr_idx} Bottom neighbor is {neighbor_bottom}, expected {l2_br_idx}"

    print("\nTest 2: Hanging Node Connectivity Passed")
