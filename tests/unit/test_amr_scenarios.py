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
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)

    qt.uniform_refine(level=1, state=state, basis=basis)
    
    h_codes = qt.block_morton_codes.numpy()
    h_levels = qt.block_levels.numpy()
    
    h_levels[0] = 2
    h_codes[0] = 16 # (1 << 4) | interleaved(0,0)
    
    qt.block_levels = wp.array(h_levels, dtype=wp.int32, device=device)
    qt.block_morton_codes = wp.array(h_codes, dtype=wp.int32, device=device)
    
    try:
        qt.build_connectivity(state)
    except Exception as e:
        pytest.fail(f"build_connectivity failed with mixed levels: {e}")

    map_keys = qt.map_keys.numpy()
    map_vals = qt.map_values.numpy()
    
    target_key = 16
    found = False
    for k, v in zip(map_keys, map_vals):
        if k == target_key and v == 0:
            found = True
            break
    
    assert found, "Did not find Block 0 with code 16 in hash map"

def test_hanging_node_connectivity(device):
    """
    Test 2: Hanging Node Connectivity
    Verify Fine -> Coarse and Coarse -> Fine neighbor lookups.
    """
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    qt.uniform_refine(level=1, state=state, basis=basis)
    # BL=4, BR=5, TL=6, TR=7
    
    h_codes = qt.block_morton_codes.numpy()
    active = state.active_block_indices.numpy()[:qt.num_blocks]
    
    idx_to_refine = -1
    code_l1_bl = morton_encode(0, 0, 1) # code 4
    for idx in active:
        if h_codes[idx] == code_l1_bl:
            idx_to_refine = idx
            break
            
    assert idx_to_refine != -1
    qt.refine_blocks([idx_to_refine], state, basis)
    
    active = state.active_block_indices.numpy()[:qt.num_blocks]
    codes = qt.block_morton_codes.numpy()
    levels = qt.block_levels.numpy()
    neighbors = state.neighbors.numpy()
    
    # L1 Block (Bottom-Right). Code: (1,0) at L1 -> 5.
    code_l1_br = morton_encode(1, 0, 1)
    l1_br_idx = -1
    for idx in active:
        if levels[idx] == 1 and codes[idx] == code_l1_br:
            l1_br_idx = idx
            break
    assert l1_br_idx != -1
    
    # L2 Block (Top-Right of BL). (ix=1, iy=1) at L2 -> 19.
    code_l2_tr = morton_encode(1, 1, 2)
    l2_tr_idx = -1
    for idx in active:
        if levels[idx] == 2 and codes[idx] == code_l2_tr:
            l2_tr_idx = idx
            break
    assert l2_tr_idx != -1
    
    # L2 TR (1,1). RIGHT (Face 1) -> (2,1) L2 -> MORTAR interface.
    assert neighbors[l2_tr_idx, 1] == -2
        
    # L1 BR (1,0). LEFT (Face 3) -> (0,0) L1 -> MORTAR interface.
    assert neighbors[l1_br_idx, 3] == -2
        
    # L2 TR (1,1). BOTTOM (Face 0) -> (1,0) L2 -> index 17.
    code_l2_br = morton_encode(1, 0, 2) # 17
    l2_br_idx = -1
    for idx in active:
        if levels[idx] == 2 and codes[idx] == code_l2_br:
            l2_br_idx = idx
            break
            
    assert neighbors[l2_tr_idx, 0] == l2_br_idx