import pytest
import warp as wp
import numpy as np
from src.geometry.quadtree import Quadtree
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.kernels.grid_kernels import morton_decode

@pytest.fixture(scope="module")
def device():
    wp.init()
    return "cpu"

def test_quadtree_refine_blocks(device):
    """
    Test refining a block in the Quadtree and checking topology updates.
    """
    # Setup
    Np = 4
    # Basis constructor: polynomial_degree, device, dtype, over_integration_order
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    # 1. Uniform Refine Level 1 (2x2 = 4 blocks)
    # Codes: 0(0,0), 1(1,0), 2(0,1), 3(1,1)
    qt.uniform_refine(level=1, state=state, basis=basis)
    
    assert qt.num_blocks == 4
    active_indices = state.active_block_indices.numpy()[:4]
    # In uniform refine, indices are 0,1,2,3
    assert sorted(active_indices) == [0, 1, 2, 3]
    
    # 2. Refine Block 0 (Indices: 0)
    # This should deactivate 0, and add 4 children (Indices 4, 5, 6, 7 probably)
    qt.refine_blocks([0], state, basis)
    
    # Check counts
    # 4 (initial) - 1 (removed) + 4 (added) = 7
    assert qt.num_blocks == 7
    
    active_indices = state.active_block_indices.numpy()[:7]
    active_set = set(active_indices)
    
    assert 0 not in active_set # Parent removed
    assert {1, 2, 3}.issubset(active_set) # Others remain
    
    # Check levels
    levels = qt.block_levels.numpy()
    # Indices in active_set that are NOT 1,2,3 must be level 2
    new_blocks = active_set - {1, 2, 3}
    assert len(new_blocks) == 4
    for idx in new_blocks:
        assert levels[idx] == 2
        
    # Check Morton Codes of new blocks
    # Parent 0 was code 0 (0,0) at level 1.
    # Children at level 2:
    # (0,0)->0, (1,0)->1, (0,1)->2, (1,1)->3
    # Wait, parent (0,0) at level 1 covers (0,0) to (1,1) at level 2?
    # No. 
    # Level 0: 1 block (0,0)
    # Level 1: 4 blocks. Parent (0,0) covers (0,0) of level 0.
    #   Coords at L1: (0,0), (1,0), (0,1), (1,1).
    #   Block 0 is (0,0) at L1.
    # Level 2: Block 0 splits into 4.
    #   Coords at L2: (0,0), (1,0), (0,1), (1,1).
    #   (These are global coords at L2).
    #   Wait, L1(0,0) -> L2 coords: 2*0 + [0,1], ...
    #   So L2 coords are (0,0), (1,0), (0,1), (1,1).
    #   Morton codes: 0, 1, 2, 3.
    #   So children codes are 0, 1, 2, 3 at Level 2.
    
    codes = qt.block_morton_codes.numpy()
    child_codes = [codes[idx] for idx in new_blocks]
    assert sorted(child_codes) == [0, 1, 2, 3]
    
    # 3. Check Connectivity (AMR Neighbor Lookup)
    # Let's check a child block. 
    # Child with code 1 (x=1, y=0 at L2).
    # Its Right neighbor is x=2.
    # At L2, x=2, y=0 is part of Block 1 (Parent was x=1, y=0 at L1 -> x=2,3 at L2).
    # So neighbor is Block 1.
    
    # Find the pool index of child 1
    child_1_idx = -1
    for idx in new_blocks:
        if codes[idx] == 1:
            child_1_idx = idx
            break
            
    assert child_1_idx != -1
    
    neighbors = state.neighbors.numpy()
    # Face 1 is Right
    right_n = neighbors[child_1_idx, 1]
    
    # Expectation: right_n should be index of Block 1
    # Block 1 was index 1 in the pool (preserved).
    assert right_n == 1 
    
    # Check Child 2 (x=0, y=1 at L2) looking Top (Face 3)
    # Neighbor y=2. Corresponds to Block 2 (Parent x=0, y=1 at L1 -> y=2,3 at L2).
    child_2_idx = -1
    for idx in new_blocks:
        if codes[idx] == 2:
            child_2_idx = idx
            break
            
    top_n = neighbors[child_2_idx, 3]
    assert top_n == 2

    print("AMR Connectivity Test Passed!")
