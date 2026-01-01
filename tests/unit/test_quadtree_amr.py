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
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    qt.uniform_refine(level=1, state=state, basis=basis)
    assert qt.num_blocks == 4
    
    qt.refine_blocks([0], state, basis)
    assert qt.num_blocks == 7
    
    active_indices = state.active_block_indices.numpy()[:7]
    active_set = set(active_indices)
    assert 0 not in active_set
    
    levels = qt.block_levels.numpy()
    new_blocks = active_set - {1, 2, 3}
    for idx in new_blocks:
        assert levels[idx] == 2
        
    codes = qt.block_morton_codes.numpy()
    child_codes = [codes[idx] for idx in new_blocks]
    assert sorted(child_codes) == [16, 17, 18, 19]
    
    # 3. Check Connectivity
    # Child with code 17 (x=1, y=0 at L2).
    child_17_idx = -1
    for idx in new_blocks:
        if codes[idx] == 17:
            child_17_idx = idx
            break
    assert child_17_idx != -1
    
    neighbors = state.neighbors.numpy()
    # Face 1 is Right (+x). Should be MORTAR_FLAG (-2).
    assert neighbors[child_17_idx, 1] == -2 
    
    # Child 18 (x=0, y=1 at L2). Top (Face 2) -> MORTAR_FLAG.
    child_18_idx = -1
    for idx in new_blocks:
        if codes[idx] == 18:
            child_18_idx = idx
            break
    assert neighbors[child_18_idx, 2] == -2