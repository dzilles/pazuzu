import pytest
import warp as wp
import numpy as np
from src.geometry.quadtree import Quadtree
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.core.config import PhysicsConfig # Need params for mortar kernel
from src.kernels.mortar_kernels import compute_mortar_fluxes

@pytest.fixture(scope="module")
def device():
    wp.init()
    return "cpu"

def test_mortar_list_population(device):
    """
    Verify that mortar interfaces are correctly identified and stored.
    """
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    # 1. Uniform L1 (4 blocks)
    qt.uniform_refine(level=1, state=state, basis=basis)
    
    # 2. Refine Block 0 (Bottom-Left)
    # 2. Refine Block 0 (Bottom-Left)
    # Set gradient to trigger refinement
    q_np = np.zeros((100, Np, 4), dtype=np.float32)
    q_np[0, 0, 0] = 10.0
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    qt.refine_marked_blocks(state, basis, threshold=0.5)
    
    # Grid: 
    # BL sector: 4 small blocks (L2).
    # BR, TL, TR sectors: 1 large block each (L1).
    #
    # Interfaces:
    # L2 TR child (x=1, y=1) -> Right neighbor is L1 BR block.
    # L2 TR child (x=1, y=1) -> Top neighbor is L1 TL block.
    # etc.
    
    # Check Mortar List
    num_mortars = qt.num_mortars.numpy()[0]
    assert num_mortars > 0, "No mortars detected!"
    
    mortar_list = qt.mortar_list.numpy()[:num_mortars]
    
    # Find L2 TR child (Code 3 at L2).
    codes = qt.block_morton_codes.numpy()
    levels = qt.block_levels.numpy()
    
    l2_tr_idx = -1
    for idx in range(qt.num_blocks):
        # We need active blocks
        pool_idx = state.active_block_indices.numpy()[idx]
        if levels[pool_idx] == 2 and codes[pool_idx] == 3: # (1,1) child
            l2_tr_idx = pool_idx
            break
            
    assert l2_tr_idx != -1
    
    # Check if this block is in mortar list
    # It should have a Right neighbor (Face 1) which is Coarse.
    # And a Top neighbor (Face 3) which is Coarse.
    
    found_right = False
    found_top = False
    
    for i in range(num_mortars):
        fine_idx = mortar_list[i, 0]
        face = mortar_list[i, 1]
        coarse_idx = mortar_list[i, 2]
        subface = mortar_list[i, 3]
        
        if fine_idx == l2_tr_idx:
            if face == 1: # Right
                found_right = True
                # Subface: y=1 (top half) -> 1
                assert subface == 1
            if face == 3: # Top
                found_top = True
                # Subface: x=1 (right half) -> 1
                assert subface == 1
                
    assert found_right, "Mortar interface Right not found for L2 TR block"
    assert found_top, "Mortar interface Top not found for L2 TR block"

def test_mortar_flux_kernel_execution(device):
    """
    Verify that the flux kernel runs.
    """
    # Setup from previous test state
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    qt.uniform_refine(level=1, state=state, basis=basis)
    # 2. Refine Block 0 (Bottom-Left)
    # Set gradient to trigger refinement
    q_np = np.zeros((100, Np, 4), dtype=np.float32)
    q_np[0, 0, 0] = 10.0
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    qt.refine_marked_blocks(state, basis, threshold=0.5)
    
    # Params
    params = PhysicsConfig()
    
    # Launch Kernel
    wp.launch(
        kernel=compute_mortar_fluxes,
        dim=qt.num_mortars.numpy()[0],
        inputs=[
            state.q,
            state.rhs,
            qt.mortar_list,
            qt.num_mortars,
            basis.face_nodes,
            basis.P_left,
            basis.P_right,
            basis.R_left,
            basis.R_right,
            params
        ],
        device=device
    )
    
    # Check that RHS is not NaN?
    # Or just that it finished.
    rhs = state.rhs.numpy()
    assert not np.any(np.isnan(rhs))
