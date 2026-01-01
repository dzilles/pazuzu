import pytest
import warp as wp
import numpy as np
from src.geometry.quadtree import Quadtree
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.core.config import PhysicsConfig
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
    
    qt.uniform_refine(level=1, state=state, basis=basis)
    
    q_np = np.zeros((100, Np, 4), dtype=np.float32)
    q_np[0, 0, 0] = 10.0
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    qt.refine_marked_blocks(state, basis, threshold=0.5)
    
    num_mortars = qt.num_mortars.numpy()[0]
    assert num_mortars > 0, "No mortars detected!"
    
    mortar_list = qt.mortar_list.numpy()[:num_mortars]
    
    # Find L2 TR child (Code 19 at L2).
    codes = qt.block_morton_codes.numpy()
    levels = qt.block_levels.numpy()
    
    l2_tr_idx = -1
    for idx in range(qt.num_blocks):
        pool_idx = state.active_block_indices.numpy()[idx]
        if levels[pool_idx] == 2 and codes[pool_idx] == 19: # (1,1) child
            l2_tr_idx = pool_idx
            break
            
    assert l2_tr_idx != -1
    
    # Check if this block is in mortar list
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
                assert subface == 1 # Top half of right edge
            if face == 2: # Top
                found_top = True
                assert subface == 1 # Right half of top edge
                
    assert found_right, "Mortar interface Right not found for L2 TR block"
    assert found_top, "Mortar interface Top not found for L2 TR block"

def test_mortar_flux_kernel_execution(device):
    """
    Verify that the flux kernel runs.
    """
    Np = 4
    basis = Basis(polynomial_degree=Np-1, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    qt.uniform_refine(level=1, state=state, basis=basis)
    q_np = np.zeros((100, Np, 4), dtype=np.float32)
    q_np[0, 0, 0] = 10.0
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    qt.refine_marked_blocks(state, basis, threshold=0.5)
    
    from src.kernels.structs import EquationParams32
    params = EquationParams32()
    params.gamma = 1.4; params.rho_floor = 1e-5; params.p_floor = 1e-5
    params.one = 1.0; params.half = 0.5; params.mu = 0.0
    params.prandtl = 0.72; params.cp = 1.0; params.gas_constant = 1.0
    params.flux_type = 0 # Rusanov
    
    num_mortars_host = int(qt.num_mortars.numpy()[0])
    wp.launch(
        kernel=compute_mortar_fluxes,
        dim=num_mortars_host,
        inputs=[
            state.q, state.rhs, qt.mortar_list, qt.num_mortars,
            basis.nodes_1d, basis.face_nodes, basis.dg_L, basis.dg_R,
            basis.P_left, basis.P_right, basis.R_left, basis.R_right,
            qt.root_bounds_wp, qt.block_levels, params
        ],
        device=device
    )
    
    rhs = state.rhs.numpy()
    assert not np.any(np.isnan(rhs))