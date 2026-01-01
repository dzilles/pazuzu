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

def compute_total_mass(state, basis, qt):
    """ Helper to compute the integral of density across the mesh. """
    active = state.active_block_indices.numpy()[:qt.num_blocks]
    q = state.q.numpy()
    levels = qt.block_levels.numpy()
    
    # 2D weights for reference element [-1, 1]^2
    w1d = basis.weights_1d.numpy()
    w2d = np.kron(w1d, w1d)
    
    total_mass = 0.0
    # Assuming ROOT_BOUNDS are (-1, -1, 1, 1) -> L=2
    L = 2.0 
    
    for idx in active:
        level = levels[idx]
        # Jacobian for Cartesian grid: J = (dx/2) * (dy/2)
        # dx = L / 2^level
        J = (L / (2.0**(level + 1)))**2
        
        block_mass = np.sum(w2d * q[idx, :, 0]) * J
        total_mass += block_mass
        
    return total_mass

def test_amr_conservation_of_mass(device):
    """
    Test 1: Conservation of Mass (Refine -> Coarsen)
    Ensures that the integral of the state is preserved across AMR operations.
    """
    # Setup
    poly_degree = 2
    basis = Basis(polynomial_degree=poly_degree, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    # 1. Initialize 2x2 grid (Level 1)
    qt.uniform_refine(level=1, state=state, basis=basis)
    
    # Set uniform density rho=1.0
    q_np = np.zeros((100, basis.Np, 4), dtype=np.float32)
    q_np[:, :, 0] = 1.0
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    # Total mass should be 4.0 (Area 2x2 * rho 1.0)
    initial_mass = compute_total_mass(state, basis, qt)
    assert pytest.approx(initial_mass) == 4.0
    
    # 2. Refine one block (Block 0)
    # Use refine_marked_blocks by setting a mock gradient
    q_np[0, 0, 0] = 10.0 # Force gradient in block 0
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    # RECOMPUTE mass after the spike to have the correct baseline for conservation
    initial_mass_with_spike = compute_total_mass(state, basis, qt)
    
    qt.refine_marked_blocks(state, basis, threshold=0.5)
    
    # Verify counts
    assert qt.num_blocks == 7
    
    # Check mass conservation
    refined_mass = compute_total_mass(state, basis, qt)
    assert pytest.approx(refined_mass, abs=1e-6) == initial_mass_with_spike
    
    # 3. Coarsen back
    qt.coarsen_marked_blocks(state, basis)
    
    # Verify counts
    assert qt.num_blocks == 4
    
    # Check mass conservation
    coarsened_mass = compute_total_mass(state, basis, qt)
    assert pytest.approx(coarsened_mass, abs=1e-6) == initial_mass_with_spike


def test_amr_data_preservation_linear(device):
    """
    Test 2: Data Preservation
    Ensures that a linear function is perfectly preserved during refinement.
    """
    # Setup
    poly_degree = 2 # P=2 supports linear exactly
    basis = Basis(polynomial_degree=poly_degree, device=device)
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    # 1. Initialize grid (L1)
    qt.uniform_refine(level=1, state=state, basis=basis)
    
    # Set f(x,y) = 10*x + 10*y
    # This has a high gradient (diff=20 across L1 block) so it triggers refinement naturally
    x_np = state.x.numpy()
    y_np = state.y.numpy()
    q_np = np.zeros((100, basis.Np, 4), dtype=np.float32)
    q_np[:, :, 0] = 10.0 * (x_np + y_np)
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    # 2. Refine Blocks
    # Threshold 1.0 will trigger refinement for all blocks
    qt.refine_marked_blocks(state, basis, threshold=1.0)
    
    # 3. Verification
    # Evaluate children: q should match 10*(x + y) at node locations
    q_refined = state.q.numpy()
    x_refined = state.x.numpy()
    y_refined = state.y.numpy()
    levels = qt.block_levels.numpy()
    active = state.active_block_indices.numpy()[:qt.num_blocks]
    
    for idx in active:
        if levels[idx] == 2:
            # This is a new child
            expected = 10.0 * (x_refined[idx, :] + y_refined[idx, :])
            actual = q_refined[idx, :, 0]
            
            np.testing.assert_allclose(actual, expected, atol=1e-6)
            
    print("\nAMR Data Preservation Test Passed!")
