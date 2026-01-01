import pytest
import warp as wp
import numpy as np
from src.geometry.quadtree import Quadtree
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.kernels.structs import EquationParams32
from src.kernels.fr_kernels import compute_fr_update
from src.kernels.mortar_kernels import compute_mortar_fluxes

@pytest.fixture(scope="module")
def device():
    wp.init()
    return "cpu"

def get_params():
    params = EquationParams32()
    params.gamma = 1.4
    params.rho_floor = 1e-5
    params.p_floor = 1e-5
    params.one = 1.0
    params.half = 0.5
    params.mu = 0.0
    params.prandtl = 0.72
    params.cp = 1.0
    params.gas_constant = 1.0
    params.flux_type = 0 # Rusanov
    return params

def test_broken_linear_patch(device):
    """
    Test 1: Constant Flow Consistency
    Verify that a constant field results in zero RHS across mortars.
    """
    poly_degree = 1
    basis = Basis(polynomial_degree=poly_degree, device=device)
    Np = basis.Np
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    # 1. Setup Non-Uniform Grid
    qt.uniform_refine(level=1, state=state, basis=basis)
    qt.refine_blocks([1, 3], state, basis)
    
    # 2. Set Constant Field q = 2.0
    params = get_params()
    q_np = np.zeros((100, Np, 4), dtype=np.float32)
    q_np[:, :, 0] = 2.0 # rho
    q_np[:, :, 1] = 2.0 # rhou (u=1)
    q_np[:, :, 2] = 0.0 # rhov
    q_np[:, :, 3] = 1.0 / (params.gamma - 1.0) + 0.5 * 2.0 * (1.0**2) # E
    
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    # 3. Compute RHS
    state.rhs.zero_()
    
    wp.launch(
        kernel=compute_fr_update,
        dim=qt.num_blocks * Np,
        inputs=[state.q, state.active_block_indices, state.neighbors, qt.num_blocks, state.rhs,
                basis.nodes_1d, basis.D1D, basis.dg_L, basis.dg_R, qt.root_bounds_wp, qt.block_levels, params, 0.0],
        device=device
    )
    
    num_mortars = int(qt.num_mortars.numpy()[0])
    if num_mortars > 0:
        wp.launch(
            kernel=compute_mortar_fluxes,
            dim=num_mortars,
            inputs=[state.q, state.rhs, qt.mortar_list, qt.num_mortars, basis.nodes_1d, basis.face_nodes,
                    basis.dg_L, basis.dg_R, basis.P_left, basis.P_right, basis.R_left, basis.R_right,
                    qt.root_bounds_wp, qt.block_levels, params],
            device=device
        )
        
    # 4. Verification
    rhs_np = state.rhs.numpy()
    active = state.active_block_indices.numpy()[:qt.num_blocks]
    
    for pool_idx in active:
        np.testing.assert_allclose(rhs_np[pool_idx, :, 0], 0.0, atol=1e-6)

    print("\nTest 1: Constant Flow Consistency Passed")


def test_mortar_conservation(device):
    """
    Test 2: Conservation Check
    Sum(rho * dA) should be preserved across mortar interfaces.
    """
    poly_degree = 1
    basis = Basis(polynomial_degree=poly_degree, device=device)
    Np = basis.Np
    qt = Quadtree(device=device, max_blocks=100)
    state = SimulationState(Np=Np, dtype=wp.vec4, device=device, max_blocks=100, basis=basis)
    
    # 1. Setup Grid
    qt.uniform_refine(level=1, state=state, basis=basis)
    qt.refine_blocks([1, 3], state, basis) # Right side fine
    
    # 2. Pulse in Fine block (Right side)
    x_np = state.x.numpy()
    y_np = state.y.numpy()
    params = get_params()
    
    # Gaussian pulse centered at x=0.5 (Fine side)
    # Moving Left (u=-1)
    q_np = np.zeros((100, Np, 4), dtype=np.float32)
    rho = 1.0 + np.exp(-100.0 * ((x_np - 0.5)**2 + y_np**2))
    u = -1.0
    v = 0.0
    p = 1.0
    
    q_np[:, :, 0] = rho
    q_np[:, :, 1] = rho * u
    q_np[:, :, 2] = rho * v
    q_np[:, :, 3] = p / (params.gamma - 1.0) + 0.5 * rho * (u**2 + v**2)
    
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    def get_mass():
        active = state.active_block_indices.numpy()[:qt.num_blocks]
        q = state.q.numpy()
        levels = qt.block_levels.numpy()
        w2d = np.kron(basis.weights_1d.numpy(), basis.weights_1d.numpy())
        total_mass = 0.0
        L = 2.0
        for idx in active:
            J = (L / (2.0**(levels[idx] + 1)))**2
            total_mass += np.sum(q[idx, :, 0] * w2d) * J
        return total_mass

    initial_mass = get_mass()
    
    # 3. Take a few steps (Euler Step for simplicity)
    dt = 1e-5
    for step in range(10):
        state.rhs.zero_()
        
        # Volume
        wp.launch(
            kernel=compute_fr_update,
            dim=qt.num_blocks * Np,
            inputs=[state.q, state.active_block_indices, state.neighbors, qt.num_blocks, state.rhs,
                    basis.nodes_1d, basis.D1D, basis.dg_L, basis.dg_R, qt.root_bounds_wp, qt.block_levels, params, 0.0],
            device=device
        )
        
        # Mortars
        num_mortars = int(qt.num_mortars.numpy()[0])
        if num_mortars > 0:
            wp.launch(
                kernel=compute_mortar_fluxes,
                dim=num_mortars,
                inputs=[state.q, state.rhs, qt.mortar_list, qt.num_mortars, basis.nodes_1d, basis.face_nodes,
                        basis.dg_L, basis.dg_R, basis.P_left, basis.P_right, basis.R_left, basis.R_right,
                        qt.root_bounds_wp, qt.block_levels, params],
                device=device
            )
            
        # q = q + dt * rhs
        q_host = state.q.numpy()
        rhs_host = state.rhs.numpy()
        q_host += dt * rhs_host
        state.q = wp.array(q_host, dtype=wp.vec4, device=device)
        
    final_mass = get_mass()
    
    # Verify conservation
    assert pytest.approx(final_mass, abs=1e-3) == initial_mass
    print("\nTest 2: Mortar Conservation Passed")

