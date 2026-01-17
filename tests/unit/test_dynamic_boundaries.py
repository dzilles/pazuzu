import warp as wp
import numpy as np
import pytest
from src.kernels import boundary_conditions as bc
from src.kernels import grid_kernels
from src.kernels import fr_kernels
from src.kernels import utils as u
from src.kernels.structs import BoundaryState32, EquationParams32

@pytest.fixture(scope="module")
def device():
    return "cuda" if wp.get_cuda_device_count() > 0 else "cpu"

def test_boundary_tagging_kernel(device):
    """
    Phase 1 Verification:
    Checks if 'tag_quadtree_boundaries' correctly populates bc_mask 
    based on the block's position in the domain.
    """
    wp.init()
    
    # --- Setup Domain ---
    # Create 4 Blocks (Level 1 Refinement)
    # Block 0: Bottom-Left (Morton 4)
    # Block 1: Bottom-Right (Morton 5)
    # Block 2: Top-Left (Morton 6)
    # Block 3: Top-Right (Morton 7)
    num_blocks = 4
    active_indices = wp.array(np.arange(num_blocks, dtype=np.int32), dtype=wp.int32, device=device)
    block_levels = wp.full(num_blocks, 1, dtype=wp.int32, device=device)
    
    # Morton Encoding for Level 1 (sentinel bit at 1<<2)
    # (0,0) -> 100|00 = 4
    # (1,0) -> 100|01 = 5
    # (0,1) -> 100|10 = 6
    # (1,1) -> 100|11 = 7
    block_morton_codes = wp.array(np.array([4, 5, 6, 7], dtype=np.int32), dtype=wp.int32, device=device)
    
    # --- BC Configuration ---
    # Indices in our imaginary bc_data array
    IDX_LEFT = 10
    IDX_RIGHT = 11
    IDX_BOTTOM = 12
    IDX_TOP = 13
    
    # Output Mask
    bc_mask = wp.full((num_blocks, 4), -1, dtype=wp.int32, device=device)
    
    # --- Run Kernel ---
    wp.launch(
        kernel=grid_kernels.tag_quadtree_boundaries,
        dim=num_blocks,
        inputs=[
            bc_mask,
            active_indices,
            num_blocks,
            block_morton_codes,
            block_levels,
            IDX_LEFT, IDX_RIGHT, IDX_BOTTOM, IDX_TOP
        ],
        device=device
    )
    
    # --- Verify Results ---
    mask = bc_mask.numpy()
    
    # Block 0 (Bottom-Left): Should have Left(3) and Bottom(0) set
    assert mask[0, 3] == IDX_LEFT, f"Block 0 Left face incorrect, got {mask[0, 3]}"
    assert mask[0, 0] == IDX_BOTTOM, f"Block 0 Bottom face incorrect, got {mask[0, 0]}"
    assert mask[0, 1] == -1, "Block 0 Right face should be internal"
    assert mask[0, 2] == -1, "Block 0 Top face should be internal"

    # Block 1 (Bottom-Right): Should have Right(1) and Bottom(0) set
    assert mask[1, 1] == IDX_RIGHT, "Block 1 Right face incorrect"
    assert mask[1, 0] == IDX_BOTTOM, "Block 1 Bottom face incorrect"
    
    # Block 3 (Top-Right): Should have Right(1) and Top(2) set
    assert mask[3, 1] == IDX_RIGHT, "Block 3 Right face incorrect"
    assert mask[3, 2] == IDX_TOP, "Block 3 Top face incorrect"

    print("\n[Pass] Boundary Tagging Kernel Correctly Identified Faces.")


def test_fr_dynamic_dispatch(device):
    """
    Phase 2 Verification:
    Checks if 'compute_fr_update' actually applies different physics
    when we switch the BC ID in bc_mask.
    """
    wp.init() 
    
    # --- Setup Single Block ---
    # We will look at the LEFT face of a single block.
    # Neighbors = -1 (Boundary)
    num_blocks = 1
    Np = 4 # 2x2 nodes
    active_indices = wp.array(np.array([0], dtype=np.int32), dtype=wp.int32, device=device)
    neighbors = wp.full((1, 4), -1, dtype=wp.int32, device=device)
    
    # State: Uniform Flow (rho=1, u=100, v=0, p=1e5)
    # This flow is going RIGHT.
    rho, u_vel, v_vel, p = 1.0, 100.0, 0.0, 100000.0
    gamma = 1.4
    E = p/(gamma-1) + 0.5*rho*(u_vel**2 + v_vel**2)
    q_init = np.array([rho, rho*u_vel, rho*v_vel, E], dtype=np.float32)
    q = wp.array(np.tile(q_init, (1, Np, 1)), dtype=wp.vec4, device=device)
    
    # --- Setup BC Data ---
    # Entry 0: INLET (Matches internal flow) -> Expect minimal flux disturbance
    # Entry 1: WALL (Slip Wall) -> Velocity becomes 0 normal to face -> High flux change
    
    bc_structs = np.zeros(2, dtype=BoundaryState32.numpy_dtype())
    
    # Setup Inlet (Index 0)
    bc_structs[0]['type'] = bc.BC_INLET
    bc_structs[0]['v0'] = rho
    bc_structs[0]['v1'] = u_vel
    bc_structs[0]['v2'] = v_vel
    bc_structs[0]['v3'] = p
    
    # Setup Wall (Index 1)
    bc_structs[1]['type'] = bc.BC_WALL # Slip Wall
    
    bc_data = wp.array(bc_structs, dtype=BoundaryState32, device=device)
    
    # --- Helper Params ---
    params = EquationParams32()
    params.gamma = gamma
    params.rho_inf = rho
    params.u_inf = u_vel
    params.v_inf = v_vel
    params.p_inf = p
    params.flux_type = 0 # Rusanov
    params.rho_floor = 1e-8
    params.p_floor = 1e-8
    params.half = 0.5
    params.one = 1.0
    
    # --- Run 1: Apply INLET to Left Face ---
    bc_mask_inlet = wp.full((1, 4), -1, dtype=wp.int32, device=device)
    # Set Left Face (Index 3) to Inlet (Index 0)
    wp.copy(bc_mask_inlet, wp.array(np.array([[-1, -1, -1, 0]], dtype=np.int32), dtype=wp.int32, device=device))
    
    rhs_inlet = wp.zeros((1, Np), dtype=wp.vec4, device=device)
    
    # Mock Geometry Arrays (Needed for kernel args)
    nodes_1d = wp.array(np.array([-1.0, 1.0], dtype=np.float32), dtype=wp.float32, device=device)
    D1D = wp.zeros((2,2), dtype=wp.float32, device=device)
    dg_L = wp.array(np.array([1.0, 1.0], dtype=np.float32), dtype=wp.float32, device=device)
    dg_R = wp.array(np.array([1.0, 1.0], dtype=np.float32), dtype=wp.float32, device=device)
    root_bounds = wp.vec4(-1.0, -1.0, 1.0, 1.0)
    block_levels = wp.zeros(1, dtype=wp.int32, device=device)
    solver_mode = wp.zeros(1, dtype=wp.int32, device=device)
    x = wp.zeros((1, Np), dtype=wp.float32, device=device)
    y = wp.zeros((1, Np), dtype=wp.float32, device=device)

    # Launch Kernel (INLET)
    wp.launch(
        kernel=fr_kernels.compute_fr_update,
        dim=1 * Np,
        inputs=[
            q, active_indices, neighbors, solver_mode,
            bc_mask_inlet, bc_data, x, y,
            1, rhs_inlet,
            nodes_1d, D1D, dg_L, dg_R, root_bounds, block_levels,
            params, 0.0
        ],
        device=device
    )
    
    # --- Run 2: Apply WALL to Left Face ---
    bc_mask_wall = wp.full((1, 4), -1, dtype=wp.int32, device=device)
    # Set Left Face (Index 3) to Wall (Index 1)
    wp.copy(bc_mask_wall, wp.array(np.array([[-1, -1, -1, 1]], dtype=np.int32), dtype=wp.int32, device=device))
    
    rhs_wall = wp.zeros((1, Np), dtype=wp.vec4, device=device)
    
    # Launch Kernel (WALL)
    wp.launch(
        kernel=fr_kernels.compute_fr_update,
        dim=1 * Np,
        inputs=[
            q, active_indices, neighbors, solver_mode,
            bc_mask_wall, bc_data, x, y,
            1, rhs_wall,
            nodes_1d, D1D, dg_L, dg_R, root_bounds, block_levels,
            params, 0.0
        ],
        device=device
    )
    
    # --- Comparison ---
    # The RHS should differ significantly. 
    # Inlet allows flow through -> Flux ~ 0 (since q_inner == q_inlet).
    # Wall reflects flow -> Stagnation -> Large pressure forces at boundary.
    
    res_inlet = rhs_inlet.numpy()
    res_wall = rhs_wall.numpy()
    
    diff = np.linalg.norm(res_inlet - res_wall)
    print(f"\nRHS Difference (Inlet vs Wall): {diff}")
    
    assert diff > 1.0, f"Changing BC from Inlet to Wall had no effect! Diff: {diff}. Dynamic dispatch failed."
    print("[Pass] Dynamic Dispatch successfully altered physics.")
