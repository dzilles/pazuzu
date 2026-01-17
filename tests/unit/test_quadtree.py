import pytest
import warp as wp
import numpy as np
from src.core.basis import Basis
from src.core.simulation_state import SimulationState
from src.geometry.quadtree import Quadtree, ROOT_BOUNDS

def test_quadtree_uniform_refine(device):
    wp.init()
    
    # 1. Setup
    N = 1
    basis = Basis(polynomial_degree=N, device=device)
    Np = basis.Np
    max_blocks = 64
    
    state = SimulationState(Np=Np, dtype=wp.vec4, scalar_dtype=wp.float32, device=device, max_blocks=max_blocks)
    quadtree = Quadtree(device=device, max_blocks=max_blocks)
    
    # 2. Generate Grid Level 1 (2x2 = 4 blocks)
    level = 1
    quadtree.uniform_refine(level, state, basis)
    
    assert quadtree.num_blocks == 4
    
    # 3. Verify Coordinates
    # Copy to host
    x_host = state.x.numpy()
    y_host = state.y.numpy()
    
    # Check Bounds
    # Root bounds are (-1, -1) to (1, 1)
    # Level 1 means 4 quadrants:
    # BL: [-1, 0] x [-1, 0]
    # BR: [0, 1] x [-1, 0]
    # TL: [-1, 0] x [0, 1]
    # TR: [0, 1] x [0, 1]
    
    # In Morton order (Z-curve):
    # 0: (0,0) -> BL
    # 1: (1,0) -> BR
    # 2: (0,1) -> TL
    # 3: (1,1) -> TR
    
    # Let's check Block 0 (BL)
    # Nodes should be in [-1, 0] range.
    assert np.all(x_host[0, :] >= -1.0) and np.all(x_host[0, :] <= 0.0)
    assert np.all(y_host[0, :] >= -1.0) and np.all(y_host[0, :] <= 0.0)
    
    # Let's check Block 3 (TR)
    assert np.all(x_host[3, :] >= 0.0) and np.all(x_host[3, :] <= 1.0)
    assert np.all(y_host[3, :] >= 0.0) and np.all(y_host[3, :] <= 1.0)
    
    # Check precise node locations for Block 0
    # GLL nodes for N=1 are [-1, 1] in reference.
    # Mapped to [-1, 0]: -1 -> -1, 1 -> 0.
    # So physical x should be -1 and 0.
    expected_x = np.array([-1.0, 0.0, -1.0, 0.0], dtype=np.float32) # Tensor product ordering
    # Ordering depends on Basis implementation. Basis.nodes_2d uses meshgrid(x, y).flatten().
    # nodes_1d = [-1, 1]
    # meshgrid x = [[-1, 1], [-1, 1]], y = [[-1, -1], [1, 1]]
    # flattened: (-1, -1), (1, -1), (-1, 1), (1, 1)
    # So x coords: -1, 1, -1, 1
    # y coords: -1, -1, 1, 1
    
    # Mapped to Block 0 (BL) [-1, 0] x [-1, 0]:
    # x: -1, 0, -1, 0
    # y: -1, -1, 0, 0
    
    # ... existing assertions ...
    assert np.allclose(x_host[0, :], [-1.0, 0.0, -1.0, 0.0], atol=1e-5)
    assert np.allclose(y_host[0, :], [-1.0, -1.0, 0.0, 0.0], atol=1e-5)

def test_quadtree_connectivity(device):
    """Test neighbor finding logic for a uniform grid."""
    wp.init()
    basis = Basis(polynomial_degree=1, device=device)
    max_blocks = 16
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=max_blocks)
    quadtree = Quadtree(device=device, max_blocks=max_blocks)
    
    # 2x2 Grid (Level 1)
    # BL=4, BR=5, TL=6, TR=7 (New Morton codes)
    # (ix, iy): (0,0)=4, (1,0)=5, (0,1)=6, (1,1)=7
    # BL (0,0) -> Neighbors: B=-1, R=5, T=6, L=-1
    # BR (1,0) -> Neighbors: B=-1, R=-1, T=7, L=4
    # TL (0,1) -> Neighbors: B=4, R=7, T=-1, L=-1
    # TR (1,1) -> Neighbors: B=5, R=-1, T=-1, L=6
    
    quadtree.uniform_refine(1, state, basis)
    
    neighbors = state.neighbors.numpy()
    codes = quadtree.block_morton_codes.numpy()
    
    # Identify pool indices
    idx_4 = -1; idx_5 = -1; idx_6 = -1; idx_7 = -1
    for i in range(4):
        if codes[i] == 4: idx_4 = i
        if codes[i] == 5: idx_5 = i
        if codes[i] == 6: idx_6 = i
        if codes[i] == 7: idx_7 = i
    
    # Block 4 Neighbors (0:Bottom, 1:Right, 2:Top, 3:Left)
    assert neighbors[idx_4, 0] == -1    # Bottom
    assert neighbors[idx_4, 1] == idx_5 # Right (BR)
    assert neighbors[idx_4, 2] == idx_6 # Top (TL)
    assert neighbors[idx_4, 3] == -1    # Left
    
    # Block 7 Neighbors (TR)
    assert neighbors[idx_7, 0] == idx_5 # Bottom (BR)
    assert neighbors[idx_7, 1] == -1    # Right
    assert neighbors[idx_7, 2] == -1    # Top
    assert neighbors[idx_7, 3] == idx_6 # Left (TL)

def test_quadtree_connectivity_level2(device):
    """Test connectivity for a Level 2 (4x4 = 16 blocks) grid with internal blocks."""
    wp.init()
    basis = Basis(polynomial_degree=1, device=device)
    max_blocks = 32
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=max_blocks)
    quadtree = Quadtree(device=device, max_blocks=max_blocks)
    
    # 4x4 Grid (Level 2)
    quadtree.uniform_refine(2, state, basis)
    
    neighbors = state.neighbors.numpy()
    codes = quadtree.block_morton_codes.numpy()
    
    # Helper to find pool index by (ix, iy) at level 2
    def get_idx(ix, iy):
        code = Quadtree.morton_encode(ix, iy, 2)
        for i in range(16):
            if codes[i] == code: return i
        return -1

    # Let's check an internal block at (ix=1, iy=1)
    idx_internal = get_idx(1, 1)
    # Neighbors (B, R, T, L):
    # Bottom: (1, 0)
    # Right: (2, 1)
    # Top: (1, 2)
    # Left: (0, 1)
    
    assert neighbors[idx_internal, 0] == get_idx(1, 0)
    assert neighbors[idx_internal, 1] == get_idx(2, 1)
    assert neighbors[idx_internal, 2] == get_idx(1, 2)
    assert neighbors[idx_internal, 3] == get_idx(0, 1)
    
    # Let's check a corner block (3, 3)
    idx_corner = get_idx(3, 3)
    assert neighbors[idx_corner, 0] == get_idx(3, 2) # Bottom
    assert neighbors[idx_corner, 1] == -1            # Right
    assert neighbors[idx_corner, 2] == -1            # Top
    assert neighbors[idx_corner, 3] == get_idx(2, 3) # Left

def test_quadtree_limits(device):
    """Test error handling for Quadtree limits."""
    wp.init()
    basis = Basis(polynomial_degree=1, device=device)
    # Small pool
    quadtree = Quadtree(device=device, max_blocks=3) 
    state = SimulationState(Np=basis.Np, dtype=wp.vec4, device=device, max_blocks=3)
    
    # Level 1 requires 4 blocks. Should fail.
    with pytest.raises(ValueError, match="MAX_BLOCKS"):
        quadtree.uniform_refine(1, state, basis)

def test_morton_encoding_logic(device):
    
    # (0, 0) Level 0 -> Root code = 1
    assert Quadtree.morton_encode(0, 0, 0) == 1
    
    # (0, 0) Level 1 -> (1 << 2) | 0 = 4
    assert Quadtree.morton_encode(0, 0, 1) == 4
    # (1, 0) Level 1 -> (1 << 2) | 1 = 5
    assert Quadtree.morton_encode(1, 0, 1) == 5
    # (0, 1) Level 1 -> (1 << 2) | 2 = 6
    assert Quadtree.morton_encode(0, 1, 1) == 6
    # (1, 1) Level 1 -> (1 << 2) | 3 = 7
    assert Quadtree.morton_encode(1, 1, 1) == 7
    
    # (2, 2) Level 2 -> (1 << 4) | interleaved(2, 2)
    # interleaved(2, 2) = 12
    # result = 16 | 12 = 28
    assert Quadtree.morton_encode(2, 2, 2) == 28

    # Verify uniqueness across levels (The original bug)
    # Level 1 Block 2:
    l1_c2 = Quadtree.morton_encode(0, 1, 1) # code 6
    # Level 2 Block 2 of Parent 0:
    l2_c2 = (Quadtree.morton_encode(0, 0, 1) << 2) | 2 # code (4 << 2) | 2 = 18
    assert l1_c2 != l2_c2
