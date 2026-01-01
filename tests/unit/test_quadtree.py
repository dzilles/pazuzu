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
    # BL=0, BR=1, TL=2, TR=3 (Morton)
    # BL (0,0) -> Neighbors: L=-1, R=1, B=-1, T=2
    # BR (1,0) -> Neighbors: L=0, R=-1, B=-1, T=3
    # TL (0,1) -> Neighbors: L=-1, R=3, B=0, T=-1
    # TR (1,1) -> Neighbors: L=2, R=-1, B=1, T=-1
    
    quadtree.uniform_refine(1, state, basis)
    
    neighbors = state.neighbors.numpy()
    
    # Check Block 0 (BL)
    # 0:Left, 1:Right, 2:Bottom, 3:Top
    # Neighbors stores POOL INDICES. Since we filled sequentially:
    # Pool 0 = Morton 0 (BL)
    # Pool 1 = Morton 1 (BR)
    # Pool 2 = Morton 2 (TL)
    # Pool 3 = Morton 3 (TR)
    
    # Block 0 Neighbors
    assert neighbors[0, 0] == -1 # Left
    assert neighbors[0, 1] == 1  # Right (BR)
    assert neighbors[0, 2] == -1 # Bottom
    assert neighbors[0, 3] == 2  # Top (TL)
    
    # Block 3 Neighbors (TR)
    assert neighbors[3, 0] == 2  # Left (TL)
    assert neighbors[3, 1] == -1 # Right
    assert neighbors[3, 2] == 1  # Bottom (BR)
    assert neighbors[3, 3] == -1 # Top

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
    
    # Let's check an internal block at (ix=1, iy=1)
    # Pool index = iy * grid_dim + ix = 1 * 4 + 1 = 5
    # Neighbors (Pool Indices):
    # Left: (0, 1)   -> index 4
    # Right: (2, 1)  -> index 6
    # Bottom: (1, 0) -> index 1
    # Top: (1, 2)    -> index 9
    
    idx_internal = 5
    assert neighbors[idx_internal, 0] == 4 # Left
    assert neighbors[idx_internal, 1] == 6 # Right
    assert neighbors[idx_internal, 2] == 1 # Bottom
    assert neighbors[idx_internal, 3] == 9 # Top
    
    # Let's check a corner block (3, 3) -> index 15
    idx_corner = 15
    assert neighbors[idx_corner, 0] == 14 # Left (2, 3)
    assert neighbors[idx_corner, 1] == -1 # Right
    assert neighbors[idx_corner, 2] == 11 # Bottom (3, 2)
    assert neighbors[idx_corner, 3] == -1 # Top

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
        
    # Default max_depth is 10. Try 11.
    quadtree_big = Quadtree(device=device, max_blocks=10000)
    with pytest.raises(ValueError, match="max_depth"):
        quadtree_big.uniform_refine(11, state, basis)
        
    # Custom max_depth
    quadtree_custom = Quadtree(device=device, max_blocks=10000, max_depth=2)
    with pytest.raises(ValueError, match="max_depth"):
        quadtree_custom.uniform_refine(3, state, basis)

def test_morton_encoding_logic(device):
    from src.geometry.quadtree import morton_encode
    
    # (0, 0) -> 0
    assert morton_encode(0, 0) == 0
    # (1, 0) -> 1
    assert morton_encode(1, 0) == 1
    # (0, 1) -> 2
    assert morton_encode(0, 1) == 2
    # (1, 1) -> 3
    assert morton_encode(1, 1) == 3
    # (2, 2) -> (10, 10)_2. Interleave: 1100_2 = 12
    # x=10, y=10. 
    # part1(x) = ...00100 -> ...00000100
    # part1(y) = ...00100 -> ...00000100
    # y << 1 = ...00001000
    # result = 1100 = 12. Correct.
    assert morton_encode(2, 2) == 12
