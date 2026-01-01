import pytest
import warp as wp
import numpy as np
from src.core.simulation_state import SimulationState

def test_simulation_state_allocation(device):
    """
    Verifies that SimulationState allocates buffers of correct size and type on the specified device.
    """
    Np = 16
    max_blocks = 100
    
    # vec4 state (Euler), float32 scalars
    state = SimulationState(
        Np=Np, 
        dtype=wp.vec4, 
        device=device, 
        scalar_dtype=wp.float32, 
        max_blocks=max_blocks,
        use_filtering=True
    )
    
    # 1. Check Primary State Buffers
    # Expect (MAX_BLOCKS, Np) of vec4
    assert state.q.shape == (max_blocks, Np)
    assert state.q.dtype == wp.vec4
    
    # 2. Check Scalar Buffers (Grid)
    # Expect (MAX_BLOCKS, Np) of float32
    assert state.x.shape == (max_blocks, Np)
    assert state.x.dtype == wp.float32
    assert state.y.shape == (max_blocks, Np)
    
    # 3. Check Block Management
    assert state.active_block_indices.shape == (max_blocks,)
    assert state.active_block_indices.dtype == wp.int32
    
    # 4. Check Optional Buffers
    assert state.filter_buffer is not None
    assert state.filter_buffer.shape == (max_blocks, Np)

@pytest.mark.parametrize("device", ["cpu"])
def test_simulation_state_zeroing(device):
    """Test utility methods like zero_rhs."""
    wp.init()
    Np = 4
    state = SimulationState(Np=Np, dtype=wp.float32, device=device, max_blocks=10)
    
    # Manually dirty the RHS
    rhs_np = np.ones((10, Np), dtype=np.float32)
    state.rhs = wp.array(rhs_np, dtype=wp.float32, device=device)
    
    state.zero_rhs()
    
    res = state.rhs.numpy()
    assert np.all(res == 0.0)
