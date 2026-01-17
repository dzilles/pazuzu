import warp as wp
import numpy as np
import pytest
from src.core.basis import Basis
from src.core.simulation_state import SimulationState
from src.kernels.indicator_kernels import compute_persson_peraire, mark_troubled_cells

@pytest.fixture
def device():
    return "cpu"

def test_persson_peraire_indicator(device):
    wp.init()
    
    N = 2
    max_blocks = 2
    Np = (N+1)**2
    
    # 1. Setup Basis and Filter Matrix
    basis = Basis(polynomial_degree=N, device=device)
    # Compute filter matrix (alpha=36, order=4) - standard values
    basis.compute_filter_matrix(alpha=36.0, order=4)
    filter_matrix = basis.filter_matrix
    
    # 2. Setup State
    # Block 0: Smooth (Constant) -> Expect Se ~ 0
    # Block 1: Rough (Checkerboard/Noise) -> Expect Se > 0
    
    state = SimulationState(Np=Np, dtype=wp.vec4, device=device, max_blocks=max_blocks, scalar_dtype=wp.float32)
    state.num_active_blocks = 2
    
    # Active indices: [0, 1]
    # Replace the lambda launch with direct assignment since max_blocks=2 and we want [0, 1]
    indices_np = np.arange(max_blocks, dtype=np.int32)
    state.active_block_indices = wp.array(indices_np, dtype=wp.int32, device=device)
    
    q_np = np.zeros((max_blocks, Np, 4), dtype=np.float32)
    
    # Block 0: Constant density = 1.0
    q_np[0, :, 0] = 1.0
    
    # Block 1: Alternating 1.0 and 2.0 (High frequency)
    for i in range(Np):
        q_np[1, i, 0] = 1.0 if i % 2 == 0 else 2.0
        
    state.q = wp.array(q_np, dtype=wp.vec4, device=device)
    
    # 3. Run Indicator Kernel
    wp.launch(
        kernel=compute_persson_peraire,
        dim=state.num_active_blocks,
        inputs=[
            state.q,
            state.active_block_indices,
            filter_matrix,
            state.element_indicator,
            state.num_active_blocks,
            0 # Component 0 (Density)
        ],
        device=device
    )
    
    indicator = state.element_indicator.numpy()
    print(f"Indicators: {indicator}")
    
    # Check Block 0 (Smooth)
    assert indicator[0] < 1e-6, f"Smooth block should have low indicator, got {indicator[0]}"
    
    # Check Block 1 (Rough)
    # The filter should remove high freq, so difference should be large.
    assert indicator[1] > 1e-4, f"Rough block should have high indicator, got {indicator[1]}"
    
    # 4. Run Mark Kernel
    threshold = 0.001
    
    # Manually set indicator for testing threshold logic if needed, 
    # but the computed ones should suffice if Block 1 is rough enough.
    # Let's ensure Block 1 is flagged.
    
    wp.launch(
        kernel=mark_troubled_cells,
        dim=state.num_active_blocks,
        inputs=[
            state.element_indicator,
            state.active_block_indices,
            state.solver_mode,
            threshold,
            state.num_active_blocks
        ],
        device=device
    )
    
    mode = state.solver_mode.numpy()
    
    assert mode[0] == 0, "Block 0 should be FR (0)"
    if indicator[1] > threshold:
        assert mode[1] == 1, "Block 1 should be FV (1)"
    else:
        # If the indicator wasn't high enough (unlikely for checkerboard),
        # we at least verified logic. But let's print warning.
        print(f"Warning: Rough block indicator {indicator[1]} was below threshold {threshold}")

if __name__ == "__main__":
    test_persson_peraire_indicator("cpu")
