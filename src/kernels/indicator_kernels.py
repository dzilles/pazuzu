import warp as wp
from typing import Any

@wp.kernel
def compute_persson_peraire(
    q: Any,
    active_indices: Any,
    filter_matrix: Any,
    element_indicator: Any,
    num_active: int,
    component: int
):
    """
    Computes the Persson-Peraire smoothness indicator (Se) for each active block.
    """
    tid = wp.tid()
    if tid >= num_active:
        return
        
    block_idx = active_indices[tid]
    
    # Get Np from filter matrix dimensions
    Np = filter_matrix.shape[0]
    
    # Calculate Energies
    # Use subtraction to get a zero of the correct type
    val_0 = q[block_idx, 0][component]
    zero = val_0 - val_0
    
    E_total = zero
    E_diff = zero
    
    # Loop over nodes i
    for i in range(Np):
        # Get value at node i
        val_i = q[block_idx, i][component]
        
        # Apply Filter: tilde_val_i = sum_j (F[i, j] * val_j)
        tilde_val_i = zero
        for j in range(Np):
            val_j = q[block_idx, j][component]
            f_ij = filter_matrix[i, j] 
            tilde_val_i += f_ij * val_j
            
        # Accumulate Energies
        diff = val_i - tilde_val_i
        
        E_total += val_i * val_i
        E_diff += diff * diff
        
    # Avoid division by zero
    Se = zero
    if E_total > 1e-20:
        Se = E_diff / E_total
        
    element_indicator[block_idx] = Se

@wp.kernel
def mark_troubled_cells(
    element_indicator: Any,
    active_indices: Any,
    solver_mode: Any,
    threshold: Any, # Generic to support float32/float64
    num_active: int
):
    """
    Marks blocks as troubled (FV) or smooth (FR) based on the indicator.
    
    Args:
        element_indicator: Array of Se values (max_blocks)
        active_indices: Array of active block indices
        solver_mode: Output array (max_blocks). 0=FR, 1=FV
        threshold: The value above which a block is marked troubled.
        num_active: Number of active blocks
    """
    tid = wp.tid()
    if tid >= num_active:
        return
        
    block_idx = active_indices[tid]
    
    se = element_indicator[block_idx]
    
    if se > threshold:
        solver_mode[block_idx] = 1 # FV
    else:
        solver_mode[block_idx] = 0 # FR
