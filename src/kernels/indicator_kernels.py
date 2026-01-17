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
    
    Se = log10( (q - F*q)^2 / q^2 )  <-- Actually usually just the ratio, log is done in threshold check or here.
    The instruction says: Se = E_diff / E.
    
    Args:
        q: State array (max_blocks, Np, 4)
        active_indices: Array of active block indices
        filter_matrix: (Np, Np) spectral filter matrix
        element_indicator: Output array (max_blocks)
        num_active: Number of active blocks
        component: The state component to analyze (e.g., 0 for Density)
    """
    tid = wp.tid()
    if tid >= num_active:
        return
        
    block_idx = active_indices[tid]
    
    # Get Np from filter matrix dimensions
    Np = filter_matrix.shape[0]
    
    # Calculate Energies
    E_total = float(0.0)
    E_diff = float(0.0)
    
    # Loop over nodes i
    for i in range(Np):
        # Get value at node i
        val_i = q[block_idx, i][component]
        
        # Apply Filter: tilde_val_i = sum_j (F[i, j] * val_j)
        tilde_val_i = float(0.0)
        for j in range(Np):
            val_j = q[block_idx, j][component]
            # filter_matrix is (Np, Np)
            # Assuming row-major storage for matrix multiplication F * q
            f_ij = filter_matrix[i, j] 
            tilde_val_i += f_ij * val_j
            
        # Accumulate Energies
        diff = val_i - tilde_val_i
        
        E_total += val_i * val_i
        E_diff += diff * diff
        
    # Avoid division by zero
    Se = 0.0
    if E_total > 1e-20:
        Se = E_diff / E_total
        
    element_indicator[block_idx] = Se

@wp.kernel
def mark_troubled_cells(
    element_indicator: Any,
    active_indices: Any,
    solver_mode: Any,
    threshold: float,
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
