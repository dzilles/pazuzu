import warp as wp
from typing import Any

@wp.kernel
def rk_stage_1(
    q: Any,
    rhs: Any,
    dt: Any,
    q_out: Any
):
    """
    Performs the first stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    q_out[e, i] = q[e, i] + dt * rhs[e, i]

@wp.kernel
def rk_stage_2(
    q: Any,      # Q_n (Initial state)
    q_1: Any,    # Q(1) from Stage 1
    rhs: Any,    # RHS(Q(1))
    dt: Any,
    q_out: Any,   # Destination for Q(2)
    c1: Any, # 0.75
    c2: Any  # 0.25
):
    """
    Performs the second stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    # 0.75 * Q_n + 0.25 * (Q_1 + dt * RHS)
    q_out[e, i] = c1 * q[e, i] + c2 * (q_1[e, i] + dt * rhs[e, i])

@wp.kernel
def rk_stage_3(
    q: Any,      # Q_n
    q_2: Any,    # Q(2) from Stage 2
    rhs: Any,    # RHS(Q(2))
    dt: Any,
    q_out: Any,   # Destination for Q(n+1)
    c1: Any, # 1/3
    c2: Any  # 2/3
    ):
    """
    Performs the third and final stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    
    # 1/3 * Q_n + 2/3 * (Q_2 + dt * RHS)
    q_out[e, i] = c1 * q[e, i] + c2 * (q_2[e, i] + dt* rhs[e, i])

@wp.kernel
def apply_filter_matrix(
    q: Any,
    filter_matrix: Any,
    q_out: Any
):
    """
    Applies the spectral filter matrix to the state vector.
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime # i is the node index (row of output)
    
    Np = filter_matrix.shape[1]
    
    # Generic zero initialization
    val = q[e, i] - q[e, i]
    
    for j in range(Np):
        f = filter_matrix[i, j]
        q_val = q[e, j]
        # f is from filter_matrix, which matches dtype.
        val = val + q_val * f
        
    q_out[e, i] = val

@wp.kernel
def rk4_stage_update(
    q_old: Any,
    rhs: Any,
    q_accum: Any,
    q_next: Any,
    dt: Any,
    weight_accum: Any,
    weight_next: Any
):
    """
    Updates the accumulator and prepares the next stage state for RK4.
    
    q_accum += weight_accum * dt * rhs
    q_next = q_old + weight_next * dt * rhs
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    term = dt * rhs[e, i]
    q_accum[e, i] = q_accum[e, i] + weight_accum * term
    q_next[e, i] = q_old[e, i] + weight_next * term

@wp.kernel
def rk4_final_update(
    rhs: Any,
    q_accum: Any,
    dt: Any,
    weight_accum: Any
):
    """
    Performs the final update to the accumulator for RK4.
    
    q_accum += weight_accum * dt * rhs
    """
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    q_accum[e, i] = q_accum[e, i] + weight_accum * dt * rhs[e, i]

@wp.kernel
def check_nan(
    q: Any,
    has_nan: Any
):
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    val = q[e, i]
    # Check each component for NaN or Inf
    for c in range(4):
        if wp.isnan(val[c]) or wp.isinf(val[c]):
            wp.atomic_add(has_nan, 0, 1)

@wp.kernel
def check_nan_indirect(
    q: Any,
    active_indices: Any,
    has_nan: Any
):
    tid = wp.tid()
    Np = q.shape[1]
    
    block_idx = tid // Np
    node_idx = tid % Np
    
    pool_idx = active_indices[block_idx]
    
    val = q[pool_idx, node_idx]
    
    for c in range(4):
        if wp.isnan(val[c]) or wp.isinf(val[c]):
            wp.atomic_add(has_nan, 0, 1)

@wp.kernel
def check_nan_vec2(
    q: Any,
    has_nan: Any
):
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    val = q[e, i]
    if wp.isnan(val[0]) or wp.isnan(val[1]) or wp.isinf(val[0]) or wp.isinf(val[1]):
        wp.atomic_add(has_nan, 0, 1)
