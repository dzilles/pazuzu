import warp as wp
from typing import Any
from src.kernels import utils as u

@wp.kernel
def rk_stage_1(
    q: Any,
    rhs: Any,
    dt: Any,
    q_out: Any,
    active_indices: Any
):
    """
    Performs the first stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    block_idx, node_idx = wp.tid()  # type: ignore
    pool_idx = active_indices[block_idx] # type: ignore
    
    q_out[pool_idx, node_idx] = q[pool_idx, node_idx] + dt * rhs[pool_idx, node_idx] # type: ignore

@wp.kernel
def rk_stage_2(
    q: Any,      # Q_n (Initial state)
    q_1: Any,    # Q(1) from Stage 1
    rhs: Any,    # RHS(Q(1))
    dt: Any,
    q_out: Any,   # Destination for Q(2)
    c1: Any, # 0.75
    c2: Any, # 0.25
    active_indices: Any
):
    """
    Performs the second stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    block_idx, node_idx = wp.tid()  # type: ignore
    pool_idx = active_indices[block_idx] # type: ignore
    
    q_out[pool_idx, node_idx] = c1 * q[pool_idx, node_idx] + c2 * (q_1[pool_idx, node_idx] + dt * rhs[pool_idx, node_idx]) # type: ignore

@wp.kernel
def rk_stage_3(
    q: Any,      # Q_n
    q_2: Any,    # Q(2) from Stage 2
    rhs: Any,    # RHS(Q(2))
    dt: Any,
    q_out: Any,   # Destination for Q(n+1)
    c1: Any, # 1/3
    c2: Any,  # 2/3
    active_indices: Any
    ):
    """
    Performs the third and final stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    block_idx, node_idx = wp.tid()  # type: ignore
    pool_idx = active_indices[block_idx] # type: ignore
    
    q_out[pool_idx, node_idx] = c1 * q[pool_idx, node_idx] + c2 * (q_2[pool_idx, node_idx] + dt * rhs[pool_idx, node_idx]) # type: ignore

@wp.kernel
def rk4_stage_update(
    q_old: Any,
    rhs: Any,
    q_accum: Any,
    q_next: Any,
    dt: Any,
    weight_accum: Any,
    weight_next: Any,
    active_indices: Any
):
    """
    Updates the accumulator and prepares the next stage state for RK4.
    
    q_accum += weight_accum * dt * rhs
    q_next = q_old + weight_next * dt * rhs
    """
    block_idx, node_idx = wp.tid()  # type: ignore
    pool_idx = active_indices[block_idx] # type: ignore
    
    term = dt * rhs[pool_idx, node_idx] # type: ignore
    q_accum[pool_idx, node_idx] = q_accum[pool_idx, node_idx] + weight_accum * term # type: ignore
    q_next[pool_idx, node_idx] = q_old[pool_idx, node_idx] + weight_next * term # type: ignore

@wp.kernel
def rk4_final_update(
    rhs: Any,
    q_accum: Any,
    dt: Any,
    weight_accum: Any,
    active_indices: Any
):
    """
    Performs the final update to the accumulator for RK4.
    
    q_accum += weight_accum * dt * rhs
    """
    block_idx, node_idx = wp.tid()  # type: ignore
    pool_idx = active_indices[block_idx] # type: ignore
    
    q_accum[pool_idx, node_idx] = q_accum[pool_idx, node_idx] + weight_accum * dt * rhs[pool_idx, node_idx] # type: ignore

@wp.kernel
def apply_filter_matrix(
    q: Any,
    F: Any,
    q_out: Any
):
    block_idx, node_idx = wp.tid()  # type: ignore
    Np = q.shape[1]
    
    template = q[0, 0][0]
    zero = template - template
    res = u.make_vec4_generic(zero, zero, zero, zero)
    
    for j in range(Np):
        res += q[block_idx, j] * F[node_idx, j] # type: ignore
        
    q_out[block_idx, node_idx] = res # type: ignore

@wp.kernel
def apply_filter_matrix_indirect(
    q: Any,
    active_indices: Any,
    F: Any,
    q_out: Any
):
    block_idx, node_idx = wp.tid()  # type: ignore
    pool_idx = active_indices[block_idx] # type: ignore
    Np = q.shape[1]
    
    template = q[0, 0][0]
    zero = template - template
    res = u.make_vec4_generic(zero, zero, zero, zero)
    
    for j in range(Np):
        res += q[pool_idx, j] * F[node_idx, j] # type: ignore
        
    q_out[pool_idx, node_idx] = res # type: ignore

@wp.kernel
def check_nan(
    q: Any,
    has_nan: Any
):
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    val = q[e, i]  # type: ignore # Warp type inference
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
    
    block_idx = tid // Np # type: ignore
    node_idx = tid % Np # type: ignore
    
    pool_idx = active_indices[block_idx] # type: ignore
    
    val = q[pool_idx, node_idx]  # type: ignore # Warp type inference
    
    for c in range(4):
        if wp.isnan(val[c]) or wp.isinf(val[c]):
            wp.atomic_add(has_nan, 0, 1)

@wp.kernel
def check_nan_vec2(
    q: Any,
    has_nan: Any
):
    e, i = wp.tid()  # type: ignore # Warp returns a tuple at runtime
    val = q[e, i]  # type: ignore # Warp type inference
    if wp.isnan(val[0]) or wp.isnan(val[1]) or wp.isinf(val[0]) or wp.isinf(val[1]):
        wp.atomic_add(has_nan, 0, 1)
