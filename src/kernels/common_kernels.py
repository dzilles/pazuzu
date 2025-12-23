import warp as wp
from typing import Any

@wp.kernel
def rk_stage_1(
    q: wp.array(dtype=Any, ndim=2),
    rhs: wp.array(dtype=Any, ndim=2),
    dt: Any,
    q_out: wp.array(dtype=Any, ndim=2)
):
    """
    Performs the first stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    e, i = wp.tid()
    q_out[e, i] = q[e, i] + dt * rhs[e, i]

@wp.kernel
def rk_stage_2(
    q: wp.array(dtype=Any, ndim=2),      # Q_n (Initial state)
    q_1: wp.array(dtype=Any, ndim=2),    # Q(1) from Stage 1
    rhs: wp.array(dtype=Any, ndim=2),    # RHS(Q(1))
    dt: Any,
    q_out: wp.array(dtype=Any, ndim=2),   # Destination for Q(2)
    c1: Any, # 0.75
    c2: Any  # 0.25
):
    """
    Performs the second stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    e, i = wp.tid()
    
    # 0.75 * Q_n + 0.25 * (Q_1 + dt * RHS)
    q_out[e, i] = c1 * q[e, i] + c2 * (q_1[e, i] + dt * rhs[e, i])

@wp.kernel
def rk_stage_3(
    q: wp.array(dtype=Any, ndim=2),      # Q_n
    q_2: wp.array(dtype=Any, ndim=2),    # Q(2) from Stage 2
    rhs: wp.array(dtype=Any, ndim=2),    # RHS(Q(2))
    dt: Any,
    q_out: wp.array(dtype=Any, ndim=2),   # Destination for Q(n+1)
    c1: Any, # 1/3
    c2: Any  # 2/3
    ):
    """
    Performs the third and final stage of the Low-Storage SSP-RK3 time integration scheme.
    """
    e, i = wp.tid()
    
    # 1/3 * Q_n + 2/3 * (Q_2 + dt * RHS)
    q_out[e, i] = c1 * q[e, i] + c2 * (q_2[e, i] + dt* rhs[e, i])

@wp.kernel
def apply_filter_matrix(
    q: wp.array(dtype=Any, ndim=2),
    filter_matrix: wp.array(dtype=Any, ndim=2),
    q_out: wp.array(dtype=Any, ndim=2)
):
    """
    Applies the spectral filter matrix to the state vector.
    """
    e, i = wp.tid() # i is the node index (row of output)
    
    Np = filter_matrix.shape[1]
    
    # Generic zero initialization
    val = q[e, i] - q[e, i]
    
    for j in range(Np):
        f = filter_matrix[i, j]
        q_val = q[e, j]
        # f is from filter_matrix, which matches dtype.
        val = val + q_val * f
        
    q_out[e, i] = val