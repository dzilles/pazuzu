import warp as wp
from typing import Any

@wp.kernel
def mark_blocks_gradient(
    q: Any,                 # (MAX_BLOCKS, Np) vec4
    active_indices: Any,    # (num_active) int32
    num_active: int,
    refine_flags: Any,      # (MAX_BLOCKS) int32, Output
    threshold: float
):
    tid = wp.tid()
    if tid >= num_active:
        return
        
    pool_idx = active_indices[tid]
    
    # Compute min/max of density (component 0)
    Np = q.shape[1]
    
    rho_min = float(1.0e20)
    rho_max = float(-1.0e20)
    
    for i in range(Np):
        rho = q[pool_idx, i][0]
        if rho < rho_min: rho_min = rho
        if rho > rho_max: rho_max = rho
        
    diff = rho_max - rho_min
    
    if diff > threshold:
        refine_flags[pool_idx] = 1
    else:
        refine_flags[pool_idx] = 0

@wp.func
def prolongate_block_func(
    parent_q: Any,     # (MAX, Np) vec4
    p_idx: int,
    child_q: Any,      # (MAX, Np) vec4
    c_idx: int,
    child_quadrant: int, # 0..3
    tid_node: int,     # 0..Np-1
    P_left: Any,       # (N1, N1)
    P_right: Any       # (N1, N1)
):
    # This function computes ONE node's value (as vec4)
    N1 = P_left.shape[0] # N+1
    row = tid_node // N1
    col = tid_node % N1
    
    use_right_row = (child_quadrant >= 2)
    use_right_col = (child_quadrant % 2 != 0)

    # Accumulate result vector
    # We assume vec4 (float32).
    res = wp.vec4(0.0, 0.0, 0.0, 0.0)
    
    # We can iterate components or do vector math.
    # Vector math is better.
    
    for k in range(N1):
        p_row_val = P_left[row, k]
        if use_right_row: p_row_val = P_right[row, k]
            
        inner_sum = wp.vec4(0.0, 0.0, 0.0, 0.0)
        for l in range(N1):
            p_col_val = P_left[col, l]
            if use_right_col: p_col_val = P_right[col, l]
            
            # parent_q index
            q_idx = k * N1 + l
            inner_sum += parent_q[p_idx, q_idx] * p_col_val
            
        res += inner_sum * p_row_val
        
    child_q[c_idx, tid_node] = res

@wp.kernel
def prolongate_batch(
    q: Any,                 # (MAX, Np) vec4
    op_list: Any,           # (num_ops, 3) -> [parent_idx, child_idx, child_quadrant]
    num_ops: int,
    P_left: Any,
    P_right: Any
):
    tid = wp.tid()
    Np = q.shape[1]
    
    op_idx = tid // Np
    node_idx = tid % Np
    
    if op_idx >= num_ops: return
    
    p_idx = op_list[op_idx, 0]
    c_idx = op_list[op_idx, 1]
    quad = op_list[op_idx, 2]
    
    prolongate_block_func(q, p_idx, q, c_idx, quad, node_idx, P_left, P_right)


@wp.func
def restrict_block_func(
    child_q: Any,      # (MAX, Np) vec4
    c_idx: int,
    parent_q: Any,     # (MAX, Np) vec4
    p_idx: int,
    child_quadrant: int,
    tid_node: int,
    R_left: Any,
    R_right: Any
):
    N1 = R_left.shape[0]
    row = tid_node // N1
    col = tid_node % N1
    
    use_right_row = (child_quadrant >= 2)
    use_right_col = (child_quadrant % 2 != 0)

    res = wp.vec4(0.0, 0.0, 0.0, 0.0)
    
    for k in range(N1):
        r_row_val = R_left[row, k]
        if use_right_row: r_row_val = R_right[row, k]
            
        inner_sum = wp.vec4(0.0, 0.0, 0.0, 0.0)
        for l in range(N1):
            r_col_val = R_left[col, l]
            if use_right_col: r_col_val = R_right[col, l]
            
            q_idx = k * N1 + l
            inner_sum += child_q[c_idx, q_idx] * r_col_val
            
        res += inner_sum * r_row_val
            
    # Atomic Add to parent
    wp.atomic_add(parent_q, p_idx, tid_node, res)

@wp.kernel
def restrict_batch(
    q: Any,                 # (MAX, Np) vec4
    op_list: Any,           # (num_ops, 3) -> [child_idx, parent_idx, child_quadrant]
    num_ops: int,
    R_left: Any,
    R_right: Any
):
    tid = wp.tid()
    Np = q.shape[1]
    
    op_idx = tid // Np
    node_idx = tid % Np
    
    if op_idx >= num_ops: return
    
    c_idx = op_list[op_idx, 0]
    p_idx = op_list[op_idx, 1]
    quad = op_list[op_idx, 2]
    
    restrict_block_func(q, c_idx, q, p_idx, quad, node_idx, R_left, R_right)

