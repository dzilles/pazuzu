import warp as wp
from typing import Any

from src.kernels import utils as u
from src.kernels.grid_kernels import morton_decode, morton_encode
from src.kernels.connectivity_kernels import map_lookup

@wp.kernel
def balance_refine_flags(
    active_indices: Any,
    num_active: int,
    morton_codes: Any,
    block_levels: Any,
    map_keys: Any,
    map_values: Any,
    map_capacity: int,
    periodic_x: int,
    periodic_y: int,
    refine_flags: Any,  # In/Out
    changed: Any        # Out (size 1)
):
    tid = wp.tid() # type: ignore # type: ignore
    if tid >= num_active:
        return
        
    pool_idx = active_indices[tid]
    if refine_flags[pool_idx] != 1:
        return
        
    level = block_levels[pool_idx]
    if level == 0:
        return
        
    code = morton_codes[pool_idx]
    ix, iy = morton_decode(code, level)
    grid_dim = 1 << level
    
    # Check 4 directions: Bottom, Right, Top, Left
    for face in range(4):
        nx = ix
        ny = iy
        
        if face == 0:
            ny = iy - 1
        elif face == 1:
            nx = ix + 1
        elif face == 2:
            ny = iy + 1
        elif face == 3:
            nx = ix - 1
        
        # Periodic Wrap / Boundary Check
        if nx < 0:
            if periodic_x != 0:
                nx = grid_dim - 1
            else:
                continue
        elif nx >= grid_dim:
            if periodic_x != 0:
                nx = 0
            else:
                continue
            
        if ny < 0:
            if periodic_y != 0:
                ny = grid_dim - 1
            else:
                continue
        elif ny >= grid_dim:
            if periodic_y != 0:
                ny = 0
            else:
                continue
            
        # If we refine level L, any level L-1 neighbor MUST also be refined
        # to maintain 2:1 balance (otherwise we'd have L+1 next to L-1).
        px = nx >> 1
        py = ny >> 1
        plevel = level - 1
        pcode = morton_encode(px, py, plevel)
        
        n_idx_coarse = map_lookup(map_keys, map_values, map_capacity, pcode)
        
        if n_idx_coarse != -1:
            # Found a neighbor at level L-1. 
            # We must mark it for refinement.
            # Refinement takes priority over coarsening/keeping.
            old_val = wp.atomic_cas(refine_flags, n_idx_coarse, 0, 1)
            if old_val == 0:
                changed[0] = 1
            elif old_val == -1:
                # If it was marked for coarsening, we override it to 1
                # Assignment is safe here as all concurrent writes will be '1'
                refine_flags[n_idx_coarse] = 1
                changed[0] = 1

@wp.kernel
def mark_blocks_gradient(
    q: Any,                 # (MAX_BLOCKS, Np) vec4
    active_indices: Any,    # (num_active) int32
    num_active: int,
    block_levels: Any,      # (MAX_BLOCKS) int32
    refine_flags: Any,      # (MAX_BLOCKS) int32, Output
    refine_threshold: float,
    coarsen_threshold: float,
    max_depth: int
):
    tid = wp.tid() # type: ignore
    if tid >= num_active:
        return
        
    pool_idx = active_indices[tid]
    level = block_levels[pool_idx]
    
    # Compute min/max of density (component 0)
    Np = q.shape[1]
    
    # Use utility to get correct precision
    rho_min = u.get_any_generic(q[0, 0][0], 1.0e20)
    rho_max = u.get_any_generic(q[0, 0][0], -1.0e20)
    
    for i in range(Np):
        rho = q[pool_idx, i][0]
        if rho < rho_min:
            rho_min = rho
        if rho > rho_max:
            rho_max = rho
        
    diff = rho_max - rho_min
    
    if diff > refine_threshold and level < max_depth:
        refine_flags[pool_idx] = 1
    elif diff < coarsen_threshold:
        refine_flags[pool_idx] = -1
    else:
        refine_flags[pool_idx] = 0

@wp.kernel
def mark_blocks_on_interface(
    phi: Any,                 # (MAX_BLOCKS, Np) scalar
    active_indices: Any,      # (num_active) int32
    num_active: int,
    refine_flags: Any         # (MAX_BLOCKS) int32 (In/Out)
):
    """
    Forces refinement if the Immersed Boundary interface (phi=0) passes through the block.
    """
    tid = wp.tid() # type: ignore
    if tid >= num_active:
        return
        
    pool_idx = active_indices[tid]
    Np = phi.shape[1]
    
    # 1. Compute Min/Max of Signed Distance in this block
    # Initialize with the first node's value
    val_0 = phi[pool_idx, 0]
    min_phi = val_0
    max_phi = val_0
    
    for i in range(1, Np):
        val = phi[pool_idx, i]
        if val < min_phi:
            min_phi = val
        if val > max_phi:
            max_phi = val
            
    # 2. Check for Interface crossing
    # If min is negative (solid) and max is positive (fluid), the interface is inside.
    # Also refine if exactly 0.0 (on surface).
    if min_phi * max_phi <= 0.0:
        # Force refinement (Override any 'coarsen' or 'keep' decision)
        refine_flags[pool_idx] = 1

@wp.kernel
def zero_blocks(
    q: Any,                 # (MAX, Np) vec4
    block_indices: Any,     # (num_blocks) int32
    num_blocks: int
):
    tid_block, tid_node = wp.tid() # type: ignore
    if tid_block >= num_blocks: # type: ignore
        return
        
    pool_idx = block_indices[tid_block] # type: ignore
    template = q[0, 0][0]
    zero = template - template
    q[pool_idx, tid_node] = u.make_vec4_generic(zero, zero, zero, zero) # type: ignore

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
    # This function computes ONE node's value (generic precision)
    N1 = P_left.shape[0] # N+1
    row = tid_node // N1
    col = tid_node % N1
    
    use_right_row = (child_quadrant >= 2)
    use_right_col = (child_quadrant % 2 != 0)

    # Accumulate result vector
    template = parent_q[0, 0][0]
    zero = template - template
    res = u.make_vec4_generic(zero, zero, zero, zero)
    
    for k in range(N1):
        p_row_val = P_left[row, k]
        if use_right_row:
            p_row_val = P_right[row, k]
            
        inner_sum = u.make_vec4_generic(zero, zero, zero, zero)
        for ll in range(N1):
            p_col_val = P_left[col, ll]
            if use_right_col:
                p_col_val = P_right[col, ll]
            
            # parent_q index
            q_idx = k * N1 + ll
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
    tid = wp.tid() # type: ignore
    Np = q.shape[1]
    
    op_idx = tid // Np # type: ignore
    node_idx = tid % Np # type: ignore
    
    if op_idx >= num_ops:
        return
    
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

    template = child_q[0, 0][0]
    zero = template - template
    res = u.make_vec4_generic(zero, zero, zero, zero)
    
    for k in range(N1):
        r_row_val = R_left[row, k]
        if use_right_row:
            r_row_val = R_right[row, k]
            
        inner_sum = u.make_vec4_generic(zero, zero, zero, zero)
        for ll in range(N1):
            r_col_val = R_left[col, ll]
            if use_right_col:
                r_col_val = R_right[col, ll]
            
            q_idx = k * N1 + ll
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
    tid = wp.tid() # type: ignore
    Np = q.shape[1]
    
    op_idx = tid // Np # type: ignore
    node_idx = tid % Np # type: ignore
    
    if op_idx >= num_ops:
        return
    
    c_idx = op_list[op_idx, 0]
    p_idx = op_list[op_idx, 1]
    quad = op_list[op_idx, 2]
    
    restrict_block_func(q, c_idx, q, p_idx, quad, node_idx, R_left, R_right)