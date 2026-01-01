import warp as wp
from src.kernels.grid_kernels import morton_decode, morton_encode
from typing import Any

# Simple Linear Probing Hash Map constants
EMPTY = -1
MORTAR_FLAG = -2

@wp.func
def hash_int(k: int):
    # MurmurHash3 integer finalizer
    k ^= k >> 16
    k *= 0x85ebca6b
    k ^= k >> 13
    k *= 0xc2b2ae35
    k ^= k >> 16
    return k

@wp.func
def map_insert(keys: Any, values: Any, capacity: int, key: int, value: int):
    """
    Thread-safe insertion using Atomic Compare-And-Swap (CAS).
    """
    slot = hash_int(key) % capacity
    start_slot = slot
    
    while True:
        # Try to acquire the slot atomically.
        # atomic_cas(array, index, compare_value, new_value)
        # Returns the OLD value at that slot.
        old_key = wp.atomic_cas(keys, slot, EMPTY, key)
        
        # Case 1: Slot was EMPTY, and we successfully claimed it.
        if old_key == EMPTY:
            values[slot] = value
            return
            
        # Case 2: Slot already contained OUR key (idempotent update).
        if old_key == key:
            values[slot] = value
            return
            
        # Case 3: Slot contained a DIFFERENT key (Collision).
        # Linear probe to next slot.
        slot = (slot + 1) % capacity
        
        # Safety break if map is completely full
        if slot == start_slot:
            return

@wp.func
def map_lookup(keys: Any, values: Any, capacity: int, key: int):
    slot = hash_int(key) % capacity
    start_slot = slot
    
    while True:
        existing = keys[slot]
        if existing == EMPTY:
            return EMPTY
        
        if existing == key:
            return values[slot]
            
        slot = (slot + 1) % capacity
        if slot == start_slot:
            return EMPTY

@wp.kernel
def init_hash_map(keys: Any, values: Any):
    tid = wp.tid()
    keys[tid] = EMPTY
    values[tid] = EMPTY

@wp.kernel
def populate_hash_map(
    morton_codes: Any,
    block_levels: Any,
    active_indices: Any, 
    num_active: int,
    map_keys: Any,
    map_values: Any,
    capacity: int
):
    tid = wp.tid()
    if tid >= num_active:
        return
        
    pool_idx = active_indices[tid]
    code = morton_codes[pool_idx]
    
    # Morton code is now unique across levels due to sentinel bit
    map_insert(map_keys, map_values, capacity, code, pool_idx)

@wp.kernel
def mark_mortar_neighbors(
    mortar_list: Any,
    num_mortars: Any,
    neighbors: Any
):
    tid = wp.tid()
    count = num_mortars[0]
    if tid >= count:
        return
        
    # mortar_list: [fine_idx, face, coarse_idx, subface]
    fine_face = mortar_list[tid, 1]
    coarse_idx = mortar_list[tid, 2]
    
    # Opposing face
    coarse_face = 0
    if fine_face == 0: coarse_face = 1
    elif fine_face == 1: coarse_face = 0
    elif fine_face == 2: coarse_face = 3
    elif fine_face == 3: coarse_face = 2
    
    # Mark coarse side as mortar too
    neighbors[coarse_idx, coarse_face] = MORTAR_FLAG
@wp.kernel
def compute_neighbors(
    morton_codes: Any,
    block_levels: Any,
    active_indices: Any,
    num_active: int,
    map_keys: Any,
    map_values: Any,
    map_capacity: int,
    out_neighbors: Any, 
    mortar_list: Any,       # (MAX_MORTARS, 4) int32
    num_mortars: Any,       # (1) int32 counter
    max_mortars: int,
    periodic_x: int, 
    periodic_y: int  
):
    tid = wp.tid()
    if tid >= num_active:
        return

    pool_idx = active_indices[tid]
    code = morton_codes[pool_idx]
    level = block_levels[pool_idx]
    
    ix, iy = morton_decode(code, level)
    grid_dim = 1 << level
    
    # Check 4 directions
    # 0: Left (-x), 1: Right (+x), 2: Bottom (-y), 3: Top (+y)
    
    for face in range(4):
        nx = ix
        ny = iy
        
        if face == 0: # Left
            nx = ix - 1
        elif face == 1: # Right
            nx = ix + 1
        elif face == 2: # Bottom
            ny = iy - 1
        elif face == 3: # Top
            ny = iy + 1
            
        # Handle Periodicity / Boundaries at CURRENT Level
        valid_n = 1
        if nx < 0:
            if periodic_x != 0: nx = grid_dim - 1
            else: valid_n = 0
        elif nx >= grid_dim:
            if periodic_x != 0: nx = 0
            else: valid_n = 0
            
        if ny < 0:
            if periodic_y != 0: ny = grid_dim - 1
            else: valid_n = 0
        elif ny >= grid_dim:
            if periodic_y != 0: ny = 0
            else: valid_n = 0
            
        n_idx = -1
        
        # 1. Look for neighbor at SAME level
        if valid_n == 1:
            code_n = morton_encode(nx, ny, level)
            n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
            
        # 2. If not found, look for COARSE neighbor (Parent)
        if n_idx == -1 and level > 0 and valid_n == 1:
            px = nx >> 1
            py = ny >> 1
            plevel = level - 1
            
            pcode = morton_encode(px, py, plevel)
            n_idx = map_lookup(map_keys, map_values, map_capacity, pcode)
            
            if n_idx != -1:
                # Found a Coarse Neighbor! This is a MORTAR interface.
                # Mark neighbor as MORTAR_FLAG
                
                # We record this interface in the mortar list.
                # Only need to record it once per face per fine block.
                
                m_idx = wp.atomic_add(num_mortars, 0, 1)
                if m_idx < max_mortars:
                    # Subface calculation
                    # If Face is Left/Right (0/1), subface depends on y (0=bottom, 1=top)
                    # If Face is Bottom/Top (2/3), subface depends on x (0=left, 1=right)
                    subface = 0
                    if face < 2:
                        subface = iy & 1
                    else:
                        subface = ix & 1
                        
                    mortar_list[m_idx, 0] = pool_idx
                    mortar_list[m_idx, 1] = face
                    mortar_list[m_idx, 2] = n_idx
                    mortar_list[m_idx, 3] = subface
                
                # Update neighbor to flag
                n_idx = MORTAR_FLAG
        
        out_neighbors[pool_idx, face] = n_idx