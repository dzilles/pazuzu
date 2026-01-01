import warp as wp
from src.kernels.grid_kernels import morton_decode, morton_encode
from typing import Any

# Simple Linear Probing Hash Map constants
EMPTY = -1

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
    level = block_levels[pool_idx]
    
    # Key includes level to distinguish overlapping codes at different depths
    # (code << 4) | level. Assumes level < 16.
    key = (code << 4) | level
    
    map_insert(map_keys, map_values, capacity, key, pool_idx)

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
    periodic_x: int, 
    periodic_y: int  
):
    tid = wp.tid()
    if tid >= num_active:
        return

    pool_idx = active_indices[tid]
    code = morton_codes[pool_idx]
    level = block_levels[pool_idx]
    
    ix, iy = morton_decode(code)
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
            code_n = morton_encode(nx, ny)
            key_n = (code_n << 4) | level
            n_idx = map_lookup(map_keys, map_values, map_capacity, key_n)
            
        # 2. If not found, look for COARSE neighbor (Parent)
        # Note: We do NOT look for Finer neighbors here. 
        # If the neighbor is finer, we point to -1 (Ghost handling will fix this later via the finer neighbor pointing to us).
        if n_idx == -1 and level > 0:
            # We must re-evaluate boundary conditions for the parent level?
            # Actually, standard approach: Project coord to parent space.
            # Parent coord: px = nx >> 1, py = ny >> 1.
            # But wait, we need the neighbor of our PARENT.
            # Let's think: "My left neighbor" might be a large block.
            # Its coordinate in the coarse grid is (ix >> 1) - 1 ??
            # No, if I am on the left edge of my parent, my neighbor is in a different parent.
            # If I am on the right edge of my parent, my neighbor is my sibling.
            
            # Simple approach: logic above computed (nx, ny) at fine level.
            # If (nx, ny) falls into a coarse block, that coarse block covers (nx>>1, ny>>1).
            # So, we check the parent key corresponding to (nx>>1, ny>>1).
            
            # Re-check boundary for parent level?
            # If we wrapped around at fine level, (nx, ny) are valid fine coords.
            # So (nx>>1, ny>>1) are valid coarse coords.
            # If we hit valid_n=0 (physical boundary), we stop anyway.
            
            if valid_n == 1:
                px = nx >> 1
                py = ny >> 1
                plevel = level - 1
                
                pcode = morton_encode(px, py)
                pkey = (pcode << 4) | plevel
                n_idx = map_lookup(map_keys, map_values, map_capacity, pkey)
        
        out_neighbors[pool_idx, face] = n_idx