import warp as wp
from src.kernels.grid_kernels import morton_decode, morton_encode

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
def map_insert(keys: wp.array(dtype=int), values: wp.array(dtype=int), capacity: int, key: int, value: int):
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
def map_lookup(keys: wp.array(dtype=int), values: wp.array(dtype=int), capacity: int, key: int):
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
def init_hash_map(keys: wp.array(dtype=int), values: wp.array(dtype=int)):
    tid = wp.tid()
    keys[tid] = EMPTY
    values[tid] = EMPTY

@wp.kernel
def populate_hash_map(
    morton_codes: wp.array(dtype=int),
    active_indices: wp.array(dtype=int), 
    num_active: int,
    map_keys: wp.array(dtype=int),
    map_values: wp.array(dtype=int),
    capacity: int
):
    tid = wp.tid()
    if tid >= num_active:
        return
        
    pool_idx = active_indices[tid]
    code = morton_codes[pool_idx]
    
    map_insert(map_keys, map_values, capacity, code, pool_idx)

@wp.kernel
def compute_neighbors(
    morton_codes: wp.array(dtype=int),
    active_indices: wp.array(dtype=int),
    num_active: int,
    map_keys: wp.array(dtype=int),
    map_values: wp.array(dtype=int),
    map_capacity: int,
    out_neighbors: wp.array(dtype=int, ndim=2), 
    level: int,
    periodic_x: int, 
    periodic_y: int  
):
    tid = wp.tid()
    if tid >= num_active:
        return

    pool_idx = active_indices[tid]
    code = morton_codes[pool_idx]
    
    ix, iy = morton_decode(code)
    grid_dim = 1 << level
    
    # --- LEFT (Face 0) ---
    nx = ix - 1
    ny = iy
    if nx < 0:
        if periodic_x != 0: nx = grid_dim - 1
        else: nx = -1 
    
    n_idx = -1
    if nx >= 0:
        code_n = morton_encode(nx, ny)
        n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
        # --- CHECK IF MISSING ---
        if n_idx == -1:
             wp.printf("Error: Block %d (x=%d, y=%d) missing LEFT neighbor at x=%d, y=%d\n", pool_idx, ix, iy, nx, ny)

    out_neighbors[pool_idx, 0] = n_idx

    # --- RIGHT (Face 1) ---
    nx = ix + 1
    ny = iy
    if nx >= grid_dim:
        if periodic_x != 0: nx = 0
        else: nx = -1
    
    n_idx = -1
    if nx >= 0:
        code_n = morton_encode(nx, ny)
        n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
        # --- CHECK IF MISSING ---
        if n_idx == -1:
             wp.printf("Error: Block %d (x=%d, y=%d) missing RIGHT neighbor at x=%d, y=%d\n", pool_idx, ix, iy, nx, ny)

    out_neighbors[pool_idx, 1] = n_idx
    
    # --- BOTTOM (Face 2) ---
    nx = ix
    ny = iy - 1
    if ny < 0:
        if periodic_y != 0: ny = grid_dim - 1
        else: ny = -1
        
    n_idx = -1
    if ny >= 0:
        code_n = morton_encode(nx, ny)
        n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
        # --- CHECK IF MISSING ---
        if n_idx == -1:
             wp.printf("Error: Block %d (x=%d, y=%d) missing BOTTOM neighbor at x=%d, y=%d\n", pool_idx, ix, iy, nx, ny)

    out_neighbors[pool_idx, 2] = n_idx

    # --- TOP (Face 3) ---
    nx = ix
    ny = iy + 1
    if ny >= grid_dim:
        if periodic_y != 0: ny = 0
        else: ny = -1
        
    n_idx = -1
    if ny >= 0:
        code_n = morton_encode(nx, ny)
        n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
        # --- CHECK IF MISSING ---
        if n_idx == -1:
             wp.printf("Error: Block %d (x=%d, y=%d) missing TOP neighbor at x=%d, y=%d\n", pool_idx, ix, iy, nx, ny)

    out_neighbors[pool_idx, 3] = n_idx