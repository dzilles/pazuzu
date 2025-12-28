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
    slot = hash_int(key) % capacity
    start_slot = slot
    
    while True:
        # Check if slot is empty or already contains the key (update)
        # We rely on atomic CAS if parallel, but for Phase 1 we can assume 
        # one thread per block insertion if distinct keys.
        # For simplicity in this "populate" kernel, we assume distinct keys (Morton codes are unique).
        
        # In a race-free scenario (unique keys, one thread per key):
        existing = keys[slot]
        if existing == EMPTY:
            # Found empty slot
            keys[slot] = key
            values[slot] = value
            return
        
        if existing == key:
            # Update existing
            values[slot] = value
            return
            
        slot = (slot + 1) % capacity
        if slot == start_slot:
            # Full
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
    active_indices: wp.array(dtype=int), # Mapping: 0..N -> PoolIndex
    num_active: int,
    map_keys: wp.array(dtype=int),
    map_values: wp.array(dtype=int),
    capacity: int
):
    """
    Inserts active blocks into the hash map.
    Key = Morton Code
    Value = Pool Index
    """
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
    out_neighbors: wp.array(dtype=int, ndim=2), # (MAX_BLOCKS, 4)
    level: int,
    periodic_x: int, # boolean 0/1
    periodic_y: int  # boolean 0/1
):
    """
    For each active block, finds neighbors by looking up computed codes in the hash map.
    Handles periodic boundaries if enabled.
    """
    tid = wp.tid()
    if tid >= num_active:
        return

    pool_idx = active_indices[tid]
    code = morton_codes[pool_idx]
    
    # Decode to get (x, y)
    ix, iy = morton_decode(code)
    
    grid_dim = 1 << level
    
    # Helper to wrap coordinate
    # LEFT (Face 0? No, indices 0:Left, 1:Right, 2:Bottom, 3:Top)
    # Wait, check mapping. 
    # Current code:
    # 0: Left (x-1)
    # 1: Right (x+1)
    # 2: Bottom (y-1)
    # 3: Top (y+1)
    
    # --- LEFT ---
    nx = ix - 1
    ny = iy
    if nx < 0:
        if periodic_x != 0: nx = grid_dim - 1
        else: nx = -1 # Invalid
    
    n_idx = -1
    if nx >= 0:
        code_n = morton_encode(nx, ny)
        n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
    out_neighbors[pool_idx, 0] = n_idx

    # --- RIGHT ---
    nx = ix + 1
    ny = iy
    if nx >= grid_dim:
        if periodic_x != 0: nx = 0
        else: nx = -1
    
    n_idx = -1
    if nx >= 0: # -1 indicates OOB non-periodic
        code_n = morton_encode(nx, ny)
        n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
    out_neighbors[pool_idx, 1] = n_idx
    
    # --- BOTTOM ---
    nx = ix
    ny = iy - 1
    if ny < 0:
        if periodic_y != 0: ny = grid_dim - 1
        else: ny = -1
        
    n_idx = -1
    if ny >= 0:
        code_n = morton_encode(nx, ny)
        n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
    out_neighbors[pool_idx, 2] = n_idx

    # --- TOP ---
    nx = ix
    ny = iy + 1
    if ny >= grid_dim:
        if periodic_y != 0: ny = 0
        else: ny = -1
        
    n_idx = -1
    if ny >= 0:
        code_n = morton_encode(nx, ny)
        n_idx = map_lookup(map_keys, map_values, map_capacity, code_n)
    out_neighbors[pool_idx, 3] = n_idx
