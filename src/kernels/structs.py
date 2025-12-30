import warp as wp

@wp.struct
class EquationParams32:
    gamma: wp.float32
    rho_floor: wp.float32
    p_floor: wp.float32
    half: wp.float32
    one: wp.float32
    # Viscous parameters
    mu: wp.float32
    prandtl: wp.float32
    cp: wp.float32
    gas_constant: wp.float32
    # Freestream / Background state
    rho_inf: wp.float32
    u_inf: wp.float32
    v_inf: wp.float32
    p_inf: wp.float32
    epsilon: wp.float32
    flux_type: wp.int32
    hllc_fallback: wp.int32

@wp.struct
class EquationParams64:
    gamma: wp.float64
    rho_floor: wp.float64
    p_floor: wp.float64
    half: wp.float64
    one: wp.float64
    # Viscous parameters
    mu: wp.float64
    prandtl: wp.float64
    cp: wp.float64
    gas_constant: wp.float64
    # Freestream / Background state
    rho_inf: wp.float64
    u_inf: wp.float64
    v_inf: wp.float64
    p_inf: wp.float64
    epsilon: wp.float64
    flux_type: wp.int32
    hllc_fallback: wp.int32

@wp.struct
class BoundaryState32:
    type: wp.int32
    v0: wp.float32
    v1: wp.float32
    v2: wp.float32
    v3: wp.float32

@wp.struct
class BoundaryState64:
    type: wp.int32
    v0: wp.float64
    v1: wp.float64
    v2: wp.float64
    v3: wp.float64

@wp.struct
class MortarInterface:
    coarse_idx: wp.int32      # Index of the coarse block
    fine_idx_1: wp.int32      # Index of the first fine block (e.g., left/bottom half)
    fine_idx_2: wp.int32      # Index of the second fine block (e.g., right/top half)
    face_idx: wp.int32        # The face index on the COARSE block (0-3)
