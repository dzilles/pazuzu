import warp as wp

@wp.struct
class EquationParams32:
    gamma: wp.float32
    rho_floor: wp.float32
    p_floor: wp.float32
    half: wp.float32
    one: wp.float32
    # Freestream / Background state
    rho_inf: wp.float32
    u_inf: wp.float32
    v_inf: wp.float32
    p_inf: wp.float32

@wp.struct
class EquationParams64:
    gamma: wp.float64
    rho_floor: wp.float64
    p_floor: wp.float64
    half: wp.float64
    one: wp.float64
    # Freestream / Background state
    rho_inf: wp.float64
    u_inf: wp.float64
    v_inf: wp.float64
    p_inf: wp.float64
