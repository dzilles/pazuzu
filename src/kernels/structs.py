import warp as wp

@wp.struct
class EquationParams32:
    gamma: wp.float32
    rho_floor: wp.float32
    p_floor: wp.float32
    half: wp.float32
    one: wp.float32

@wp.struct
class EquationParams64:
    gamma: wp.float64
    rho_floor: wp.float64
    p_floor: wp.float64
    half: wp.float64
    one: wp.float64
