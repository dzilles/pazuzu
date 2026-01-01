import warp as wp

@wp.func
def make_vec2_generic(x: wp.float32, y: wp.float32):
    return wp.vec2(x, y)

@wp.func # type: ignore # Warp handles overloading
def make_vec2_generic(x: wp.float64, y: wp.float64):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.vec2d(x, y)

@wp.func
def make_vec4_generic(x: wp.float32, y: wp.float32, z: wp.float32, w: wp.float32):
    return wp.vec4(x, y, z, w)

@wp.func # type: ignore # Warp handles overloading
def make_vec4_generic(x: wp.float64, y: wp.float64, z: wp.float64, w: wp.float64):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.vec4d(x, y, z, w)

@wp.func
def get_half_generic(template: wp.float32):
    return wp.float32(0.5)

@wp.func # type: ignore # Warp handles overloading
def get_half_generic(template: wp.float64):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.float64(0.5)

@wp.func
def get_one_generic(template: wp.float32):
    return wp.float32(1.0)

@wp.func # type: ignore # Warp handles overloading
def get_one_generic(template: wp.float64):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.float64(1.0)

@wp.func
def get_any_generic(template: wp.float32, val: wp.float32):
    return wp.float32(val)

@wp.func # type: ignore # Warp handles overloading
def get_any_generic(template: wp.float32, val: wp.float64):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.float32(val)

@wp.func # type: ignore # Warp handles overloading
def get_any_generic(template: wp.float32, val: int):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.float32(val)

@wp.func # type: ignore # Warp handles overloading
def get_any_generic(template: wp.float64, val: wp.float32):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.float64(val)

@wp.func # type: ignore # Warp handles overloading
def get_any_generic(template: wp.float64, val: wp.float64):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.float64(val)

@wp.func # type: ignore # Warp handles overloading
def get_any_generic(template: wp.float64, val: int):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    return wp.float64(val)

@wp.func
def set_vec4_generic(v: wp.vec4, i: int, val: wp.float32):
    res = v
    if i == 0:
        res = wp.vec4(val, v[1], v[2], v[3])
    elif i == 1:
        res = wp.vec4(v[0], val, v[2], v[3])
    elif i == 2:
        res = wp.vec4(v[0], v[1], val, v[3])
    elif i == 3:
        res = wp.vec4(v[0], v[1], v[2], val)
    return res

@wp.func
def set_vec4_generic(v: wp.vec4d, i: int, val: wp.float64):  # noqa: F811 # type: ignore # Warp handles overloading
    # type: ignore
    res = v
    if i == 0:
        res = wp.vec4d(val, v[1], v[2], v[3])
    elif i == 1:
        res = wp.vec4d(v[0], val, v[2], v[3])
    elif i == 2:
        res = wp.vec4d(v[0], v[1], val, v[3])
    elif i == 3:
        res = wp.vec4d(v[0], v[1], v[2], val)
    return res

@wp.func
def get_zero_vec4_generic(template: wp.vec4):
    return wp.vec4(0.0, 0.0, 0.0, 0.0)

@wp.func # type: ignore
def get_zero_vec4_generic(template: wp.vec4d): # noqa: F811
    # type: ignore
    return wp.vec4d(wp.float64(0.0), wp.float64(0.0), wp.float64(0.0), wp.float64(0.0))
