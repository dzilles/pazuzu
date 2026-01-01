import pytest
import warp as wp
import numpy as np
from src.kernels.grid_kernels import morton_decode, part1by1

# Wrap the warp function in a kernel to test it
@wp.kernel
def kernel_decode_test(codes: wp.array(dtype=wp.int32), 
                       levels: wp.array(dtype=wp.int32),
                       out_x: wp.array(dtype=int), 
                       out_y: wp.array(dtype=int)):
    tid = wp.tid()
    c = codes[tid]
    l = levels[tid]
    x, y = morton_decode(c, l)
    out_x[tid] = x
    out_y[tid] = y

@pytest.mark.parametrize("device", ["cpu"])
def test_morton_decode_implementation(device):
    """
    Exhaustively test Morton decoding for a small range of coordinates.
    """
    wp.init()
    
    # (code, level, expected_x, expected_y)
    test_cases = [
        (1, 0, 0, 0),
        (4, 1, 0, 0),
        (5, 1, 1, 0),
        (6, 1, 0, 1),
        (7, 1, 1, 1),
        (28, 2, 2, 2), # 16 | 12
        (31, 2, 3, 3)  # 16 | 15
    ]
    
    codes_np = np.array([c[0] for c in test_cases], dtype=np.int32)
    levels_np = np.array([c[1] for c in test_cases], dtype=np.int32)
    
    # Warp Arrays
    codes_wp = wp.array(codes_np, dtype=wp.int32, device=device)
    levels_wp = wp.array(levels_np, dtype=wp.int32, device=device)
    out_x_wp = wp.zeros_like(codes_wp)
    out_y_wp = wp.zeros_like(codes_wp)
    
    wp.launch(
        kernel=kernel_decode_test,
        dim=len(codes_np),
        inputs=[codes_wp, levels_wp, out_x_wp, out_y_wp],
        device=device
    )
    # ...
    
    res_x = out_x_wp.numpy()
    res_y = out_y_wp.numpy()
    
    for i, (code, level, ex_x, ex_y) in enumerate(test_cases):
        assert res_x[i] == ex_x, f"Code {code} decoded X={res_x[i]}, expected {ex_x}"
        assert res_y[i] == ex_y, f"Code {code} decoded Y={res_y[i]}, expected {ex_y}"
