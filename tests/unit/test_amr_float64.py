import pytest
import warp as wp
import numpy as np
from src.core.basis import Basis
from src.kernels.amr_kernels import prolongate_batch

@pytest.fixture(scope="module")
def device():
    wp.init()
    return "cpu"

def test_amr_accuracy_prolongation_float64(device):
    dtype = wp.vec4d
    scalar_dtype = np.float64
    N = 2
    basis = Basis(polynomial_degree=N, device=device, dtype=wp.float64) 
    Np = basis.Np
    
    # f(x, y) = x^2 + y
    def func(x, y):
        return x**2 + y
        
    # Parent Nodes
    nodes = basis.nodes_1d.numpy()
    x_p, y_p = np.meshgrid(nodes, nodes)
    q_parent_vals = func(x_p, y_p).flatten().astype(scalar_dtype)
    
    # Global Array: 0=Parent, 1-4=Children
    q_global_np = np.zeros((5, Np, 4), dtype=scalar_dtype)
    q_global_np[0, :, 0] = q_parent_vals # Set component 0
    
    q_global = wp.array(q_global_np, dtype=dtype, device=device)
    
    # Op List for Prolongation: [parent, child, quadrant]
    ops = []
    for i in range(4):
        ops.append([0, i + 1, i])
    
    ops_wp = wp.array(np.array(ops, dtype=np.int32), dtype=wp.int32, device=device)
    
    wp.launch(
        kernel=prolongate_batch,
        dim=4 * Np,
        inputs=[q_global, ops_wp, 4, basis.P_left, basis.P_right],
        device=device
    )
    
    # Verify
    res = q_global.numpy()
    
    for child_idx in range(4):
        if child_idx in [0, 2]: cx = 0.5 * (nodes - 1.0)
        else: cx = 0.5 * (nodes + 1.0)
        if child_idx in [0, 1]: cy = 0.5 * (nodes - 1.0)
        else: cy = 0.5 * (nodes + 1.0)
            
        xc, yc = np.meshgrid(cx, cy)
        q_expected = func(xc, yc).flatten().astype(scalar_dtype)
        
        q_child = res[child_idx + 1, :, 0]
        max_err = np.max(np.abs(q_child - q_expected))
        
        print(f"Child {child_idx} max error: {max_err}")
        assert max_err < 1e-10

def test_amr_conservation_restriction_float64(device):
    dtype = wp.vec4d
    scalar_dtype = np.float64
    N = 2
    basis = Basis(polynomial_degree=N, device=device, dtype=wp.float64)
    Np = basis.Np
    
    # Global Array: 0=Parent, 1-4=Children
    q_global_np = np.zeros((5, Np, 4), dtype=scalar_dtype)
    
    # Initialize children with random noise
    np.random.seed(42)
    for i in range(1, 5):
        q_global_np[i, :, 0] = np.random.rand(Np).astype(scalar_dtype)
        
    q_global = wp.array(q_global_np, dtype=dtype, device=device)
    
    # Op List: [child, parent, quadrant]
    ops = []
    for i in range(4):
        ops.append([i + 1, 0, i]) # Children 1..4 restrict to Parent 0
        
    ops_wp = wp.array(np.array(ops, dtype=np.int32), dtype=wp.int32, device=device)
    
    from src.kernels.amr_kernels import restrict_batch, zero_blocks
    
    # Zero parent
    wp.launch(kernel=zero_blocks, dim=(1, Np), inputs=[q_global, wp.array([0], dtype=wp.int32, device=device), 1], device=device)

    wp.launch(
        kernel=restrict_batch,
        dim=4 * Np,
        inputs=[q_global, ops_wp, 4, basis.R_left, basis.R_right],
        device=device
    )
    
    # Verify Conservation
    res = q_global.numpy()
    w1d = basis.weights_1d.numpy()
    w2d = np.kron(w1d, w1d)
    
    int_parent = np.sum(res[0, :, 0] * w2d)
    
    int_children_sum = 0.0
    for i in range(4):
        int_children_sum += np.sum(res[i + 1, :, 0] * w2d) * 0.25
        
    err = abs(int_parent - int_children_sum)
    print(f"Restriction conservation error: {err}")
    assert err < 1e-12
