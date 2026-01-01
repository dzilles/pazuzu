import numpy as np
import pytest
import warp as wp
import os
import sys

# Add src to path to allow imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from core.basis import Basis

def test_initialization(device):
    N = 2
    basis = Basis(polynomial_degree=N, device=device)
    
    assert basis.N == N
    assert basis.N1 == N + 1
    assert basis.Np == (N + 1) ** 2
    assert basis.Nfp == N + 1

@pytest.mark.parametrize("N", [1, 2, 4, 8])
def test_nodes_1d_range(N, device):
    basis = Basis(polynomial_degree=N, device=device)
    nodes = basis.nodes_1d.numpy()
    
    assert len(nodes) == N + 1
    assert pytest.approx(nodes[0]) == -1.0
    assert pytest.approx(nodes[-1]) == 1.0
    # Check symmetry
    np.testing.assert_allclose(nodes, -nodes[::-1], atol=1e-14)

@pytest.mark.parametrize("N", [1, 2, 4, 8])
def test_weights_sum(N, device):
    """Weights in 1D should sum to 2 (integral of 1 from -1 to 1)"""
    basis = Basis(polynomial_degree=N, device=device)
    weights = basis.weights_1d.numpy()
    
    assert pytest.approx(np.sum(weights)) == 2.0

@pytest.mark.parametrize("N", [1, 2, 4, 8])
def test_differentiation_matrix(N, device):
    """Test differentiation matrix properties: D*x = 1 and D*1 = 0"""
    basis = Basis(polynomial_degree=N, device=device)
    D1D = basis.D1D.numpy()
    nodes = basis.nodes_1d.numpy()
    
    # 1. Test derivative of constant function should be 0: D * 1 = 0
    f_const = np.ones_like(nodes)
    df_const = D1D @ f_const
    expected_const = np.zeros_like(nodes)
    np.testing.assert_allclose(df_const, expected_const, atol=1e-6)

    # 2. Test derivative of linear function: D * x = 1
    f_linear = nodes.copy()
    df_linear = D1D @ f_linear
    expected_linear = np.ones_like(nodes)
    np.testing.assert_allclose(df_linear, expected_linear, atol=1e-6)

def test_2d_nodes_bounds(device):
    N = 1
    basis = Basis(polynomial_degree=N, device=device)
    nodes_2d = basis.nodes_2d.numpy()
    
    # Expect 4 nodes for N=1 (corners)
    assert nodes_2d.shape == (4, 2)
    assert np.all(nodes_2d >= -1.0)
    assert np.all(nodes_2d <= 1.0)

def test_face_nodes_indices(device):
    """Check if face nodes map to valid indices"""
    N = 2
    basis = Basis(polynomial_degree=N, device=device)
    face_nodes = basis.face_nodes.numpy()
    
    # 4 faces
    assert face_nodes.shape == (4, N+1)
    
    # Indices should be within [0, Np-1]
    assert np.all(face_nodes >= 0)
    assert np.all(face_nodes < basis.Np)
    
    # Specific check for Face 0 (Bottom, y=-1)
    # Nodes should be 0, 1, ..., N
    expected_face0 = np.arange(N + 1)
    np.testing.assert_array_equal(face_nodes[0], expected_face0)

@pytest.mark.parametrize("N", [1, 2, 4, 8])
def test_filter_matrix_computation(N, device):
    basis = Basis(polynomial_degree=N, device=device)
    basis.compute_filter_matrix(alpha=10.0, order=2)
    
    assert basis.filter_matrix is not None
    filter_mat = basis.filter_matrix.numpy()
    
    # Filter matrix should be Np x Np
    assert filter_mat.shape == (basis.Np, basis.Np)

def test_lift_matrix_shape(device):
    N = 2
    basis = Basis(polynomial_degree=N, device=device)
    LIFT = basis.LIFT.numpy()
    
    # LIFT should be (Np, 4 * Nfp)
    expected_shape = (basis.Np, 4 * basis.Nfp)
    assert LIFT.shape == expected_shape
    
    # It shouldn't be all zeros
    assert not np.all(LIFT == 0)