import unittest
import numpy as np
import warp as wp
import sys
import os

# Add src to path to allow imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src')))

from core.basis import Basis

class TestBasis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Initialize Warp once
        if not wp.is_cuda_available():
            cls.device = "cpu"
        else:
            cls.device = "cpu" # Force CPU for unit tests to be safe/simpler
        wp.init()

    def test_initialization(self):
        N = 2
        basis = Basis(polynomial_degree=N, device=self.device)
        
        self.assertEqual(basis.N, N)
        self.assertEqual(basis.N1, N + 1)
        self.assertEqual(basis.Np, (N + 1) ** 2)
        self.assertEqual(basis.Nfp, N + 1)

    def test_nodes_1d_range(self):
        N = 3
        basis = Basis(polynomial_degree=N, device=self.device)
        nodes = basis.nodes_1d.numpy()
        
        self.assertEqual(len(nodes), N + 1)
        self.assertAlmostEqual(nodes[0], -1.0)
        self.assertAlmostEqual(nodes[-1], 1.0)
        # Check symmetry
        np.testing.assert_allclose(nodes, -nodes[::-1], atol=1e-14)

    def test_weights_sum(self):
        """Weights in 1D should sum to 2 (integral of 1 from -1 to 1)"""
        N = 4
        basis = Basis(polynomial_degree=N, device=self.device)
        weights = basis.weights_1d.numpy()
        
        self.assertAlmostEqual(np.sum(weights), 2.0)

    def test_differentiation_matrix(self):
        """Test derivative of f(x) = x should be 1"""
        N = 2
        basis = Basis(polynomial_degree=N, device=self.device)
        D1D = basis.D1D.numpy()
        nodes = basis.nodes_1d.numpy()
        
        # f = x
        f = nodes.copy()
        # df/dx = D * f
        df = D1D @ f
        
        expected = np.ones_like(nodes)
        np.testing.assert_allclose(df, expected, atol=1e-12)

    def test_differentiation_matrix_constant(self):
        """Test derivative of constant function should be 0"""
        N = 2
        basis = Basis(polynomial_degree=N, device=self.device)
        D1D = basis.D1D.numpy()
        nodes = basis.nodes_1d.numpy()
        
        # f = 1
        f = np.ones_like(nodes)
        df = D1D @ f
        
        expected = np.zeros_like(nodes)
        np.testing.assert_allclose(df, expected, atol=1e-12)

    def test_2d_nodes_bounds(self):
        N = 1
        basis = Basis(polynomial_degree=N, device=self.device)
        nodes_2d = basis.nodes_2d.numpy()
        
        # Expect 4 nodes for N=1 (corners)
        self.assertEqual(nodes_2d.shape, (4, 2))
        self.assertTrue(np.all(nodes_2d >= -1.0))
        self.assertTrue(np.all(nodes_2d <= 1.0))

    def test_face_nodes_indices(self):
        """Check if face nodes map to valid indices"""
        N = 2
        basis = Basis(polynomial_degree=N, device=self.device)
        face_nodes = basis.face_nodes.numpy()
        
        # 4 faces
        self.assertEqual(face_nodes.shape, (4, N+1))
        
        # Indices should be within [0, Np-1]
        self.assertTrue(np.all(face_nodes >= 0))
        self.assertTrue(np.all(face_nodes < basis.Np))
        
        # Specific check for Face 0 (Bottom, y=-1)
        # Nodes should be 0, 1, ..., N
        expected_face0 = np.arange(N + 1)
        np.testing.assert_array_equal(face_nodes[0], expected_face0)

    def test_filter_matrix_computation(self):
        N = 3
        basis = Basis(polynomial_degree=N, device=self.device)
        basis.compute_filter_matrix(alpha=10.0, order=2)
        
        self.assertIsNotNone(basis.filter_matrix)
        filter_mat = basis.filter_matrix.numpy()
        
        # Filter matrix should be Np x Np
        self.assertEqual(filter_mat.shape, (basis.Np, basis.Np))

    def test_lift_matrix_shape(self):
        N = 2
        basis = Basis(polynomial_degree=N, device=self.device)
        LIFT = basis.LIFT.numpy()
        
        # LIFT should be (Np, 4 * Nfp)
        expected_shape = (basis.Np, 4 * basis.Nfp)
        self.assertEqual(LIFT.shape, expected_shape)
        
        # It shouldn't be all zeros
        self.assertFalse(np.all(LIFT == 0))

if __name__ == '__main__':
    unittest.main()
