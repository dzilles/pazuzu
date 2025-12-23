import numpy as np
from scipy.special import legendre
import warp as wp

class Basis:
    """
    Manages the nodal Discontinuous Galerkin (DG) basis functions and operators on a reference quadrilateral element.
    
    This class handles the creation of Gauss-Lobatto-Legendre (GLL) nodes, differentiation matrices,
    mass matrices, and the LIFT operator used to project surface fluxes into the volume.
    
    Attributes:
        N (int): Polynomial degree of the basis.
        N1 (int): Number of nodes in 1D (N + 1).
        Np (int): Total number of nodes in 2D (N1 * N1).
        Nfp (int): Number of face points (N1).
        device (str): Compute device for Warp arrays ("cpu" or "cuda").
        dtype (wp.dtype): Data type for Warp arrays (e.g., wp.float32).
        nodes_1d (wp.array): 1D GLL node coordinates.
        weights_1d (wp.array): 1D GLL integration weights.
        nodes_2d (wp.array): 2D tensor product node coordinates.
        face_nodes (wp.array): Indices of nodes belonging to each of the 4 faces.
        D1D (wp.array): 1D differentiation matrix.
        Dr (wp.array): 2D differentiation matrix with respect to r (reference x).
        Ds (wp.array): 2D differentiation matrix with respect to s (reference y).
        LIFT (wp.array): Lift operator matrix for flux reconstruction.
        min_node_dist (float): Minimum distance between 1D nodes (used for CFL condition).
    """
    def __init__(self, polynomial_degree, device="cuda", dtype=wp.float32):
        """
        Initializes the DG Basis.

        Args:
            polynomial_degree (int): The order of the polynomial basis (N).
            device (str, optional): The Warp device to allocate arrays on. Defaults to "cuda".
            dtype (wp.dtype, optional): Floating point precision. Defaults to wp.float32.
        """
        self.N = polynomial_degree
        self.N1 = self.N + 1
        self.Np = self.N1 * self.N1
        self.Nfp = self.N1
        self.device = device
        self.dtype = dtype

        # Perform calculations in float64 on the host for numerical stability
        nodes_1d, weights_1d = self._gauss_lobatto_quadrature(self.N)
        
        nodes_2d = self._create_2d_nodes(nodes_1d)
        face_nodes = self._create_face_node_maps()
        
        # 1D Differentiation Matrix
        D1D = self._differentiation_matrix_1d(nodes_1d)
        
        # 2D Differentiation Matrices
        Dr, Ds = self._differentiation_matrices_2d(D1D)

        # Mass Matrices
        # The 2D weights are the tensor product of 1D weights.
        # Since the Mass Matrix M is diagonal for the GLL basis, M_inv is simply 1/weights.
        weights_2d = np.kron(weights_1d, weights_1d)
        inv_mass_matrix_diag = 1.0 / weights_2d

        # Surface Mass Matrix (for integration along edges)
        M_face = np.diag(weights_1d)

        # Create LIFT Matrix (includes inverse mass matrix application)
        LIFT = self._create_lift_matrix(face_nodes, M_face, inv_mass_matrix_diag)
        
        # Transfer data to Warp arrays (cast to target dtype here)
        vec2_type = wp.vec2d if dtype == wp.float64 else wp.vec2
        
        self.nodes_1d = wp.array(nodes_1d, dtype=dtype, device=device)
        self.weights_1d = wp.array(weights_1d, dtype=dtype, device=device)
        self.nodes_2d = wp.array(nodes_2d, dtype=vec2_type, device=device)
        self.face_nodes = wp.array(face_nodes, dtype=wp.int32, device=device)
        
        self.D1D = wp.array(D1D, dtype=dtype, device=device)
        self.Dr = wp.array(Dr, dtype=dtype, device=device)
        self.Ds = wp.array(Ds, dtype=dtype, device=device)
        self.LIFT = wp.array(LIFT, dtype=dtype, device=device)
        
        # Useful for CFL condition estimation
        self.min_node_dist = np.min(np.diff(nodes_1d))
        
        # Filter Matrix (lazy initialization)
        self.filter_matrix = None

    def compute_filter_matrix(self, alpha, order):
        """
        Computes and initializes the spectral viscosity filter matrix.
        
        Constructs the 1D filter matrix F1d = V * Lambda * V^(-1), 
        where V is the Vandermonde matrix of Legendre polynomials,
        and Lambda is the diagonal filter matrix.
        Then computes the 2D tensor product F2d = F1d (x) F1d.
        
        Args:
            alpha (float): Filter strength parameter.
            order (int): Filter order parameter.
        """
        # 1. Construct 1D Vandermonde Matrix V
        # V[i, j] = P_j(x_i) where x_i are the GLL nodes
        nodes = self.nodes_1d.numpy()
        V = np.zeros((self.N1, self.N1))
        
        # Use unnormalized Legendre polynomials as the basis functions P_j
        for j in range(self.N1):
            P_j = legendre(j)
            V[:, j] = P_j(nodes)
            
        # 2. Construct Diagonal Filter Matrix Lambda
        # sigma_k = exp(-alpha * (k/N)^order)
        Lambda = np.zeros((self.N1, self.N1))
        for k in range(self.N1):
            sigma = np.exp(-alpha * ((k / self.N) ** order))
            Lambda[k, k] = sigma
            
        # 3. Compute 1D Filter Matrix F1d = V * Lambda * V^(-1)
        V_inv = np.linalg.inv(V)
        F1d = V @ Lambda @ V_inv
        
        # 4. Compute 2D Filter Matrix via Tensor Product
        # F2d = F1d (kron) F1d
        F2d = np.kron(F1d, F1d)
        
        # 5. Upload to Device
        self.filter_matrix = wp.array(F2d, dtype=self.dtype, device=self.device)

    def _gauss_lobatto_quadrature(self, N):
        """
        Computes Gauss-Lobatto-Legendre (GLL) nodes and weights.

        Args:
            N (int): Polynomial degree.

        Returns:
            tuple: A tuple (nodes, weights) containing numpy arrays of the quadrature points and weights.
        """
        if N == 0: return np.array([0.0]), np.array([2.0])
        if N == 1: roots = np.array([])
        else: roots = np.roots(legendre(N).deriv(1))
        
        nodes = np.concatenate(([-1.0], np.sort(roots), [1.0]))
        weights = 2 / (N * self.N1 * legendre(N)(nodes)**2)
        return nodes, weights

    def _create_2d_nodes(self, nodes_1d):
        """
        Creates 2D nodal coordinates via tensor product.

        Args:
            nodes_1d (np.array): 1D GLL nodes.

        Returns:
            np.array: A (Np, 2) array containing (r, s) coordinates for all 2D nodes.
        """
        # Meshgrid is simpler and less error-prone
        x, y = np.meshgrid(nodes_1d, nodes_1d)
        # Warp expects Array of Structs (x, y), so we stack them: (Np, 2)
        return np.column_stack((x.flatten(), y.flatten()))

    def _create_face_node_maps(self):
        """
        Generates mappings from face indices to global 2D node indices.
        
        Assumes standard tensor-product ordering (row-major flattening).
        Face 0: Bottom (y=-1)
        Face 1: Right (x=+1)
        Face 2: Top (y=+1)
        Face 3: Left (x=-1)

        Returns:
            np.array: A (4, Nfp) integer array containing node indices for each face.
        """
        face_nodes = np.zeros((4, self.Nfp), dtype=np.int32)
        face_nodes[0, :] = np.arange(self.N1) # Bottom
        face_nodes[1, :] = np.arange(self.N1-1, self.Np, self.N1) # Right
        face_nodes[2, :] = np.arange(self.Np - self.N1, self.Np) # Top
        face_nodes[3, :] = np.arange(0, self.Np - self.N1 + 1, self.N1) # Left
        return face_nodes

    def _differentiation_matrix_1d(self, nodes):
        """
        Computes the 1D Lagrange differentiation matrix for the GLL points.

        Args:
            nodes (np.array): 1D GLL nodes.

        Returns:
            np.array: The (N1, N1) differentiation matrix D, where D[i,j] = dL_j/dx(x_i).
        """
        D = np.zeros((self.N1, self.N1))
        for i in range(self.N1):
            for j in range(self.N1):
                if i != j:
                    Ln_i = legendre(self.N)(nodes[i])
                    Ln_j = legendre(self.N)(nodes[j])
                    D[i, j] = Ln_i / (Ln_j * (nodes[i] - nodes[j]))
                else:
                    if i == 0:
                        D[i, j] = -self.N * self.N1 / 4.0
                    elif i == self.N:
                        D[i, j] = self.N * self.N1 / 4.0
                    else:
                        D[i, j] = 0.0 # Important: 0 for Legendre GLL internal nodes
        return D

    def _differentiation_matrices_2d(self, D1D):
        """
        Constructs 2D differentiation matrices using the Kronecker product.

        Args:
            D1D (np.array): The 1D differentiation matrix.

        Returns:
            tuple: (Dr, Ds), the differentiation matrices for r and s directions respectively.
        """
        I = np.eye(self.N1)
        # Note: The order depends on the flattening strategy.
        # Since nodes_2d was created via meshgrid(x, y) then flattened:
        # x varies fast (inner loop), y varies slow (outer loop).
        # Dr acts on x -> Kron(I, D)
        # Ds acts on y -> Kron(D, I)
        Dr = np.kron(I, D1D) 
        Ds = np.kron(D1D, I)
        return Dr, Ds
        
    def _create_lift_matrix(self, face_nodes, M_face, inv_mass_matrix_diag):
        """
        Constructs the LIFT matrix, which maps surface flux terms to volume residuals.
        
        The operator is defined as LIFT = M^(-1) * E, where E is the flux accumulation matrix.

        Args:
            face_nodes (np.array): Map of face node indices.
            M_face (np.array): 1D Surface Mass Matrix (diagonal).
            inv_mass_matrix_diag (np.array): Inverse of the 2D diagonal Mass Matrix.

        Returns:
            np.array: The (Np, 4*Nfp) LIFT matrix.
        """
        # E Matrix: accumulates flux contributions from faces to volume nodes
        E = np.zeros((self.Np, 4 * self.Nfp))
        
        for i in range(4):
            # Indices of nodes on Face i
            f_idx = face_nodes[i, :]
            
            # For Nodal DG, projection is simplified:
            # We write the 1D mass matrix weights into the corresponding rows of E
            for k in range(self.Nfp):
                row = f_idx[k]
                col = i * self.Nfp + k
                # Diagonal entry of the face mass matrix (or simply weights_1d[k])
                E[row, col] = M_face[k, k] 

        # Apply M^(-1). Since M is diagonal, we simply scale the rows.
        # Broadcasting: (Np, 1) * (Np, 4*Nfp)
        LIFT = inv_mass_matrix_diag[:, None] * E
        
        return LIFT
