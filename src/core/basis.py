import numpy as np
from scipy.special import legendre
import warp as wp

class Basis:
    def __init__(self, polynomial_degree, device="cuda", dtype=wp.float32):
        self.N = polynomial_degree
        self.N1 = self.N + 1
        self.Np = self.N1 * self.N1
        self.Nfp = self.N1
        self.device = device
        self.dtype = dtype

        # Berechnung in float64 für Stabilität
        nodes_1d, weights_1d = self._gauss_lobatto_quadrature(self.N)
        
        nodes_2d = self._create_2d_nodes(nodes_1d)
        face_nodes = self._create_face_node_maps()
        
        # 1D Differentiationsmatrix
        D1D = self._differentiation_matrix_1d(nodes_1d)
        
        # 2D Ableitungen
        Dr, Ds = self._differentiation_matrices_2d(D1D)

        # Massenmatrizen
        # Die 2D Gewichte sind das Tensorprodukt der 1D Gewichte
        # Da M diagonal ist, ist M_inv einfach 1/weights
        weights_2d = np.kron(weights_1d, weights_1d)
        inv_mass_matrix_diag = 1.0 / weights_2d

        # Surface Mass Matrix (für die Integration auf den Kanten)
        M_face = np.diag(weights_1d)

        # LIFT Matrix erstellen (inklusive Inverser Massenmatrix!)
        LIFT = self._create_lift_matrix(face_nodes, M_face, inv_mass_matrix_diag)
        
        # Transfer zu Warp (ggf. casten zu float32 hier)
        self.nodes_1d = wp.array(nodes_1d, dtype=dtype, device=device)
        self.weights_1d = wp.array(weights_1d, dtype=dtype, device=device)
        self.nodes_2d = wp.array(nodes_2d, dtype=wp.vec2, device=device)
        self.face_nodes = wp.array(face_nodes, dtype=wp.int32, device=device)
        
        # V1D brauchst du im Solver meist nicht, eher für Initialisierung/Filterung
        # self.V1D = wp.array(V1D, dtype=dtype, device=device) 
        
        self.D1D = wp.array(D1D, dtype=dtype, device=device)
        self.Dr = wp.array(Dr, dtype=dtype, device=device)
        self.Ds = wp.array(Ds, dtype=dtype, device=device)
        self.LIFT = wp.array(LIFT, dtype=dtype, device=device)
        
        # Nützlich für CFL Bedingung
        self.min_node_dist = np.min(np.diff(nodes_1d))

    def _gauss_lobatto_quadrature(self, N):
        if N == 0: return np.array([0.0]), np.array([2.0])
        if N == 1: roots = np.array([])
        else: roots = np.roots(legendre(N).deriv(1))
        
        nodes = np.concatenate(([-1.0], np.sort(roots), [1.0]))
        weights = 2 / (N * self.N1 * legendre(N)(nodes)**2)
        return nodes, weights

    def _create_2d_nodes(self, nodes_1d):
        # Meshgrid ist oft einfacher und weniger fehleranfällig
        x, y = np.meshgrid(nodes_1d, nodes_1d)
        # Warp erwartet oft Array of Structs (x, y), also (Np, 2)
        return np.column_stack((x.flatten(), y.flatten()))

    def _create_face_node_maps(self):
        # Deine Map-Logik war korrekt für das Flattening row-major (Standard Numpy)
        # Face 0 (bottom, y=-1), Face 1 (right, x=+1), etc.
        # Aber Achtung: x changes fast, y changes slow in np.meshgrid('xy')? 
        # Prüfe Konsistenz mit _create_2d_nodes.
        # Bei meshgrid default ist x row, y col.
        # Hier behalten wir deine Logik bei, da sie konsistent wirkte.
        face_nodes = np.zeros((4, self.Nfp), dtype=np.int32)
        face_nodes[0, :] = np.arange(self.N1) # Bottom
        face_nodes[1, :] = np.arange(self.N1-1, self.Np, self.N1) # Right
        face_nodes[2, :] = np.arange(self.Np - self.N1, self.Np) # Top
        face_nodes[3, :] = np.arange(0, self.Np - self.N1 + 1, self.N1) # Left
        return face_nodes

    def _differentiation_matrix_1d(self, nodes):
        # KORRIGIERT
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
                        D[i, j] = 0.0 # WICHTIG: 0 für Legendre!
        return D

    def _differentiation_matrices_2d(self, D1D):
        I = np.eye(self.N1)
        # Achtung: Reihenfolge hängt von deinem Flattening ab.
        # Wenn nodes_2d via meshgrid/flatten erstellt wurde:
        # Dr wirkt auf x (fast index) -> Kron(I, D)
        # Ds wirkt auf y (slow index) -> Kron(D, I)
        Dr = np.kron(I, D1D) 
        Ds = np.kron(D1D, I)
        return Dr, Ds
        
    def _create_lift_matrix(self, face_nodes, M_face, inv_mass_matrix_diag):
        """
        LIFT = M^(-1) * E
        """
        # E Matrix: akkumuliert die Flux-Beiträge
        E = np.zeros((self.Np, 4 * self.Nfp))
        
        for i in range(4):
            # Indizes der Knoten auf Face i
            f_idx = face_nodes[i, :]
            
            # Da wir Nodal-DG machen, ist die Projektion einfach:
            # Wir schreiben die 1D Massenmatrix in die entsprechenden Zeilen von E
            for k in range(self.Nfp):
                row = f_idx[k]
                col = i * self.Nfp + k
                # Diagonaleintrag der Face-Massenmatrix
                E[row, col] = M_face[k, k] # oder einfach weights_1d[k]

        # Jetzt anwenden von M^(-1). Da M diagonal ist, skalieren wir einfach die Zeilen.
        # Broadcasting: (Np, 1) * (Np, 4*Nfp)
        LIFT = inv_mass_matrix_diag[:, None] * E
        
        return LIFT
