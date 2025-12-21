from .base_solver import BaseSolver
from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.physics import equations as eq
from src.kernels import euler_kernels as wk # Renamed import
import numpy as np
import warp as wp

class Euler2DSolver(BaseSolver):
    """
    Solves the 2D Euler equations using a Nodal Discontinuous Galerkin method.

    This solver manages the setup and execution of the spatial discretization for the Euler equations.
    It handles the initialization of the state vector Q, computation of geometric factors for the mesh,
    setup of boundary conditions, and the computation of the right-hand side (RHS) using Warp kernels.

    Attributes:
        geo_factors_host (np.array): Geometric factors for volume integration (host).
        face_geo_factors_host (np.array): Geometric factors for surface integration (host).
        x_host (np.array): Physical x-coordinates of nodes (host).
        y_host (np.array): Physical y-coordinates of nodes (host).
        bc_mask_host (np.array): Boundary condition masks (host).
        mesh (Mesh): The computational mesh.
        basis (Basis): The DG basis.
        Q (wp.array): The state vector [rho, rho*u, rho*v, E] on device.
        rhs (wp.array): The RHS buffer on device.
        face_geo_factors (wp.array): Geometric factors for faces on device.
        bc_mask (wp.array): Boundary condition masks on device.
        x (wp.array): Physical x-coordinates on device.
        y (wp.array): Physical y-coordinates on device.
        M_inv_diag (wp.array): Inverse of the diagonal mass matrix on device.
        max_wave_speed (wp.array): Buffer for maximum wave speed reduction.
    """
    def __init__(self, mesh, basis, config=None):
        """
        Initializes the Euler 2D Solver.

        Args:
            mesh (Mesh): The computational mesh.
            basis (Basis): The DG basis.
            config (dict, optional): Configuration dictionary. Defaults to None.
        """
        super().__init__(mesh, basis, config)
        
        # --- Host-side data (for setup) ---
        # Volume metrics: dx/dr, dx/ds, dy/dr, dy/ds, J
        self.geo_factors_host = np.zeros((mesh.num_elements, 5), dtype=np.float32)
        # Face metrics: nx, ny, J_surf
        self.face_geo_factors_host = np.zeros((mesh.num_elements, 4, 3), dtype=np.float32)
        # Physical coordinates
        self.x_host = np.zeros((mesh.num_elements, basis.Np), dtype=np.float32)
        self.y_host = np.zeros((mesh.num_elements, basis.Np), dtype=np.float32)
        
        # Calculate geometric factors (Metrics & Jacobians)
        self._calculate_geometric_factors()
        
        # Map reference nodes to physical space
        self._map_nodes_to_physical_space()

        # Setup Boundary Conditions mask
        self._setup_boundary_conditions()

        # --- Device-side data (for computation) ---
        
        # Update Mesh Metric Arrays on Device (ensure they are synced)
        self.mesh.rx = wp.array(self.mesh.rx_host, dtype=wp.float32, device=self.device)
        self.mesh.ry = wp.array(self.mesh.ry_host, dtype=wp.float32, device=self.device)
        self.mesh.sx = wp.array(self.mesh.sx_host, dtype=wp.float32, device=self.device)
        self.mesh.sy = wp.array(self.mesh.sy_host, dtype=wp.float32, device=self.device)
        self.mesh.J  = wp.array(self.mesh.J_host, dtype=wp.float32, device=self.device)
        
        self.face_geo_factors = wp.array(self.face_geo_factors_host, dtype=wp.float32, device=self.device)
        self.bc_mask = wp.array(self.bc_mask_host, dtype=wp.int32, device=self.device)
        
        self.x = wp.array(self.x_host, dtype=wp.float32, device=self.device)
        self.y = wp.array(self.y_host, dtype=wp.float32, device=self.device)

        # Inverse of Mass Matrix (diagonal)
        weights_2d = np.kron(self.basis.weights_1d.numpy(), self.basis.weights_1d.numpy())
        M_inv_diag_host = 1.0 / weights_2d
        self.M_inv_diag = wp.array(M_inv_diag_host, dtype=wp.float32, device=self.device)
        
        # Buffer for global max wave speed reduction
        self.max_wave_speed = wp.zeros(1, dtype=wp.float32, device=self.device)

        # Ramping parameter for inlet BC
        self.ramp_time = 1.0
        if config and 'simulation' in config:
            self.ramp_time = config['simulation'].get('ramp_time', 1.0)
            
    def initialize(self, initial_condition_func):
        """
        Initializes the state vector Q using the provided initial condition function.

        Args:
            initial_condition_func (callable): Function f(x, y) -> (rho, u, v, p).
        """
        Q_host = np.zeros((self.mesh.num_elements, self.basis.Np, 4), dtype=np.float32)
        self._set_initial_conditions(Q_host, initial_condition_func)
        self.Q = wp.array(Q_host, dtype=wp.vec4, device=self.device)
        self.rhs = wp.zeros_like(self.Q)

    def compute_rhs(self, t, dt, q, rhs):
        """
        Computes the Right-Hand Side (RHS) of the semi-discrete Euler equations.
        
        RHS = -M^(-1) * (VolumeIntegral + SurfaceIntegral)
        However, in this strong form implementation, we compute:
        RHS = - (div(F)) + LIFT(F* - F_n)
        And the inverse mass matrix is applied implicitly or within the LIFT.
        Actually, the kernels currently compute terms directly. 
        Volume term: -div(F)
        Surface term: LIFT * (FluxJump) / J
        
        Args:
            t (float): Current simulation time.
            dt (float): Current time step.
            q (wp.array): Input state vector.
            rhs (wp.array): Output RHS buffer.
        """
        rhs.zero_()
        
        # --- Volume Integral ---
        # Computes -div(F)
        wp.launch(
            kernel=wk.compute_volume_term,
            dim=(self.mesh.num_elements, self.basis.Np),
            inputs=[
                q,
                rhs,
                self.basis.Dr,
                self.basis.Ds,
                self.mesh.rx,
                self.mesh.ry,
                self.mesh.sx,
                self.mesh.sy,
                self.basis.Np
            ],
            device=self.device
        )
        
        # --- Surface Integral ---
        # Adds LIFT * (NumericalFlux - NormalFlux)
        wp.launch(
            kernel=wk.compute_surface_term,
            dim=self.mesh.num_elements,
            inputs=[
                q,
                rhs,
                self.mesh.connectivity,
                self.mesh.connectivity_face_indices,
                self.basis.face_nodes,
                self.basis.LIFT,
                self.face_geo_factors, 
                self.mesh.J,
                self.bc_mask,
                self.basis.Nfp,
                t,
                self.ramp_time
            ],
            device=self.device
        )

    def calculate_dt(self, CFL):
        """
        Calculates the maximum stable time step size based on the CFL condition.

        Args:
            CFL (float): The desired CFL number.

        Returns:
            float: The time step size dt.
        """
        self.max_wave_speed.zero_()
        wp.launch(
            kernel=wk.compute_max_wave_speed,
            dim=(self.mesh.num_elements, self.basis.Np),
            inputs=[self.Q, self.max_wave_speed],
            device=self.device
        )
        max_speed = self.max_wave_speed.numpy()[0]
        
        # Safety check to avoid division by zero if fluid is at rest and T=0 (unlikely)
        if max_speed < 1.0e-6:
            max_speed = 1.0 

        # Characteristic length / Max Speed
        dt = CFL * self.mesh.dx / max_speed
        return dt

    def _setup_boundary_conditions(self):
        """
        Parses the configuration to setup the boundary condition mask.
        Maps physical tags from the mesh to solver-specific BC IDs.
        """
        self.bc_mask_host = np.zeros((self.mesh.num_elements, 4), dtype=np.int32)
        
        if self.config and 'boundaries' in self.config:
            BC_WALL = 1
            BC_FARFIELD = 2
            BC_INLET = 3
            BC_OUTLET = 4
            
            tag_to_bc = {}
            for name, bc_conf in self.config['boundaries'].items():
                tag = -1
                # Try to map name to tag using mesh info
                if hasattr(self.mesh, 'physical_groups') and name in self.mesh.physical_groups:
                    tag = self.mesh.physical_groups[name]
                else:
                    try: tag = int(name) 
                    except: pass
                
                if tag != -1:
                    bc_type = bc_conf.get('type')
                    if bc_type == "slip_wall": tag_to_bc[tag] = BC_WALL
                    elif bc_type in ["farfield", "outflow"]: tag_to_bc[tag] = BC_FARFIELD
                    elif bc_type == "inlet": tag_to_bc[tag] = BC_INLET
                    elif bc_type == "outlet": tag_to_bc[tag] = BC_OUTLET
            
            # Apply to mask
            for e in range(self.mesh.num_elements):
                for f in range(4):
                    tag = self.mesh.boundary_tags_host[e, f]
                    if tag > 0:
                        # Default to Farfield if tag exists but type not specified
                        self.bc_mask_host[e, f] = tag_to_bc.get(tag, BC_FARFIELD)

    def _calculate_geometric_factors(self):
        """
        Calculates geometric metrics (Jacobians) for all elements.
        Currently assumes affine (linear) quadrilateral elements.
        """
        # Derivatives of bilinear shape functions at the center (r=0, s=0)
        dNd_r_center = 0.25 * np.array([-1, 1, 1, -1])
        dNd_s_center = 0.25 * np.array([-1, -1, 1, 1])

        # Precomputed derivative vectors for faces to handle general quads correctly
        # Face 0 (Bottom, s=-1): dN/dr at r=0, s=-1 -> [-0.5, 0.5, 0.0, 0.0]
        dNd_r_f0 = np.array([-0.5, 0.5, 0.0, 0.0])
        
        # Face 1 (Right, r=1): dN/ds at r=1, s=0 -> [0.0, -0.5, 0.5, 0.0]
        dNd_s_f1 = np.array([0.0, -0.5, 0.5, 0.0])
        
        # Face 2 (Top, s=1): dN/dr at r=0, s=1 -> [0.0, 0.0, 0.5, -0.5]
        dNd_r_f2 = np.array([0.0, 0.0, 0.5, -0.5])
        
        # Face 3 (Left, r=-1): dN/ds at r=-1, s=0 -> [-0.5, 0.0, 0.0, 0.5]
        dNd_s_f3 = np.array([-0.5, 0.0, 0.0, 0.5])
        
        for i in range(self.mesh.num_elements):
            v = self.mesh.vertices_host[i, :, :] 
            dx_dr = dNd_r_center @ v[:, 0]
            dx_ds = dNd_s_center @ v[:, 0]
            dy_dr = dNd_r_center @ v[:, 1]
            dy_ds = dNd_s_center @ v[:, 1]
            
            # Jacobian Determinant J = det(dx/dr dx/ds; dy/dr dy/ds)
            # Actually J = dx/dr * dy/ds - dx/ds * dy/dr
            J = dx_dr * dy_ds - dx_ds * dy_dr
            self.mesh.J_host[i] = J
            
            # Inverse Jacobian components
            dr_dx =  dy_ds / J
            dr_dy = -dx_ds / J
            ds_dx = -dy_dr / J
            ds_dy =  dx_dr / J
            
            self.mesh.rx_host[i] = dr_dx
            self.mesh.ry_host[i] = dr_dy
            self.mesh.sx_host[i] = ds_dx
            self.mesh.sy_host[i] = ds_dy

            self.geo_factors_host[i, 0] = dx_dr
            self.geo_factors_host[i, 1] = dx_ds
            self.geo_factors_host[i, 2] = dy_dr
            self.geo_factors_host[i, 3] = dy_ds
            self.geo_factors_host[i, 4] = J
            
            # --- Face Normals and Surface Jacobians ---
            # Normals are outward pointing.
            # J_face is the scaling factor (length of the edge relative to reference interval length 2)
            
            # Face 0 (bottom): s=-1, r in [-1, 1]. Tangent vector T = dx/dr. Normal N = (dy/dr, -dx/dr).
            dx_dr_0 = dNd_r_f0 @ v[:, 0]
            dy_dr_0 = dNd_r_f0 @ v[:, 1]
            nx, ny = dy_dr_0, -dx_dr_0
            J_face = np.sqrt(nx**2 + ny**2) # Length / 2
            self.face_geo_factors_host[i, 0, :] = [nx/J_face, ny/J_face, J_face]

            # Face 1 (right): r=1, s in [-1, 1]. Tangent T = dx/ds. Normal N = (dy/ds, -dx/ds).
            dx_ds_1 = dNd_s_f1 @ v[:, 0]
            dy_ds_1 = dNd_s_f1 @ v[:, 1]
            nx, ny = dy_ds_1, -dx_ds_1
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 1, :] = [nx/J_face, ny/J_face, J_face]
            
            # Face 2 (top): s=1, r in [1, -1] (reverse!). Tangent T = -dx/dr. Normal N = (-dy/dr, dx/dr).
            # This points outwards (Opposite to bottom normal direction relative to derivative)
            dx_dr_2 = dNd_r_f2 @ v[:, 0]
            dy_dr_2 = dNd_r_f2 @ v[:, 1]
            nx, ny = -dy_dr_2, dx_dr_2 
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 2, :] = [nx/J_face, ny/J_face, J_face]

            # Face 3 (left): r=-1, s in [1, -1] (reverse!). Tangent T = -dx/ds. Normal N = (-dy/ds, dx/ds).
            dx_ds_3 = dNd_s_f3 @ v[:, 0]
            dy_ds_3 = dNd_s_f3 @ v[:, 1]
            nx, ny = -dy_ds_3, dx_ds_3 
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 3, :] = [nx/J_face, ny/J_face, J_face]

    def _map_nodes_to_physical_space(self):
        """
        Maps the reference GLL nodes to physical coordinates using bilinear interpolation.
        """
        r = self.basis.nodes_2d.numpy()[:, 0]; s = self.basis.nodes_2d.numpy()[:, 1]
        N1 = 0.25*(1-r)*(1-s); N2 = 0.25*(1+r)*(1-s); N3 = 0.25*(1+r)*(1+s); N4 = 0.25*(1-r)*(1+s)
        N = np.vstack([N1, N2, N3, N4])
        for i in range(self.mesh.num_elements):
            v = self.mesh.vertices_host[i, :, :]
            self.x_host[i, :] = N.T @ v[:, 0]
            self.y_host[i, :] = N.T @ v[:, 1]

    def _set_initial_conditions(self, Q_host, func):
        """
        Evaluates the initial condition function at all nodes.
        """
        for i in range(self.mesh.num_elements):
            for j in range(self.basis.Np):
                rho, u, v, p = func(self.x_host[i, j], self.y_host[i, j])
                Q_host[i, j, :] = eq.primitive_to_conservative([rho, u, v, p])
