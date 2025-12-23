from .base_solver import BaseSolver
from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.physics import equations as eq
from src.kernels import euler_kernels as wk # Renamed import
from src.kernels import common_kernels as ck
import numpy as np
import warp as wp
from src.core.config import PazuzuConfig, FluxType, LimiterType

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
        cfg (PazuzuConfig): Typed configuration object.
    """
    def __init__(self, mesh, basis, config: PazuzuConfig):
        """
        Initializes the Euler 2D Solver.

        Args:
            mesh (Mesh): The computational mesh.
            basis (Basis): The DG basis.
            config (PazuzuConfig): Configuration object.
        """
        super().__init__(mesh, basis, config)
        
        self.cfg = config

        # --- Precision Setup ---
        if config.numerics.precision == "double":
            self.dtype_np = np.float64
            self.dtype_warp = wp.float64
            self.dtype_vec4 = wp.vec4d
        else:
            self.dtype_np = np.float32
            self.dtype_warp = wp.float32
            self.dtype_vec4 = wp.vec4

        # --- Host-side data (for setup) ---
        # Volume metrics: dx/dr, dx/ds, dy/dr, dy/ds, J
        self.geo_factors_host = np.zeros((mesh.num_elements, basis.Np, 5), dtype=self.dtype_np)
        # Face metrics: nx, ny, J_surf
        self.face_geo_factors_host = np.zeros((mesh.num_elements, 4, 3), dtype=self.dtype_np)
        # Physical coordinates
        self.x_host = np.zeros((mesh.num_elements, basis.Np), dtype=self.dtype_np)
        self.y_host = np.zeros((mesh.num_elements, basis.Np), dtype=self.dtype_np)
        
        # Calculate geometric factors (Metrics & Jacobians)
        self._calculate_geometric_factors()
        
        # Map reference nodes to physical space
        self._map_nodes_to_physical_space()

        # Setup Boundary Conditions mask
        self._setup_boundary_conditions()

        # --- Device-side data (for computation) ---
        
        # Update Mesh Metric Arrays on Device (ensure they are synced)
        self.mesh.rx = wp.array(self.mesh.rx_host, dtype=self.dtype_warp, device=self.device)
        self.mesh.ry = wp.array(self.mesh.ry_host, dtype=self.dtype_warp, device=self.device)
        self.mesh.sx = wp.array(self.mesh.sx_host, dtype=self.dtype_warp, device=self.device)
        self.mesh.sy = wp.array(self.mesh.sy_host, dtype=self.dtype_warp, device=self.device)
        self.mesh.J  = wp.array(self.mesh.J_host, dtype=self.dtype_warp, device=self.device)
        
        self.face_geo_factors = wp.array(self.face_geo_factors_host, dtype=self.dtype_warp, device=self.device)
        self.bc_mask = wp.array(self.bc_mask_host, dtype=wp.int32, device=self.device)
        
        self.x = wp.array(self.x_host, dtype=self.dtype_warp, device=self.device)
        self.y = wp.array(self.y_host, dtype=self.dtype_warp, device=self.device)

        # Inverse of Mass Matrix (diagonal)
        weights_2d = np.kron(self.basis.weights_1d.numpy(), self.basis.weights_1d.numpy())
        M_inv_diag_host = 1.0 / weights_2d
        self.M_inv_diag = wp.array(M_inv_diag_host, dtype=self.dtype_warp, device=self.device)
        
        # Buffer for global max wave speed reduction
        self.max_wave_speed = wp.zeros(1, dtype=self.dtype_warp, device=self.device)

        # Ramping parameter for inlet BC
        self.ramp_time = self.dtype_warp(config.simulation.ramp_time)
        
        # Filter Setup
        self.filter_buffer = None
        if config.numerics.use_filtering:
            self.basis.compute_filter_matrix(config.numerics.filter_alpha, config.numerics.filter_order)
            
        # Physics Constants (Typed)
        self.gamma_val = self.dtype_warp(config.physics.gamma)
        self.rho_floor_val = self.dtype_warp(1.0e-5)
        self.p_floor_val = self.dtype_warp(1.0e-5)
        
        self.half_val = self.dtype_warp(0.5)
        self.one_val = self.dtype_warp(1.0)
            
    def initialize(self, initial_condition_func):

        """
        Initializes the state vector Q using the provided initial condition function.

        Args:
            initial_condition_func (callable): Function f(x, y) -> (rho, u, v, p).
        """
        Q_host = np.zeros((self.mesh.num_elements, self.basis.Np, 4), dtype=self.dtype_np)
        self._set_initial_conditions(Q_host, initial_condition_func)
        self.Q = wp.array(Q_host, dtype=self.dtype_vec4, device=self.device)
        self.rhs = wp.zeros_like(self.Q)
        
        if self.cfg.numerics.use_filtering:
            self.filter_buffer = wp.zeros_like(self.Q)

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
        t_val = self.dtype_warp(t)
        
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
                self.basis.Np,
                self.gamma_val,
                self.rho_floor_val,
                self.p_floor_val,
                self.half_val,
                self.one_val
            ],
            device=self.device
        )
        
        # Determine Flux Type ID
        flux_type_id = wk.FLUX_RUSANOV
        if self.cfg.numerics.flux_type == FluxType.HLLC:
            flux_type_id = wk.FLUX_HLLC
        
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
                self.x, # Physical X coordinates
                self.y, # Physical Y coordinates
                self.basis.Nfp,
                t_val,
                self.ramp_time,
                flux_type_id,
                self.gamma_val,
                self.rho_floor_val,
                self.p_floor_val,
                self.half_val,
                self.one_val
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
            inputs=[
                self.Q, 
                self.max_wave_speed,
                self.gamma_val,
                self.rho_floor_val,
                self.p_floor_val,
                self.half_val,
                self.one_val
            ],
            device=self.device
        )
        max_speed = self.max_wave_speed.numpy()[0]
        
        # Safety check to avoid division by zero if fluid is at rest and T=0 (unlikely)
        if max_speed < 1.0e-6:
            max_speed = 1.0 

        # Characteristic length / Max Speed
        # Effective grid size for stability (scales with 1/N^2)
        h_eff = self.mesh.dx * (self.basis.min_node_dist / 2.0)
        dt = CFL * h_eff / max_speed
        return dt

    def _setup_boundary_conditions(self):
        """
        Parses the configuration to setup the boundary condition mask.
        Maps physical tags from the mesh to solver-specific BC IDs.
        """
        self.bc_mask_host = np.zeros((self.mesh.num_elements, 4), dtype=np.int32)
        
        if self.config.boundaries:
            BC_WALL = 1
            BC_FARFIELD = 2      # Freestream / Characteristic
            BC_INLET = 3
            BC_OUTLET = 4
            BC_EXTRAPOLATION = 5 # 0th Order Extrapolation
            BC_CYLINDER_WALL = 6 # Slip Wall with analytical normals
            
            tag_to_bc = {}
            for name, bc_conf in self.config.boundaries.items():
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
                    elif bc_type == "cylinder_wall": tag_to_bc[tag] = BC_CYLINDER_WALL
                    elif bc_type == "farfield": tag_to_bc[tag] = BC_FARFIELD
                    elif bc_type in ["outflow", "extrapolation"]: tag_to_bc[tag] = BC_EXTRAPOLATION
                    elif bc_type == "inlet": tag_to_bc[tag] = BC_INLET
                    elif bc_type == "outlet": tag_to_bc[tag] = BC_OUTLET
            
            # Apply to mask
            for e in range(self.mesh.num_elements):
                for f in range(4):
                    tag = self.mesh.boundary_tags_host[e, f]
                    if tag > 0:
                        # Default to Farfield (Freestream) if tag exists but type not specified
                        self.bc_mask_host[e, f] = tag_to_bc.get(tag, BC_FARFIELD)

    def _calculate_geometric_factors(self):
        """
        Calculates geometric metrics (Jacobians) for all elements.
        Computes metrics at every node to support general unstructured quadrilaterals.
        """
        # Get reference nodes (r, s) for the basis
        # Shape: (Np, 2)
        nodes_2d = self.basis.nodes_2d.numpy()
        r = nodes_2d[:, 0]
        s = nodes_2d[:, 1]
        
        # Derivatives of bilinear shape functions w.r.t r and s at all Np nodes
        # Shape functions:
        # N0 = 0.25(1-r)(1-s)
        # N1 = 0.25(1+r)(1-s)
        # N2 = 0.25(1+r)(1+s)
        # N3 = 0.25(1-r)(1+s)
        
        # dN/dr
        # dN0_dr = -0.25(1-s)
        # dN1_dr =  0.25(1-s)
        # dN2_dr =  0.25(1+s)
        # dN3_dr = -0.25(1+s)
        
        # dN/ds
        # dN0_ds = -0.25(1-r)
        # dN1_ds = -0.25(1+r)
        # dN2_ds =  0.25(1+r)
        # dN3_ds =  0.25(1-r)
        
        # Precompute shape function derivatives at all nodes
        # Shape: (Np, 4)
        dN_dr = np.zeros((self.basis.Np, 4))
        dN_ds = np.zeros((self.basis.Np, 4))
        
        dN_dr[:, 0] = -0.25 * (1 - s)
        dN_dr[:, 1] =  0.25 * (1 - s)
        dN_dr[:, 2] =  0.25 * (1 + s)
        dN_dr[:, 3] = -0.25 * (1 + s)
        
        dN_ds[:, 0] = -0.25 * (1 - r)
        dN_ds[:, 1] = -0.25 * (1 + r)
        dN_ds[:, 2] =  0.25 * (1 + r)
        dN_ds[:, 3] =  0.25 * (1 - r)

        # Derivatives for face normals (still constant per face for bilinear elements)
        # Face 0 (Bottom, s=-1): dN/dr at r=0, s=-1 -> [-0.5, 0.5, 0.0, 0.0]
        dNd_r_f0 = np.array([-0.5, 0.5, 0.0, 0.0])
        # Face 1 (Right, r=1): dN/ds at r=1, s=0 -> [0.0, -0.5, 0.5, 0.0]
        dNd_s_f1 = np.array([0.0, -0.5, 0.5, 0.0])
        # Face 2 (Top, s=1): dN/dr at r=0, s=1 -> [0.0, 0.0, 0.5, -0.5]
        dNd_r_f2 = np.array([0.0, 0.0, 0.5, -0.5])
        # Face 3 (Left, r=-1): dN/ds at r=-1, s=0 -> [-0.5, 0.0, 0.0, 0.5]
        dNd_s_f3 = np.array([-0.5, 0.0, 0.0, 0.5])
        
        # Allocate/Reset Mesh Metric Arrays (Host)
        # These now need to be (NumElements, Np)
        self.mesh.rx_host = np.zeros((self.mesh.num_elements, self.basis.Np), dtype=self.dtype_np)
        self.mesh.ry_host = np.zeros((self.mesh.num_elements, self.basis.Np), dtype=self.dtype_np)
        self.mesh.sx_host = np.zeros((self.mesh.num_elements, self.basis.Np), dtype=self.dtype_np)
        self.mesh.sy_host = np.zeros((self.mesh.num_elements, self.basis.Np), dtype=self.dtype_np)
        self.mesh.J_host  = np.zeros((self.mesh.num_elements, self.basis.Np), dtype=self.dtype_np)

        for i in range(self.mesh.num_elements):
            v = self.mesh.vertices_host[i, :, :] 
            
            # Compute Jacobians at all Np nodes
            # x = sum(N_k * x_k) -> dx/dr = sum(dN_k/dr * x_k)
            # v[:, 0] is x-coords of vertices (4,)
            # dN_dr is (Np, 4)
            # dx_dr becomes (Np,)
            dx_dr = dN_dr @ v[:, 0]
            dx_ds = dN_ds @ v[:, 0]
            dy_dr = dN_dr @ v[:, 1]
            dy_ds = dN_ds @ v[:, 1]
            
            # Jacobian Determinant J = dx/dr * dy/ds - dx/ds * dy/dr
            J = dx_dr * dy_ds - dx_ds * dy_dr
            self.mesh.J_host[i, :] = J
            
            # Inverse Jacobian components
            # dr/dx =  dy/ds / J
            # dr/dy = -dx/ds / J
            # ds/dx = -dy/dr / J
            # ds/dy =  dx/dr / J
            
            # Avoid division by zero (though J should be positive for valid meshes)
            # J_inv = 1.0 / J
            
            self.mesh.rx_host[i, :] =  dy_ds / J
            self.mesh.ry_host[i, :] = -dx_ds / J
            self.mesh.sx_host[i, :] = -dy_dr / J
            self.mesh.sy_host[i, :] =  dx_dr / J

            self.geo_factors_host[i, :, 0] = dx_dr
            self.geo_factors_host[i, :, 1] = dx_ds
            self.geo_factors_host[i, :, 2] = dy_dr
            self.geo_factors_host[i, :, 3] = dy_ds
            self.geo_factors_host[i, :, 4] = J
            
            # --- Face Normals and Surface Jacobians ---
            # Normals are outward pointing.
            # For bilinear quads with straight edges, these are constant per face.
            
            # Face 0 (bottom)
            dx_dr_0 = dNd_r_f0 @ v[:, 0]
            dy_dr_0 = dNd_r_f0 @ v[:, 1]
            nx, ny = dy_dr_0, -dx_dr_0
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 0, :] = [nx/J_face, ny/J_face, J_face]

            # Face 1 (right)
            dx_ds_1 = dNd_s_f1 @ v[:, 0]
            dy_ds_1 = dNd_s_f1 @ v[:, 1]
            nx, ny = dy_ds_1, -dx_ds_1
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 1, :] = [nx/J_face, ny/J_face, J_face]
            
            # Face 2 (top)
            dx_dr_2 = dNd_r_f2 @ v[:, 0]
            dy_dr_2 = dNd_r_f2 @ v[:, 1]
            nx, ny = -dy_dr_2, dx_dr_2 
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 2, :] = [nx/J_face, ny/J_face, J_face]

            # Face 3 (left)
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

    def filter_solution(self):
        """
        Applies exponential filtering (spectral viscosity) to the solution state.
        This suppresses aliasing errors and stabilizes high-order simulations.
        """
        if self.basis.filter_matrix is not None and self.filter_buffer is not None:
            # Apply filter: Q -> FilterBuffer
            wp.launch(
                kernel=ck.apply_filter_matrix,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[self.Q, self.basis.filter_matrix, self.filter_buffer],
                device=self.device
            )
            # Copy back: FilterBuffer -> Q
            wp.copy(self.Q, self.filter_buffer)

    def post_step(self):
        """
        Hook called by the driver after each full time step.
        """
        if self.cfg.numerics.use_filtering:
            self.filter_solution()
