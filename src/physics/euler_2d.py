from .base_solver import BaseSolver
from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.core.boundary_condition_manager import BoundaryConditionManager
from src.physics import equations as eq
from src.kernels import euler_kernels as wk # Renamed import
from src.kernels import common_kernels as ck
from src.kernels.structs import EquationParams32, EquationParams64
from src.core.state import SimulationState
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
        state (SimulationState): The simulation state container.
        face_geo_factors (wp.array): Geometric factors for faces on device.
        bc_mask (wp.array): Boundary condition masks on device.
        x (wp.array): Physical x-coordinates on device.
        y (wp.array): Physical y-coordinates on device.
        M_inv_diag (wp.array): Inverse of the diagonal mass matrix on device.
        max_wave_speed (wp.array): Buffer for maximum wave speed reduction.
        cfg (PazuzuConfig): Typed configuration object.
        params (EquationParams32/64): Struct containing physics constants.
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
            
            self.params = EquationParams64()
            self.params.gamma = config.physics.gamma
            self.params.rho_floor = config.physics.rho_floor
            self.params.p_floor = config.physics.p_floor
            self.params.half = 0.5
            self.params.one = 1.0
            self.params.mu = 0.0
            self.params.prandtl = 0.0
            self.params.cp = 0.0
            self.params.gas_constant = config.physics.gas_constant
            # Freestream
            self.params.rho_inf = config.physics.rho_inf
            self.params.u_inf = config.physics.u_inf
            self.params.v_inf = config.physics.v_inf
            self.params.p_inf = config.physics.p_inf
            self.params.epsilon = config.numerics.hllc_epsilon
        else:
            self.dtype_np = np.float32
            self.dtype_warp = wp.float32
            self.dtype_vec4 = wp.vec4

            self.params = EquationParams32()
            self.params.gamma = config.physics.gamma
            self.params.rho_floor = config.physics.rho_floor
            self.params.p_floor = config.physics.p_floor
            self.params.half = 0.5
            self.params.one = 1.0
            self.params.mu = 0.0
            self.params.prandtl = 0.0
            self.params.cp = 0.0
            self.params.gas_constant = config.physics.gas_constant
            # Freestream
            self.params.rho_inf = config.physics.rho_inf
            self.params.u_inf = config.physics.u_inf
            self.params.v_inf = config.physics.v_inf
            self.params.p_inf = config.physics.p_inf
            self.params.epsilon = config.numerics.hllc_epsilon

        # Compute Geometric Factors (Metrics & Jacobians) & Physical Coordinates
        # This now resides in the Mesh class to prevent duplication.
        self.mesh.compute_geometry(self.basis, dtype_warp=self.dtype_warp, dtype_np=self.dtype_np)
        
        # Setup Boundary Conditions mask and data
        bc_manager = BoundaryConditionManager(self.mesh, config)
        self.bc_mask_host, bc_data_list = bc_manager.setup_boundary_conditions()

        # --- Device-side data (for computation) ---
        
        self.bc_mask = wp.array(self.bc_mask_host, dtype=wp.int32, device=self.device)
        
        # Create and populate BC data array
        num_bcs = len(bc_data_list)
        if num_bcs > 0:
            from src.kernels.structs import BoundaryState32, BoundaryState64
            from src.kernels import boundary_conditions as bc_module
            bc_struct_type = BoundaryState64 if config.numerics.precision == "double" else BoundaryState32
            
            # Use Warp's internal numpy_dtype to ensure correct alignment/padding
            bc_data_host = np.zeros(num_bcs, dtype=bc_struct_type.numpy_dtype())
            
            for i, data in enumerate(bc_data_list):
                bc_id = data['type']
                params = data['params']
                bc_data_host[i]['type'] = bc_id
                
                if bc_id == bc_module.BC_INLET:
                    bc_data_host[i]['v0'] = params['rho']
                    bc_data_host[i]['v1'] = params['u']
                    bc_data_host[i]['v2'] = params['v']
                    bc_data_host[i]['v3'] = params['p']
                elif bc_id == bc_module.BC_OUTLET:
                    bc_data_host[i]['v0'] = params['p_back']
                elif bc_id == bc_module.BC_FARFIELD:
                    bc_data_host[i]['v0'] = params['rho']
                    bc_data_host[i]['v1'] = params['u']
                    bc_data_host[i]['v2'] = params['v']
                    bc_data_host[i]['v3'] = params['p']
            
            self.bc_data = wp.array(bc_data_host, dtype=bc_struct_type, device=self.device)
        else:
            self.bc_data = None

        # Inverse of Mass Matrix (diagonal)
        weights_2d = np.kron(self.basis.weights_1d.numpy(), self.basis.weights_1d.numpy())
        M_inv_diag_host = 1.0 / weights_2d
        self.M_inv_diag = wp.array(M_inv_diag_host, dtype=self.dtype_warp, device=self.device)
        self.weights_2d = wp.array(weights_2d, dtype=self.dtype_warp, device=self.device)
        
        # Buffer for global max wave speed reduction
        self.max_wave_speed = wp.zeros(1, dtype=self.dtype_warp, device=self.device)

        # Ramping parameter for inlet BC
        self.ramp_time = self.dtype_warp(config.simulation.ramp_time)
        
        # Filter Setup
        if config.numerics.use_filtering:
            self.basis.compute_filter_matrix(config.numerics.filter_alpha, config.numerics.filter_order)
            
    def initialize(self, initial_condition_func):
        """
        Initializes the state vector Q using the provided initial condition function and configuration.

        Args:
            initial_condition_func (callable): Function f(x, y, **params) -> (rho, u, v, p).
        """
        # Create SimulationState
        shape = (self.mesh.num_elements, self.basis.Np)
        self.state = SimulationState(
            shape, 
            self.dtype_vec4, 
            self.device, 
            use_filtering=self.cfg.numerics.use_filtering,
            basis=self.basis
        )
        
        Q_host = np.zeros((self.mesh.num_elements, self.basis.Np, 4), dtype=self.dtype_np)
        
        # Extract IC parameters from config
        ic_params = {}
        if not isinstance(self.cfg.initial_condition, str):
            ic_params = self.cfg.initial_condition.params
            
        self._set_initial_conditions(Q_host, initial_condition_func, ic_params)
        
        # Transfer to device state
        self.state.q = wp.array(Q_host, dtype=self.dtype_vec4, device=self.device)

    def compute_rhs(self, t, dt, q, rhs):
        """
        Computes the Right-Hand Side (RHS) of the semi-discrete Euler equations.
        """
        rhs.zero_()
        t_val = self.dtype_warp(t)
        
        # --- Volume Integral ---
        if hasattr(self.basis, 'Nq') and self.basis.Nq > 0:
            # Over-integration (Quadrature Projection)
            wp.launch(
                kernel=wk.interpolate_to_quadrature,
                dim=(self.mesh.num_elements, self.basis.Nq),
                inputs=[q, self.state.q_q, self.basis.Interp_q, self.basis.Np],
                device=self.device
            )
            wp.launch(
                kernel=wk.compute_projected_fluxes,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[self.state.q_q, self.state.f_x_n, self.state.f_y_n, self.basis.Proj_q, self.basis.Nq, self.params],
                device=self.device
            )
            wp.launch(
                kernel=wk.compute_volume_term,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[
                    self.state.f_x_n, self.state.f_y_n,
                    rhs,
                    self.basis.Dr, self.basis.Ds,
                    self.mesh.rx, self.mesh.ry,
                    self.mesh.sx, self.mesh.sy,
                    self.basis.Np,
                    self.params
                ],
                device=self.device
            )
        else:
            # Standard Collocation
            # Note: We need temporary flux buffers. SimulationState allocates them if basis is passed.
            # If basis.Nq was 0, we might need fallback buffers or just compute inside a combined kernel.
            # For robustness, we'll ensure f_x_n, f_y_n exist or use a combined kernel.
            # Since SimulationState now always checks for basis, let's assume they are there if needed.
            # Actually, let's add them to SimulationState for standard case too if needed, 
            # or just use a combined kernel to avoid allocation.
            # To keep it simple, I'll update SimulationState to always have these flux buffers.
            wp.launch(
                kernel=wk.compute_nodal_fluxes,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[q, self.state.f_x_n, self.state.f_y_n, self.params],
                device=self.device
            )
            wp.launch(
                kernel=wk.compute_volume_term,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[
                    self.state.f_x_n, self.state.f_y_n,
                    rhs,
                    self.basis.Dr, self.basis.Ds,
                    self.mesh.rx, self.mesh.ry,
                    self.mesh.sx, self.mesh.sy,
                    self.basis.Np,
                    self.params
                ],
                device=self.device
            )
        
        # Determine Flux Type ID
        flux_type_id = wk.FLUX_RUSANOV
        if self.cfg.numerics.flux_type == FluxType.HLLC:
            flux_type_id = wk.FLUX_HLLC
        
        # --- Surface Integral ---
        # Adds LIFT * (NumericalFlux - NormalFlux)
        # Parallelize over (Elements, Faces, FaceNodes)
        wp.launch(
            kernel=wk.compute_surface_term,
            dim=(self.mesh.num_elements, 4, self.basis.Nfp),
            inputs=[
                q,
                self.state.f_x_n,
                self.state.f_y_n,
                rhs,
                self.mesh.connectivity,
                self.mesh.connectivity_face_indices,
                self.basis.face_nodes,
                self.basis.LIFT,
                self.mesh.face_geo_factors, 
                self.mesh.J,
                self.bc_mask,
                self.bc_data, # Pass validated BC data
                self.mesh.x, # Physical X coordinates
                self.mesh.y, # Physical Y coordinates
                self.basis.Nfp,
                t_val,
                self.ramp_time,
                flux_type_id,
                self.params
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
                self.state.q, 
                self.max_wave_speed,
                self.params
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

    def _set_initial_conditions(self, Q_host, func, params):
        """
        Evaluates the initial condition function at all nodes.
        """
        for i in range(self.mesh.num_elements):
            for j in range(self.basis.Np):
                rho, u, v, p = func(self.mesh.x_host[i, j], self.mesh.y_host[i, j], **params)
                Q_host[i, j, :] = eq.primitive_to_conservative([rho, u, v, p])

    def filter_solution(self):
        """
        Applies exponential filtering (spectral viscosity) to the solution state.
        This suppresses aliasing errors and stabilizes high-order simulations.
        """
        if self.basis.filter_matrix is not None and self.state.filter_buffer is not None:
            # Apply filter: Q -> FilterBuffer
            wp.launch(
                kernel=ck.apply_filter_matrix,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[self.state.q, self.basis.filter_matrix, self.state.filter_buffer],
                device=self.device
            )
            # Copy back: FilterBuffer -> Q
            wp.copy(self.state.q, self.state.filter_buffer)

    def apply_limiter(self):
        """
        Applies the configured limiter to the solution state.
        Used for shock capturing and ensuring physical bounds (positivity).
        """
        if self.cfg.numerics.limiter == LimiterType.NONE:
            return

        # 1. Compute Cell Averages
        wp.launch(
            kernel=wk.compute_cell_averages,
            dim=self.mesh.num_elements,
            inputs=[self.state.q, self.state.q_avg, self.weights_2d, self.mesh.J, self.basis.Np, self.params],
            device=self.device
        )

        # 2. Compute Neighbor Min/Max (required for both limiters)
        wp.launch(
            kernel=wk.compute_neighbor_min_max,
            dim=self.mesh.num_elements,
            inputs=[self.state.q_avg, self.state.q_min, self.state.q_max, self.mesh.connectivity, self.params],
            device=self.device
        )

        if self.cfg.numerics.limiter == LimiterType.BARTH_JESPERSEN:
            # 3. Apply Barth-Jespersen
            wp.launch(
                kernel=wk.apply_barth_jespersen_limiter,
                dim=self.mesh.num_elements,
                inputs=[self.state.q, self.state.q_avg, self.state.q_min, self.state.q_max, self.basis.Np, self.params],
                device=self.device
            )
        elif self.cfg.numerics.limiter == LimiterType.MINMOD:
            # 2. Compute Gradients using Green-Gauss
            wp.launch(
                kernel=wk.compute_gradients_green_gauss,
                dim=self.mesh.num_elements,
                inputs=[
                    self.state.q_avg, 
                    self.mesh.connectivity, 
                    self.mesh.face_geo_factors, 
                    self.mesh.vol,
                    self.state.grad_x,
                    self.state.grad_y,
                    self.params
                ],
                device=self.device
            )

            # 3. Apply Gradient-based Minmod
            wp.launch(
                kernel=wk.apply_minmod_limiter,
                dim=self.mesh.num_elements,
                inputs=[
                    self.state.q, 
                    self.state.q_avg, 
                    self.state.q_min, 
                    self.state.q_max,
                    self.state.grad_x,
                    self.state.grad_y,
                    self.mesh.centroid,
                    self.mesh.x,
                    self.mesh.y,
                    self.basis.Np, 
                    self.params
                ],
                device=self.device
            )

    def post_step(self):
        """
        Hook called by the driver after each full time step.
        """
        if self.cfg.numerics.use_filtering:
            self.filter_solution()
            
        if self.cfg.numerics.limiter != LimiterType.NONE:
            self.apply_limiter()