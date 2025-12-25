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
            self.params.rho_floor = 1.0e-5
            self.params.p_floor = 1.0e-5
            self.params.half = 0.5
            self.params.one = 1.0
            # Freestream
            self.params.rho_inf = config.physics.rho_inf
            self.params.u_inf = config.physics.u_inf
            self.params.v_inf = config.physics.v_inf
            self.params.p_inf = config.physics.p_inf
        else:
            self.dtype_np = np.float32
            self.dtype_warp = wp.float32
            self.dtype_vec4 = wp.vec4

            self.params = EquationParams32()
            self.params.gamma = config.physics.gamma
            self.params.rho_floor = 1.0e-5
            self.params.p_floor = 1.0e-5
            self.params.half = 0.5
            self.params.one = 1.0
            # Freestream
            self.params.rho_inf = config.physics.rho_inf
            self.params.u_inf = config.physics.u_inf
            self.params.v_inf = config.physics.v_inf
            self.params.p_inf = config.physics.p_inf

        # Compute Geometric Factors (Metrics & Jacobians) & Physical Coordinates
        # This now resides in the Mesh class to prevent duplication.
        self.mesh.compute_geometry(self.basis, dtype_warp=self.dtype_warp, dtype_np=self.dtype_np)
        
        # Setup Boundary Conditions mask
        bc_manager = BoundaryConditionManager(self.mesh, config)
        self.bc_mask_host = bc_manager.setup_boundary_conditions()

        # --- Device-side data (for computation) ---
        
        self.bc_mask = wp.array(self.bc_mask_host, dtype=wp.int32, device=self.device)
        
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
        Initializes the state vector Q using the provided initial condition function.

        Args:
            initial_condition_func (callable): Function f(x, y) -> (rho, u, v, p).
        """
        # Create SimulationState
        shape = (self.mesh.num_elements, self.basis.Np)
        self.state = SimulationState(
            shape, 
            self.dtype_vec4, 
            self.device, 
            use_filtering=self.cfg.numerics.use_filtering
        )
        
        Q_host = np.zeros((self.mesh.num_elements, self.basis.Np, 4), dtype=self.dtype_np)
        self._set_initial_conditions(Q_host, initial_condition_func)
        
        # Transfer to device state
        self.state.q = wp.array(Q_host, dtype=self.dtype_vec4, device=self.device)

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
                self.mesh.face_geo_factors, 
                self.mesh.J,
                self.bc_mask,
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

    def _set_initial_conditions(self, Q_host, func):
        """
        Evaluates the initial condition function at all nodes.
        """
        for i in range(self.mesh.num_elements):
            for j in range(self.basis.Np):
                rho, u, v, p = func(self.mesh.x_host[i, j], self.mesh.y_host[i, j])
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

        if self.cfg.numerics.limiter == LimiterType.BARTH_JESPERSEN:
            # 2. Compute Neighbor Min/Max
            wp.launch(
                kernel=wk.compute_neighbor_min_max,
                dim=self.mesh.num_elements,
                inputs=[self.state.q_avg, self.state.q_min, self.state.q_max, self.mesh.connectivity, self.params],
                device=self.device
            )
            # 3. Apply Barth-Jespersen
            wp.launch(
                kernel=wk.apply_barth_jespersen_limiter,
                dim=self.mesh.num_elements,
                inputs=[self.state.q, self.state.q_avg, self.state.q_min, self.state.q_max, self.basis.Np, self.params],
                device=self.device
            )
        elif self.cfg.numerics.limiter == LimiterType.MINMOD:
            # Apply Minmod
            wp.launch(
                kernel=wk.apply_minmod_limiter,
                dim=self.mesh.num_elements,
                inputs=[self.state.q, self.state.q_avg, self.mesh.connectivity, self.basis.Np, self.params],
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
