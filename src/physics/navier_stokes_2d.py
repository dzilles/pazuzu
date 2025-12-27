from .base_solver import BaseSolver
from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.core.boundary_condition_manager import BoundaryConditionManager
from src.physics import equations as eq
from src.kernels import euler_kernels as ek
from src.kernels import navier_stokes_kernels as nsk
from src.kernels import common_kernels as ck
from src.kernels.structs import EquationParams32, EquationParams64
from src.core.state import SimulationState
import numpy as np
import warp as wp
from src.core.config import PazuzuConfig, FluxType, LimiterType

class NavierStokes2DSolver(BaseSolver):
    """
    Solves the 2D Navier-Stokes equations using a Nodal Discontinuous Galerkin method.
    Uses the BR1 (Bassi-Rebay 1) scheme for viscous terms.
    """
    def __init__(self, mesh, basis, config: PazuzuConfig):
        super().__init__(mesh, basis, config)
        self.cfg = config

        # --- Precision Setup ---
        if config.numerics.precision == "double":
            self.dtype_np = np.float64
            self.dtype_warp = wp.float64
            self.dtype_vec4 = wp.vec4d
            self.dtype_vec2 = wp.vec2d
            
            self.params = EquationParams64()
            self.params.gamma = config.physics.gamma
            self.params.mu = config.physics.mu
            self.params.prandtl = config.physics.prandtl
            self.params.cp = config.physics.cp
            self.params.gas_constant = config.physics.gas_constant
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
            self.dtype_vec2 = wp.vec2

            self.params = EquationParams32()
            self.params.gamma = config.physics.gamma
            self.params.mu = config.physics.mu
            self.params.prandtl = config.physics.prandtl
            self.params.cp = config.physics.cp
            self.params.gas_constant = config.physics.gas_constant
            self.params.rho_floor = 1.0e-5
            self.params.p_floor = 1.0e-5
            self.params.half = 0.5
            self.params.one = 1.0
            # Freestream
            self.params.rho_inf = config.physics.rho_inf
            self.params.u_inf = config.physics.u_inf
            self.params.v_inf = config.physics.v_inf
            self.params.p_inf = config.physics.p_inf

        self.mesh.compute_geometry(self.basis, dtype_warp=self.dtype_warp, dtype_np=self.dtype_np)
        
        bc_manager = BoundaryConditionManager(self.mesh, config)
        self.bc_mask_host, bc_data_list = bc_manager.setup_boundary_conditions()
        self.bc_mask = wp.array(self.bc_mask_host, dtype=wp.int32, device=self.device)
        
        num_bcs = len(bc_data_list)
        if num_bcs > 0:
            from src.kernels.structs import BoundaryState32, BoundaryState64
            from src.kernels import boundary_conditions as bc_module
            bc_struct_type = BoundaryState64 if config.numerics.precision == "double" else BoundaryState32
            bc_data_host = np.zeros(num_bcs, dtype=bc_struct_type.numpy_dtype())
            
            for i, data in enumerate(bc_data_list):
                bc_id = data['type']
                params = data['params']
                bc_data_host[i]['type'] = bc_id
                
                if bc_id == bc_module.BC_INLET or bc_id == bc_module.BC_FARFIELD:
                    bc_data_host[i]['v0'] = params['rho']
                    bc_data_host[i]['v1'] = params['u']
                    bc_data_host[i]['v2'] = params['v']
                    bc_data_host[i]['v3'] = params['p']
                elif bc_id == bc_module.BC_OUTLET:
                    bc_data_host[i]['v0'] = params['p_back']
                elif bc_id == bc_module.BC_ISOTHERMAL_WALL:
                    bc_data_host[i]['v0'] = params['T_wall']
            
            self.bc_data = wp.array(bc_data_host, dtype=bc_struct_type, device=self.device)
        else:
            self.bc_data = None

        weights_2d = np.kron(self.basis.weights_1d.numpy(), self.basis.weights_1d.numpy())
        self.weights_2d = wp.array(weights_2d, dtype=self.dtype_warp, device=self.device)
        self.max_wave_speed = wp.zeros(1, dtype=self.dtype_warp, device=self.device)
        self.ramp_time = self.dtype_warp(config.simulation.ramp_time)
        
        # --- Gradient Buffers ---
        shape = (self.mesh.num_elements, self.basis.Np)
        self.grad_u = wp.zeros(shape, dtype=self.dtype_vec2, device=self.device)
        self.grad_v = wp.zeros(shape, dtype=self.dtype_vec2, device=self.device)
        self.grad_T = wp.zeros(shape, dtype=self.dtype_vec2, device=self.device)
        self.has_nan = wp.zeros(1, dtype=wp.int32, device=self.device)

        if config.numerics.use_filtering:
            self.basis.compute_filter_matrix(config.numerics.filter_alpha, config.numerics.filter_order)
            
    def initialize(self, initial_condition_func):
        shape = (self.mesh.num_elements, self.basis.Np)
        self.state = SimulationState(
            shape, 
            self.dtype_vec4, 
            self.device, 
            use_filtering=self.cfg.numerics.use_filtering,
            basis=self.basis
        )
        Q_host = np.zeros((self.mesh.num_elements, self.basis.Np, 4), dtype=self.dtype_np)
        ic_params = {}
        if not isinstance(self.cfg.initial_condition, str):
            ic_params = self.cfg.initial_condition.params
        self._set_initial_conditions(Q_host, initial_condition_func, ic_params)
        self.state.q = wp.array(Q_host, dtype=self.dtype_vec4, device=self.device)

    def compute_rhs(self, t, dt, q, rhs):
        rhs.zero_()
        t_val = self.dtype_warp(t)
        
        # --- PASS 1: Primitive Gradients ---
        self.grad_u.zero_()
        self.grad_v.zero_()
        self.grad_T.zero_()
        
        wp.launch(
            kernel=nsk.compute_primitive_gradients_volume,
            dim=(self.mesh.num_elements, self.basis.Np),
            inputs=[q, self.grad_u, self.grad_v, self.grad_T, 
                    self.basis.Dr, self.basis.Ds, self.mesh.rx, self.mesh.ry, self.mesh.sx, self.mesh.sy, 
                    self.basis.Np, self.params],
            device=self.device
        )
        wp.launch(
            kernel=nsk.compute_primitive_gradients_surface,
            dim=self.mesh.num_elements,
            inputs=[q, self.grad_u, self.grad_v, self.grad_T,
                    self.mesh.connectivity, self.mesh.connectivity_face_indices, self.basis.face_nodes, self.basis.LIFT,
                    self.mesh.face_geo_factors, self.mesh.J, self.bc_mask, self.bc_data,
                    self.mesh.x, self.mesh.y, self.basis.Nfp, t_val, self.ramp_time, self.params],
            device=self.device
        )
        
        # --- PASS 2: Inviscid RHS ---
        if hasattr(self.basis, 'Nq') and self.basis.Nq > 0:
            # Over-integration (Quadrature Projection)
            wp.launch(
                kernel=ek.interpolate_to_quadrature,
                dim=(self.mesh.num_elements, self.basis.Nq),
                inputs=[q, self.state.q_q, self.basis.Interp_q, self.basis.Np],
                device=self.device
            )
            wp.launch(
                kernel=ek.compute_projected_fluxes,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[self.state.q_q, self.state.f_x_n, self.state.f_y_n, self.basis.Proj_q, self.basis.Nq, self.params],
                device=self.device
            )
            wp.launch(
                kernel=ek.compute_volume_term,
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
            wp.launch(
                kernel=ek.compute_nodal_fluxes,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[q, self.state.f_x_n, self.state.f_y_n, self.params],
                device=self.device
            )
            wp.launch(
                kernel=ek.compute_volume_term,
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
        flux_type_id = ek.FLUX_RUSANOV
        if self.cfg.numerics.flux_type == FluxType.HLLC:
            flux_type_id = ek.FLUX_HLLC
        wp.launch(
            kernel=ek.compute_surface_term,
            dim=self.mesh.num_elements,
            inputs=[q, self.state.f_x_n, self.state.f_y_n, rhs, self.mesh.connectivity, self.mesh.connectivity_face_indices, self.basis.face_nodes, self.basis.LIFT,
                    self.mesh.face_geo_factors, self.mesh.J, self.bc_mask, self.bc_data,
                    self.mesh.x, self.mesh.y, self.basis.Nfp, t_val, self.ramp_time, flux_type_id, self.params],
            device=self.device
        )
        
        # --- PASS 3: Viscous RHS ---
        wp.launch(
            kernel=nsk.compute_viscous_volume_term,
            dim=(self.mesh.num_elements, self.basis.Np),
            inputs=[q, self.grad_u, self.grad_v, self.grad_T, rhs,
                    self.basis.Dr, self.basis.Ds, self.mesh.rx, self.mesh.ry, self.mesh.sx, self.mesh.sy, self.basis.Np, self.params],
            device=self.device
        )
        
        wp.launch(
            kernel=nsk.compute_viscous_surface_term,
            dim=self.mesh.num_elements,
            inputs=[q, self.grad_u, self.grad_v, self.grad_T, rhs,
                    self.mesh.connectivity, self.mesh.connectivity_face_indices, self.basis.face_nodes, self.basis.LIFT,
                    self.mesh.face_geo_factors, self.mesh.J, self.bc_mask, self.bc_data,
                    self.mesh.x, self.mesh.y, self.basis.Nfp, t_val, self.ramp_time, self.params],
            device=self.device
        )

    def calculate_dt(self, CFL):
        self.max_wave_speed.zero_()
        wp.launch(
            kernel=ek.compute_max_wave_speed,
            dim=(self.mesh.num_elements, self.basis.Np),
            inputs=[self.state.q, self.max_wave_speed, self.params],
            device=self.device
        )
        max_speed = self.max_wave_speed.numpy()[0]
        h_eff = self.mesh.dx * (self.basis.min_node_dist / 2.0)
        rho_min = self.cfg.physics.rho_inf
        if rho_min < 1e-10: rho_min = 1.0
        visc_speed = (2.0 * self.cfg.physics.mu) / (rho_min * h_eff)
        total_speed = max_speed + visc_speed
        if total_speed < 1e-6:
            total_speed = 1.0
        dt = CFL * h_eff / total_speed
        return dt

    def _set_initial_conditions(self, Q_host, func, params):
        for i in range(self.mesh.num_elements):
            for j in range(self.basis.Np):
                rho, u, v, p = func(self.mesh.x_host[i, j], self.mesh.y_host[i, j], **params)
                Q_host[i, j, :] = eq.primitive_to_conservative([rho, u, v, p])

    def post_step(self):
        if self.cfg.numerics.use_filtering:
            if self.basis.filter_matrix is not None and self.state.filter_buffer is not None:
                wp.launch(
                    kernel=ck.apply_filter_matrix,
                    dim=(self.mesh.num_elements, self.basis.Np),
                    inputs=[self.state.q, self.basis.filter_matrix, self.state.filter_buffer],
                    device=self.device
                )
                wp.copy(self.state.q, self.state.filter_buffer)
            
        if self.cfg.numerics.limiter != LimiterType.NONE:
            wp.launch(
                kernel=ek.compute_cell_averages,
                dim=self.mesh.num_elements,
                inputs=[self.state.q, self.state.q_avg, self.weights_2d, self.mesh.J, self.basis.Np, self.params],
                device=self.device
            )
            if self.cfg.numerics.limiter == LimiterType.BARTH_JESPERSEN:
                wp.launch(
                    kernel=ek.compute_neighbor_min_max,
                    dim=self.mesh.num_elements,
                    inputs=[self.state.q_avg, self.state.q_min, self.state.q_max, self.mesh.connectivity, self.params],
                    device=self.device
                )
                wp.launch(
                    kernel=ek.apply_barth_jespersen_limiter,
                    dim=self.mesh.num_elements,
                    inputs=[self.state.q, self.state.q_avg, self.state.q_min, self.state.q_max, self.basis.Np, self.params],
                    device=self.device
                )
            elif self.cfg.numerics.limiter == LimiterType.MINMOD:
                wp.launch(
                    kernel=ek.apply_minmod_limiter,
                    dim=self.mesh.num_elements,
                    inputs=[self.state.q, self.state.q_avg, self.mesh.connectivity, self.basis.Np, self.params],
                    device=self.device
                )