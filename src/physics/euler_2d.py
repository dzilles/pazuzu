from .base_solver import BaseSolver
from src.geometry.mesh import Mesh
from src.core.basis import Basis
from src.physics import equations as eq
from src.kernels import euler_kernels as wk # Renamed import
import numpy as np
import warp as wp

class Euler2DSolver(BaseSolver):
    def __init__(self, mesh, basis, config=None):
        super().__init__(mesh, basis, config)
        
        # --- Host-side data (for setup) ---
        self.geo_factors_host = np.zeros((mesh.num_elements, 5), dtype=np.float32) # dx_dr, dx_ds, dy_dr, dy_ds, J
        self.face_geo_factors_host = np.zeros((mesh.num_elements, 4, 3), dtype=np.float32) # nx, ny, J_face
        self.x_host = np.zeros((mesh.num_elements, basis.Np), dtype=np.float32)
        self.y_host = np.zeros((mesh.num_elements, basis.Np), dtype=np.float32)
        
        # Calculate geometric factors (Metrics & Jacobians)
        self._calculate_geometric_factors()
        
        # Map nodes to physical space
        self._map_nodes_to_physical_space()

        # BC Setup
        self._setup_boundary_conditions()

        # --- Device-side data (for computation) ---
        
        # Update Mesh Metric Arrays on Device
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
        
        self.max_wave_speed = wp.zeros(1, dtype=wp.float32, device=self.device)

    def initialize(self, initial_condition_func):
        Q_host = np.zeros((self.mesh.num_elements, self.basis.Np, 4), dtype=np.float32)
        self._set_initial_conditions(Q_host, initial_condition_func)
        self.Q = wp.array(Q_host, dtype=wp.vec4, device=self.device)
        self.rhs = wp.zeros_like(self.Q)

    def compute_rhs(self, t, dt, q, rhs):
        """ Computes the right-hand side of the ODE using Warp kernels. """
        rhs.zero_()
        
        # --- Volume Integral ---
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
        wp.launch(
            kernel=wk.compute_surface_term,
            dim=self.mesh.num_elements,
            inputs=[
                q,
                rhs,
                self.mesh.connectivity,
                self.basis.face_nodes,
                self.basis.LIFT,
                self.face_geo_factors, 
                self.mesh.J,
                self.bc_mask,
                self.basis.Nfp
            ],
            device=self.device
        )

    def calculate_dt(self, CFL):
        self.max_wave_speed.zero_()
        wp.launch(
            kernel=wk.compute_max_wave_speed,
            dim=(self.mesh.num_elements, self.basis.Np),
            inputs=[self.Q, self.max_wave_speed],
            device=self.device
        )
        max_speed = self.max_wave_speed.numpy()[0]
        
        if max_speed < 1.0e-6:
            max_speed = 1.0 

        dt = CFL * self.mesh.dx / max_speed
        return dt

    def _setup_boundary_conditions(self):
        self.bc_mask_host = np.zeros((self.mesh.num_elements, 4), dtype=np.int32)
        
        if self.config and 'boundaries' in self.config:
            BC_WALL = 1
            BC_FARFIELD = 2
            
            tag_to_bc = {}
            for name, bc_conf in self.config['boundaries'].items():
                tag = -1
                if hasattr(self.mesh, 'physical_groups') and name in self.mesh.physical_groups:
                    tag = self.mesh.physical_groups[name]
                else:
                    try: tag = int(name) 
                    except: pass
                
                if tag != -1:
                    bc_type = bc_conf.get('type')
                    if bc_type == "slip_wall": tag_to_bc[tag] = BC_WALL
                    elif bc_type in ["farfield", "outflow"]: tag_to_bc[tag] = BC_FARFIELD
            
            for e in range(self.mesh.num_elements):
                for f in range(4):
                    tag = self.mesh.boundary_tags_host[e, f]
                    if tag > 0:
                        self.bc_mask_host[e, f] = tag_to_bc.get(tag, BC_FARFIELD)

    def _calculate_geometric_factors(self):
        dNd_r_center = 0.25 * np.array([-1, 1, 1, -1])
        dNd_s_center = 0.25 * np.array([-1, -1, 1, 1])
        
        for i in range(self.mesh.num_elements):
            v = self.mesh.vertices_host[i, :, :] 
            dx_dr = dNd_r_center @ v[:, 0]
            dx_ds = dNd_s_center @ v[:, 0]
            dy_dr = dNd_r_center @ v[:, 1]
            dy_ds = dNd_s_center @ v[:, 1]
            
            J = dx_dr * dy_ds - dx_ds * dy_dr
            self.mesh.J_host[i] = J
            
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
            
            nx, ny = dy_dr, -dx_dr
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 0, :] = [nx/J_face, ny/J_face, J_face]

            nx, ny = dy_ds, -dx_ds
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 1, :] = [nx/J_face, ny/J_face, J_face]
            
            nx, ny = -dy_dr, dx_dr 
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 2, :] = [nx/J_face, ny/J_face, J_face]

            nx, ny = -dy_ds, dx_ds 
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 3, :] = [nx/J_face, ny/J_face, J_face]

    def _map_nodes_to_physical_space(self):
        r = self.basis.nodes_2d.numpy()[:, 0]; s = self.basis.nodes_2d.numpy()[:, 1]
        N1 = 0.25*(1-r)*(1-s); N2 = 0.25*(1+r)*(1-s); N3 = 0.25*(1+r)*(1+s); N4 = 0.25*(1-r)*(1+s)
        N = np.vstack([N1, N2, N3, N4])
        for i in range(self.mesh.num_elements):
            v = self.mesh.vertices_host[i, :, :]
            self.x_host[i, :] = N.T @ v[:, 0]
            self.y_host[i, :] = N.T @ v[:, 1]

    def _set_initial_conditions(self, Q_host, func):
        for i in range(self.mesh.num_elements):
            for j in range(self.basis.Np):
                rho, u, v, p = func(self.x_host[i, j], self.y_host[i, j])
                Q_host[i, j, :] = eq.primitive_to_conservative([rho, u, v, p])
