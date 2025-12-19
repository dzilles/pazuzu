from .mesh import Mesh
from .basis import Basis
from . import equations as eq
from . import warp_kernels as wk

import numpy as np
import warp as wp
import os

class DGSolver:
    def __init__(self, mesh: Mesh, basis: Basis, initial_conditions_func):
        """
        Initializes the Discontinuous Galerkin solver.
        """
        self.mesh = mesh
        self.basis = basis
        self.device = mesh.device
        
        # --- Host-side data (for setup) ---
        self.geo_factors_host = np.zeros((mesh.num_elements, 5), dtype=np.float32) # dx_dr, dx_ds, dy_dr, dy_ds, J
        self.face_geo_factors_host = np.zeros((mesh.num_elements, 4, 3), dtype=np.float32) # nx, ny, J_face
        self.x_host = np.zeros((mesh.num_elements, basis.Np), dtype=np.float32)
        self.y_host = np.zeros((mesh.num_elements, basis.Np), dtype=np.float32)
        
        self._calculate_geometric_factors()
        self._map_nodes_to_physical_space()

        # --- Device-side data (for computation) ---
        self.geo_factors = wp.array(self.geo_factors_host, dtype=wp.float32, device=self.device)
        self.face_geo_factors = wp.array(self.face_geo_factors_host, dtype=wp.float32, device=self.device)
        self.x = wp.array(self.x_host, dtype=wp.float32, device=self.device)
        self.y = wp.array(self.y_host, dtype=wp.float32, device=self.device)

        # Inverse of Mass Matrix (diagonal)
        weights_2d = np.kron(self.basis.weights_1d.numpy(), self.basis.weights_1d.numpy())
        M_inv_diag_host = 1.0 / weights_2d
        self.M_inv_diag = wp.array(M_inv_diag_host, dtype=wp.float32, device=self.device)

        # Solution field
        Q_host = np.zeros((mesh.num_elements, basis.Np, 4), dtype=np.float32)
        self._set_initial_conditions(Q_host, initial_conditions_func)
        self.Q = wp.array(Q_host, dtype=wp.vec4, device=self.device)
        self.Q_stage1 = wp.zeros_like(self.Q)
        self.Q_stage2 = wp.zeros_like(self.Q)
        self.Q_rk = self.Q # Input for the RHS computation
        self.rhs = wp.zeros_like(self.Q)

        self.max_wave_speed = wp.zeros(1, dtype=wp.float32, device=self.device)


    def _calculate_geometric_factors(self):
        """
        Calculates geometric factors for affine quadrilateral elements.
        The transformation is x(r,s) = sum(N_i(r,s) * x_i), where N_i are bilinear shape functions.
        For affine elements, the derivatives (and thus the Jacobian) are constant.
        """
        # Derivatives of bilinear shape functions at the center (r=0, s=0)
        dNd_r_center = 0.25 * np.array([-1, 1, 1, -1])
        dNd_s_center = 0.25 * np.array([-1, -1, 1, 1])
        
        for i in range(self.mesh.num_elements):
            # vertices are on host, so we use the host array
            v = self.mesh.vertices_host[i, :, :] 
            dx_dr = dNd_r_center @ v[:, 0]
            dx_ds = dNd_s_center @ v[:, 0]
            dy_dr = dNd_r_center @ v[:, 1]
            dy_ds = dNd_s_center @ v[:, 1]
            J = dx_dr * dy_ds - dx_ds * dy_dr

            self.geo_factors_host[i, 0] = dx_dr
            self.geo_factors_host[i, 1] = dx_ds
            self.geo_factors_host[i, 2] = dy_dr
            self.geo_factors_host[i, 3] = dy_ds
            self.geo_factors_host[i, 4] = J
            
            # --- Face geometric factors ---
            # Face 0 (bottom): s=-1, tangent (dx/dr, dy/dr)
            nx, ny = dy_dr, -dx_dr
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 0, :] = [nx/J_face, ny/J_face, J_face]

            # Face 1 (right): r=1, tangent (dx/ds, dy/ds)
            nx, ny = dy_ds, -dx_ds
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 1, :] = [nx/J_face, ny/J_face, J_face]
            
            # Face 2 (top): s=1, tangent (dx/dr, dy/dr)
            nx, ny = -dy_dr, dx_dr # Normal points outwards
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 2, :] = [nx/J_face, ny/J_face, J_face]

            # Face 3 (left): r=-1, tangent (dx/ds, dy/ds)
            nx, ny = -dy_ds, dx_ds # Normal points outwards
            J_face = np.sqrt(nx**2 + ny**2)
            self.face_geo_factors_host[i, 3, :] = [nx/J_face, ny/J_face, J_face]


    def _map_nodes_to_physical_space(self):
        """ Maps nodal points from the reference to the physical element. """
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

    def _compute_rhs(self):
        """ Computes the right-hand side of the ODE using Warp kernels. """
        self.rhs.zero_()
        
        # --- Volume Integral ---
        wp.launch(
            kernel=wk.compute_volume_term,
            dim=(self.mesh.num_elements, self.basis.Np),
            inputs=[
                self.Q_rk,
                self.rhs,
                self.basis.Dr,
                self.basis.Ds,
                self.mesh.rx,
                self.mesh.sy,
                self.basis.Np
            ],
            device=self.device
        )
        
        # --- Surface Integral ---
        wp.launch(
            kernel=wk.compute_surface_term,
            dim=self.mesh.num_elements,  # Nur 1 Thread pro Element (Kernel loopt über Faces)
            inputs=[
                self.Q_rk,
                self.rhs,
                self.mesh.connectivity,
                self.basis.face_nodes,
                self.basis.LIFT,
                self.mesh.Js_x,
                self.mesh.Js_y,
                self.mesh.J,
                self.basis.Nfp
            ],
            device=self.device
        )

    def solve(self, t_final, CFL, log_frequency=10):
        """ The main time-stepping loop. """
        t = 0.0
        step = 0
        
        # --- Debugging: Directory to save step-by-step results ---
        output_steps_dir = "data/steps"
        if not os.path.exists(output_steps_dir):
            os.makedirs(output_steps_dir)
        
        use_stream = self.device == "cuda"
        if use_stream:
            stream = wp.Stream(device=self.device)

        while t < t_final:
            
            # --- Adaptive dt calculation ---
            self.max_wave_speed.zero_()
            wp.launch(
                kernel=wk.compute_max_wave_speed,
                dim=(self.mesh.num_elements, self.basis.Np),
                inputs=[self.Q, self.max_wave_speed],
                device=self.device
            )
            max_speed = self.max_wave_speed.numpy()[0]
            print(f"Max speed: {max_speed:.4e}")
            
            # Ensure max_speed is not zero to avoid division by zero
            # Also, handle the case where the simulation might start with zero velocity
            if max_speed < 1.0e-6:
                max_speed = 1.0 

            # Smallest distance between any two nodes in the mesh (characteristic length)
            # A simpler but less robust way is to use dx = mesh.lx / mesh.nx
            dx = self.mesh.dx 
            dt = CFL * dx / max_speed
            print(f"Calculated dt: {dt:.4e}")

            # --- End of adaptive dt calculation ---

            if use_stream:
                with stream:
                    # Standard SSP-RK3 scheme
                    self.Q_rk = self.Q
                    self._compute_rhs()
                    wp.launch(kernel=wk.rk_stage_1, dim=self.Q.shape, inputs=[self.Q, self.rhs, dt], outputs=[self.Q_stage1], device=self.device)
                    
                    self.Q_rk = self.Q_stage1
                    self._compute_rhs()
                    wp.launch(kernel=wk.rk_stage_2, dim=self.Q.shape, inputs=[self.Q, self.Q_stage1, self.rhs, dt], outputs=[self.Q_stage2], device=self.device)

                    self.Q_rk = self.Q_stage2
                    self._compute_rhs()
                    wp.launch(kernel=wk.rk_stage_3, dim=self.Q.shape, inputs=[self.Q, self.Q_stage2, self.rhs, dt], outputs=[self.Q], device=self.device)
                
                stream.synchronize()
            else: # Synchronous execution for CPU
                # 1st stage
                self.Q_rk = self.Q
                self._compute_rhs()
                wp.launch(kernel=wk.rk_stage_1, dim=self.Q.shape, inputs=[self.Q, self.rhs, dt], outputs=[self.Q_stage1], device=self.device)
                
                # 2nd stage
                self.Q_rk = self.Q_stage1
                self._compute_rhs()
                wp.launch(kernel=wk.rk_stage_2, dim=self.Q.shape, inputs=[self.Q, self.Q_stage1, self.rhs, dt], outputs=[self.Q_stage2], device=self.device)

                # 3rd stage
                self.Q_rk = self.Q_stage2
                self._compute_rhs()
                wp.launch(kernel=wk.rk_stage_3, dim=self.Q.shape, inputs=[self.Q, self.Q_stage2, self.rhs, dt], outputs=[self.Q], device=self.device)

                # For CPU, a manual synchronize is good practice to ensure completion
                wp.synchronize()


            t += dt
            step += 1

            # --- Debugging: Check for NaNs and save step data ---
            q_np = self.Q.numpy()
            if np.isnan(q_np).any():
                print(f"!!! Simulation became unstable with NaN values at step {step} (t={t:.4f}).")
                print(f"!!! The last valid data was saved at step {step-1}.")
                break
            
            # Save step data
            step_filename = os.path.join(output_steps_dir, f"step_{step:04d}.npz")
            np.savez(step_filename, q=q_np, x=self.x.numpy(), y=self.y.numpy(), vertices=self.mesh.vertices_host)

            if step % log_frequency == 0:
              print(f"Step: {step}, t = {t:.4f} / {t_final}, dt = {dt:.3e}")
