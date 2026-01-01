import warp as wp
import numpy as np
import os
from typing import Union
from src.core.config import PazuzuConfig, SolverType
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.geometry.quadtree import Quadtree
from src.numerics.time_steppers import TimeIntegrator
from src.kernels.structs import EquationParams32, EquationParams64
from src.kernels.fr_kernels import compute_fr_update
from src.kernels.mortar_kernels import compute_mortar_fluxes
from src.kernels.initial_conditions import init_isentropic_vortex
from src.kernels.common_kernels import check_nan_indirect
from src.io.data_writer import HDF5Writer

from src.kernels.time_step_kernels import compute_max_wave_speed

class PazuzuSolver:
    def __init__(self, config: Union[str, PazuzuConfig]):
        if isinstance(config, str):
            self.config = PazuzuConfig.from_yaml(config)
        else:
            self.config = config
            
        self.device = self.config.simulation.device
        wp.init()

        # 0. Set Precision
        if self.config.numerics.precision == "double":
            self.scalar_dtype = wp.float64
            self.state_dtype = wp.vec4d
            self.params_struct = EquationParams64
        else:
            self.scalar_dtype = wp.float32
            self.state_dtype = wp.vec4
            self.params_struct = EquationParams32
        
        # 1. Initialize Basis
        self.basis = Basis(
            polynomial_degree=self.config.numerics.polynomial_order,
            device=self.device,
            dtype=self.scalar_dtype
        )
        
        # 2. Initialize State
        self.state = SimulationState(
            Np=self.basis.Np,
            dtype=self.state_dtype,
            device=self.device,
            max_blocks=self.config.amr.max_blocks,
            scalar_dtype=self.scalar_dtype
        )
        
        # 3. Initialize Geometry
        bounds = (
            self.config.mesh.x_min,
            self.config.mesh.y_min,
            self.config.mesh.x_max,
            self.config.mesh.y_max
        )
        self.quadtree = Quadtree(
            device=self.device, 
            max_blocks=self.config.amr.max_blocks,
            root_bounds=bounds,
            periodic_x=self.config.mesh.periodic_x,
            periodic_y=self.config.mesh.periodic_y,
            dtype=self.scalar_dtype,
            max_depth=self.config.amr.max_depth
        )
        
        # 4. Initialize Time Integrator
        self.integrator = TimeIntegrator(self.state)
        
        # 5. Physics Parameters
        self._init_physics()
        
        # 6. Initial Condition (Mesh Generation)
        self._init_mesh()
        self._apply_initial_condition()
        
        # Initial Adaptation
        refine_threshold = self.config.amr.refinement_threshold
        if refine_threshold is not None:
            print(f"Performing initial adaptation (threshold={refine_threshold})...")
            # If coarsening_threshold is not specified, use a default ratio (e.g., 0.5x refine)
            coarsen_threshold = self.config.amr.coarsening_threshold
            if coarsen_threshold is None:
                coarsen_threshold = refine_threshold * 0.5
                
            self.quadtree.adapt_mesh(self.state, self.basis, refine_threshold, coarsen_threshold)

        # 7. Initialize Writer
        output_dir = self.config.io.output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        self.writer = HDF5Writer(os.path.join(output_dir, "results.h5"), self)
        
        # 8. Detection Flags and DT Buffer
        self.nan_flag = wp.zeros(1, dtype=wp.int32, device=self.device)
        self.max_inv_dt = wp.zeros(1, dtype=self.scalar_dtype, device=self.device)
        
        # Write Initial State (Step 0)
        print("Saving initial state...")
        self.writer.write_step(0, 0.0)

    def _init_physics(self):
        # Create the struct instance directly, NOT a wp.array
        p = self.config.physics
        self.params = self.params_struct()
        
        # Populate fields
        self.params.gamma = float(p.gamma)
        self.params.rho_inf = float(p.rho_inf)
        self.params.u_inf = float(p.u_inf)
        self.params.v_inf = float(p.v_inf)
        self.params.p_inf = float(p.p_inf)
        self.params.gas_constant = float(p.gas_constant) if hasattr(p, 'gas_constant') else 1.0
        
        # Floors and constants
        self.params.rho_floor = float(p.rho_floor) if hasattr(p, 'rho_floor') else 1e-8
        self.params.p_floor = float(p.p_floor) if hasattr(p, 'p_floor') else 1e-8
        self.params.half = 0.5
        self.params.one = 1.0
        
        # Zero out others (viscosity etc.) if not used, or set defaults
        self.params.mu = 0.0
        self.params.prandtl = 0.72
        self.params.cp = 1.0
        self.params.epsilon = 1.0e-10

        # Flux Type
        if hasattr(self.config.numerics, "flux") and self.config.numerics.flux == "hllc":
             self.params.flux_type = 1 # HLLC
        else:
             self.params.flux_type = 0 # Rusanov

        # HLLC Fallback
        self.params.hllc_fallback = int(self.config.numerics.hllc_fallback)

        # Note: We do NOT wrap this in wp.array(). 
        # We pass the 'self.params' object directly to wp.launch inputs.

    def _init_mesh(self):
        # Uniform Refinement to start
        initial_level = self.config.amr.initial_depth
        self.quadtree.uniform_refine(initial_level, self.state, self.basis)

    def _apply_initial_condition(self):
        ic_name = ""
        params = {}
        
        if isinstance(self.config.initial_condition, str):
            ic_name = self.config.initial_condition
        else:
            ic_name = self.config.initial_condition.name
            params = self.config.initial_condition.params

        if ic_name == "vortex":
            beta = params.get("beta", 5.0)
            radius = params.get("radius", 1.0)
            center_x = params.get("center_x", 0.0)
            center_y = params.get("center_y", 0.0)
            print(f"Applying Isentropic Vortex IC (beta={beta}, radius={radius}, center={center_x},{center_y})...")
            wp.launch(
                kernel=init_isentropic_vortex,
                dim=(self.quadtree.num_blocks, self.basis.Np),
                inputs=[
                    self.state.x,
                    self.state.y,
                    self.state.q,
                    self.state.active_block_indices,
                    self.quadtree.num_blocks,
                    self.params,
                    self.scalar_dtype(0.0),
                    self.scalar_dtype(beta),
                    self.scalar_dtype(radius),
                    self.scalar_dtype(center_x),
                    self.scalar_dtype(center_y)
                ],
                device=self.device
            )
        else:
            print(f"Warning: Unknown IC '{ic_name}'. State zeroed.")

    def compute_rhs(self, t, q_in, rhs_out):
        """
        Callback for the Time Integrator.
        """
        # 0. Zero out RHS for accumulation
        rhs_out.zero_()
        
        # 1. Volume and Standard Interface Fluxes
        wp.launch(
            kernel=compute_fr_update,
            dim=self.quadtree.num_blocks * self.basis.Np,
            inputs=[
                q_in,
                self.state.active_block_indices,
                self.state.neighbors,
                self.quadtree.num_blocks,
                rhs_out,
                self.basis.nodes_1d,
                self.basis.D1D,
                self.basis.dg_L,
                self.basis.dg_R,
                self.quadtree.root_bounds_wp,
                self.quadtree.block_levels,
                self.params,
                self.scalar_dtype(t)
            ],
            device=self.device
        )
        
        # 2. Mortar Interface Corrections (AMR)
        # We need to use the counter on device.
        num_mortars_host = int(self.quadtree.num_mortars.numpy()[0])
        if num_mortars_host > 0:
            wp.launch(
                kernel=compute_mortar_fluxes,
                dim=num_mortars_host,
                inputs=[
                    q_in,
                    rhs_out,
                    self.quadtree.mortar_list,
                    self.quadtree.num_mortars,
                    self.basis.nodes_1d,
                    self.basis.face_nodes,
                    self.basis.dg_L,
                    self.basis.dg_R,
                    self.basis.P_left,
                    self.basis.P_right,
                    self.basis.R_left,
                    self.basis.R_right,
                    self.quadtree.root_bounds_wp,
                    self.quadtree.block_levels,
                    self.params
                ],
                device=self.device
            )

    def compute_dt(self):
        if self.config.numerics.dt_static is not None:
            return self.config.numerics.dt_static

        self.max_inv_dt.zero_()
        wp.launch(
            kernel=compute_max_wave_speed,
            dim=self.quadtree.num_blocks * self.basis.Np,
            inputs=[
                self.state.q,
                self.state.active_block_indices,
                self.quadtree.num_blocks,
                self.quadtree.block_levels,
                self.quadtree.root_bounds_wp,
                self.params,
                self.max_inv_dt
            ],
            device=self.device
        )
        
        max_wave_metric = self.max_inv_dt.numpy()[0]
        if max_wave_metric < 1e-12:
            return self.config.numerics.dt_init # Fallback or Initial
            
        # Scaling CFL based on polynomial order P
        # For DG methods, stability requires CFL proportional to 1/(2P+1)
        p_order = self.basis.N
        scaling = 1.0 / (2.0 * float(p_order) + 1.0)
        effective_cfl = self.config.numerics.cfl * scaling

        dt = effective_cfl / max_wave_metric
        return dt

    def run(self):
        print(f"Starting simulation: {self.config.case_name}")
        t = 0.0
        
        log_freq = self.config.io.write_interval
        t_final = self.config.simulation.t_final
        max_steps = self.config.simulation.max_steps
        
        # Tolerance to prevent one extra step due to FP errors
        tol = 1e-7

        while t < t_final - tol:
            if self.state.step >= max_steps:
                print(f"Reached maximum steps ({max_steps}). Saving final state and exiting.")
                if self.state.step % log_freq != 0:
                    self.writer.write_step(self.state.step, t)
                return

            dt = self.compute_dt()
            
            # Clamp dt to hit t_final exactly
            if t + dt > t_final:
                dt = t_final - t
            
            # Stability Check
            if dt < self.config.numerics.dt_min:
                raise RuntimeError(f"Aborting: Computed timestep {dt:.2e} is smaller than dt_min {self.config.numerics.dt_min:.2e}. The simulation may be unstable.")

            if self.config.numerics.time_integrator == "rk4":
                self.integrator.step_rk4(
                    self.compute_rhs,
                    dt,
                    t,
                    self.quadtree.num_blocks
                )
            else:
                self.integrator.step_ssp_rk3(
                    self.compute_rhs,
                    dt,
                    t,
                    self.quadtree.num_blocks
                )
            
            t += dt
            self.state.t = t
            self.state.step += 1
            
            # Dynamic AMR
            refine_interval = self.config.amr.refine_interval
            refine_threshold = self.config.amr.refinement_threshold
            if refine_interval > 0 and self.state.step % refine_interval == 0 and refine_threshold is not None:
                print(f"Adapting mesh at step {self.state.step}...")
                coarsen_threshold = self.config.amr.coarsening_threshold
                if coarsen_threshold is None:
                    coarsen_threshold = refine_threshold * 0.5
                self.quadtree.adapt_mesh(self.state, self.basis, refine_threshold, coarsen_threshold)
            
            if self.state.step % log_freq == 0:
                # NaN Check
                self.nan_flag.zero_()
                wp.launch(
                    kernel=check_nan_indirect,
                    dim=self.quadtree.num_blocks * self.basis.Np,
                    inputs=[self.state.q, self.state.active_block_indices, self.nan_flag],
                    device=self.device
                )
                
                if self.nan_flag.numpy()[0] > 0:
                    print(f"FATAL: NaN detected at step {self.state.step}, t={t:.6f}")
                    self.writer.write_step(self.state.step, t)
                    break

                print(f"Step {self.state.step}, Time {t:.4f}, dt {dt:.6f}")
                self.writer.write_step(self.state.step, t)

        # Write final state if not already written
        if self.state.step % log_freq != 0:
            print(f"Final Step {self.state.step}, Time {t:.4f}")
            self.writer.write_step(self.state.step, t)

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        solver = PazuzuSolver(sys.argv[1])
        solver.run()
    else:
        print("Usage: python solver.py <config.yaml>")
