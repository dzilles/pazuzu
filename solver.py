import warp as wp
import numpy as np
import os
from typing import Union
from src.core.config import PazuzuConfig, SolverType
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.core.boundary_condition_manager import BoundaryConditionManager
from src.geometry.quadtree import Quadtree
from src.geometry.ibm import IBMManager
from src.numerics.time_steppers import TimeIntegrator
from src.kernels.structs import EquationParams32, EquationParams64, BoundaryState32, BoundaryState64
from src.kernels.fr_kernels import compute_fr_update
from src.kernels.mortar_kernels import compute_mortar_fluxes
from src.kernels.initial_conditions import init_isentropic_vortex, init_uniform
from src.kernels.common_kernels import check_nan_indirect
from src.kernels.ibm_kernels import generate_sdf_cylinder, generate_sdf_from_mesh, apply_ibm_forcing
from src.kernels.grid_kernels import tag_quadtree_boundaries
from src.io.data_writer import HDF5Writer
import src.kernels.boundary_conditions as bc

from src.kernels.time_step_kernels import compute_max_wave_speed

class PazuzuSolver:
    def __init__(self, config: Union[str, PazuzuConfig]):
        # [NEW] Determine Project Root (folder containing solver.py)
        # This ensures we always know where "pazuzu/" is, regardless of where we run python from.
        self.project_root = os.path.dirname(os.path.abspath(__file__))

        if isinstance(config, str):
            # Pass project_root to resolve paths relative to installation
            self.config = PazuzuConfig.from_yaml(config, project_root=self.project_root)
        else:
            self.config = config
            
        wp.init()
        
        # Resolve automatic device selection
        requested_device = self.config.simulation.device
        if requested_device == "automatic":
            if len(wp.get_cuda_devices()) > 0:
                self.device = "cuda"
            else:
                self.device = "cpu"
        else:
            self.device = requested_device
            
        print(f"Using device: {self.device}")

        # 0. Set Precision
        if self.config.numerics.precision == "double":
            self.scalar_dtype = wp.float64
            self.state_dtype = wp.vec4d
            self.params_struct = EquationParams64
            self.bc_struct = BoundaryState64
        else:
            self.scalar_dtype = wp.float32
            self.state_dtype = wp.vec4
            self.params_struct = EquationParams32
            self.bc_struct = BoundaryState32
        
        # 1. Initialize Basis
        self.basis = Basis(
            polynomial_degree=self.config.numerics.polynomial_order,
            device=self.device,
            dtype=self.scalar_dtype
        )
        
        # Determine max_blocks
        # If AMR is disabled, we must ensure max_blocks is large enough for the static mesh.
        # Often users might leave max_blocks small or default (10000) while requesting a deep uniform grid.
        required_blocks_static = (1 << self.config.mesh.initial_depth) ** 2
        
        if not self.config.amr.enabled:
            # If AMR is disabled, force max_blocks to be exactly what is needed (or slightly more)
            # regardless of what the user put in amr.max_blocks
            self.max_blocks = required_blocks_static
            print(f"AMR disabled: Overriding max_blocks to {self.max_blocks} for level {self.config.mesh.initial_depth}")
        else:
            # If AMR is enabled, respect the limit, but warn if it's too small for start
            self.max_blocks = self.config.amr.max_blocks
            if self.max_blocks < required_blocks_static:
                print(f"Warning: amr.max_blocks ({self.max_blocks}) is less than required for initial depth {self.config.mesh.initial_depth} ({required_blocks_static}). Simulation will likely fail.")

        # 2. Initialize State
        self.state = SimulationState(
            Np=self.basis.Np,
            dtype=self.state_dtype,
            device=self.device,
            max_blocks=self.max_blocks,
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
            max_blocks=self.max_blocks,
            root_bounds=bounds,
            periodic_x=self.config.mesh.periodic_x,
            periodic_y=self.config.mesh.periodic_y,
            dtype=self.scalar_dtype,
            max_depth=self.config.amr.max_depth
        )
        
        # 3.5 Initialize IBM Manager
        self.ibm = IBMManager(self.config.ibm, self.device)
        
        # 4. Initialize Time Integrator
        self.integrator = TimeIntegrator(self.state)
        
        # 5. Physics Parameters
        self._init_physics()
        
        # 5.5 Initialize Boundaries
        self._setup_bcs()
        
        # 6. Initial Condition (Mesh Generation)
        self._init_mesh()
        self._apply_initial_condition()
        
        # 6.5 Initialize IBM SDF
        self._init_ibm_sdf()
        
        # Initial Adaptation
        if self.config.amr.enabled and self.config.amr.refinement_threshold is not None:
            refine_threshold = self.config.amr.refinement_threshold
            print(f"Performing initial adaptation (threshold={refine_threshold})...")
            # If coarsening_threshold is not specified, use a default ratio (e.g., 0.5x refine)
            coarsen_threshold = self.config.amr.coarsening_threshold
            if coarsen_threshold is None:
                coarsen_threshold = refine_threshold * 0.5
                
            self.quadtree.adapt_mesh(self.state, self.basis, refine_threshold, coarsen_threshold, ibm_enabled=self.config.ibm.enabled)
            # Retag boundaries and re-init SDF after adaptation
            self._tag_boundaries()
            self._init_ibm_sdf()

        # 7. Initialize Writer
        # self.config.io.output_dir is now guaranteed to be an absolute path 
        # (e.g. /home/user/pazuzu/output/vortex_test)
        output_dir = self.config.io.output_dir
        
        if output_dir is None:
             output_dir = os.path.join(self.project_root, "output", self.config.case_name)

        print(f"Output Directory: {output_dir}") # Helpful log

        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            
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
        self.params.ramp_up_time = float(self.config.simulation.ramp_up_time)

        # Flux Type
        if hasattr(self.config.numerics, "flux") and self.config.numerics.flux == "hllc":
             self.params.flux_type = 1 # HLLC
        else:
             self.params.flux_type = 0 # Rusanov

        # HLLC Fallback
        self.params.hllc_fallback = int(self.config.numerics.hllc_fallback)

    def _setup_bcs(self):
        """Initializes dynamic boundary conditions."""
        self.bc_manager = BoundaryConditionManager(None, self.config)
        bounds = (
            self.config.mesh.x_min,
            self.config.mesh.y_min,
            self.config.mesh.x_max,
            self.config.mesh.y_max
        )
        
        # Choose precision for BC array
        precision = "double" if self.config.numerics.precision == "double" else "single"
        
        bc_data, bc_indices = self.bc_manager.setup_quadtree_bcs(
            bounds, 
            device=self.device, 
            precision=precision
        )
        
        self.state.bc_data = bc_data
        self.bc_indices = bc_indices
        
    def _tag_boundaries(self):
        """Updates boundary masks for active blocks."""
        if self.bc_indices is None:
            return

        wp.launch(
            kernel=tag_quadtree_boundaries,
            dim=self.quadtree.num_blocks,
            inputs=[
                self.state.bc_mask,
                self.state.active_block_indices,
                self.quadtree.num_blocks,
                self.quadtree.block_morton_codes,
                self.quadtree.block_levels,
                self.bc_indices["left"],
                self.bc_indices["right"],
                self.bc_indices["bottom"],
                self.bc_indices["top"]
            ],
            device=self.device
        )

    def _init_mesh(self):
        # Uniform Refinement to start
        initial_level = self.config.mesh.initial_depth
        self.quadtree.uniform_refine(initial_level, self.state, self.basis)
        self._tag_boundaries()

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
        elif ic_name == "uniform" or ic_name == "rest":
            rho = params.get("rho", self.params.rho_inf)
            u = params.get("u", self.params.u_inf)
            v = params.get("v", self.params.v_inf)
            p = params.get("p", self.params.p_inf)
            
            if ic_name == "rest":
                u = 0.0
                v = 0.0
                print(f"Applying Rest IC (rho={rho}, p={p})...")
            else:
                print(f"Applying Uniform IC (rho={rho}, u={u}, v={v}, p={p})...")
                
            wp.launch(
                kernel=init_uniform,
                dim=(self.quadtree.num_blocks, self.basis.Np),
                inputs=[
                    self.state.q,
                    self.state.active_block_indices,
                    self.quadtree.num_blocks,
                    self.params,
                    self.scalar_dtype(rho),
                    self.scalar_dtype(u),
                    self.scalar_dtype(v),
                    self.scalar_dtype(p)
                ],
                device=self.device
            )
        else:
            print(f"Warning: Unknown IC '{ic_name}'. State zeroed.")

    def _init_ibm_sdf(self):
        """Initializes the Signed Distance Field (phi) based on configuration."""
        if not self.config.ibm.enabled:
            return

        print("Initializing IBM Signed Distance Field...")
        
        # Launch params
        dim = (self.quadtree.num_blocks, self.basis.Np)
        
        if self.config.ibm.mode == "analytical":
            params = self.config.ibm.geometric_params
            cx = params.get("center_x", 0.0)
            cy = params.get("center_y", 0.0)
            r = params.get("radius", 0.5)
            
            wp.launch(
                kernel=generate_sdf_cylinder,
                dim=dim,
                inputs=[
                    self.state.x,
                    self.state.y,
                    self.state.phi,
                    self.state.active_block_indices,
                    self.scalar_dtype(cx),
                    self.scalar_dtype(cy),
                    self.scalar_dtype(r)
                ],
                device=self.device
            )
            
        elif self.config.ibm.mode == "stl_file":
            if self.ibm.mesh is None:
                raise RuntimeError("IBM mode is 'stl_file' but no mesh is loaded.")
                
            invert_flag = 1 if self.config.ibm.invert_inside_outside else 0
            max_dist = 100.0 # Large enough
            
            wp.launch(
                kernel=generate_sdf_from_mesh,
                dim=dim,
                inputs=[
                    self.state.x,
                    self.state.y,
                    self.state.phi,
                    self.state.active_block_indices,
                    self.ibm.mesh.id,
                    self.scalar_dtype(max_dist),
                    invert_flag
                ],
                device=self.device
            )

    def _apply_ibm(self, q_current):
        """Applies IBM Ghost-Cell Forcing if enabled."""
        if not self.config.ibm.enabled:
            return
        
        if self.config.ibm.mode == "stl_file" and self.ibm.mesh is None:
            return

        # Mesh ID: if analytical, we don't have a mesh_id for query.
        # Current implementation of `apply_ibm_forcing` REQUIRES a mesh_id for normal query.
        # TODO: Implement analytical normal query for cylinder mode.
        # For now, only apply forcing if mesh is available (stl_file).
        if self.config.ibm.mode == "analytical":
            # Skipping forcing for analytical mode as kernel requires mesh
            # (Or we need a separate kernel for analytical forcing)
            return

        cutoff = 100.0 # Force ALL solid nodes
        invert_flag = 1 if self.config.ibm.invert_inside_outside else 0
        
        # Default to No-Slip (1) for stability unless Slip (0) is explicitly requested
        if self.config.ibm.boundary_type == "slip":
            boundary_type_id = 0
        else:
            boundary_type_id = 1
            if self.config.ibm.boundary_type != "no_slip":
                print("IBM: Defaulting to No-Slip forcing for numerical stability.")
        
        wp.launch(
            kernel=apply_ibm_forcing,
            dim=(self.quadtree.num_blocks, self.basis.Np),
            inputs=[
                q_current,
                self.state.x,
                self.state.y,
                self.state.phi,
                self.state.active_block_indices,
                wp.uint64(self.ibm.mesh.id),
                self.basis.nodes_1d,
                self.basis.N1,
                self.scalar_dtype(cutoff),
                # New Arguments
                self.quadtree.map_keys,
                self.quadtree.map_values,
                self.quadtree.map_capacity,
                self.quadtree.root_bounds_wp,
                self.quadtree.block_levels,
                invert_flag,
                boundary_type_id,
                self.params  # NEW: Pass params struct
            ],
            device=self.device
        )

    def compute_rhs(self, t, q_in, rhs_out):
        """
        Callback for the Time Integrator.
        """
        # 0. Apply IBM Forcing (Ghost Cells)
        self._apply_ibm(q_in)
        
        # 1. Zero out RHS for accumulation
        rhs_out.zero_()
        
        # 2. Volume and Standard Interface Fluxes
        wp.launch(
            kernel=compute_fr_update,
            dim=self.quadtree.num_blocks * self.basis.Np,
            inputs=[
                q_in,
                self.state.active_block_indices,
                self.state.neighbors,
                self.state.solver_mode,
                self.state.bc_mask,
                self.state.bc_data,
                self.state.x,
                self.state.y,
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
        
        # 3. Mortar Interface Corrections (AMR)
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
            if self.config.amr.enabled:
                refine_interval = self.config.amr.refine_interval
                refine_threshold = self.config.amr.refinement_threshold
                if refine_interval > 0 and self.state.step % refine_interval == 0 and refine_threshold is not None:
                    print(f"Adapting mesh at step {self.state.step}...")
                    coarsen_threshold = self.config.amr.coarsening_threshold
                    if coarsen_threshold is None:
                        coarsen_threshold = refine_threshold * 0.5
                    self.quadtree.adapt_mesh(self.state, self.basis, refine_threshold, coarsen_threshold, ibm_enabled=self.config.ibm.enabled)
                    # Retag boundaries and re-init SDF after adaptation
                    self._tag_boundaries()
                    self._init_ibm_sdf()
            
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