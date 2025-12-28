import warp as wp
import numpy as np
import os
from src.core.config import PazuzuConfig, SolverType
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.geometry.quadtree import Quadtree
from src.core.time_integrator import TimeIntegrator
from src.physics.laws.euler import EquationParams32
from src.kernels.fr_kernels import compute_fr_update
from src.kernels.structs import EquationParams32 as ParamsStruct
from src.kernels.initial_conditions import init_isentropic_vortex
from src.io.data_writer import HDF5Writer

class PazuzuSolver:
    def __init__(self, config_path: str):
        self.config = PazuzuConfig.from_yaml(config_path)
        self.device = self.config.simulation.device
        wp.init()
        
        # 1. Initialize Basis
        self.basis = Basis(
            polynomial_degree=self.config.numerics.polynomial_order,
            device=self.device
        )
        
        # 2. Initialize State
        self.state = SimulationState(
            Np=self.basis.Np,
            dtype=wp.vec4,
            device=self.device,
            max_blocks=self.config.amr.max_blocks,
            scalar_dtype=wp.float32 # Assume single precision for now
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
            periodic_y=self.config.mesh.periodic_y
        )
        
        # 4. Initialize Time Integrator
        self.integrator = TimeIntegrator(self.state)
        
        # 5. Physics Parameters
        self._init_physics()
        
        # 6. Initial Condition (Mesh Generation)
        self._init_mesh()
        self._apply_initial_condition()

        # 7. Initialize Writer
        output_dir = self.config.io.output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        self.writer = HDF5Writer(os.path.join(output_dir, "results.h5"), self)
        
        # Write Initial State (Step 0)
        print("Saving initial state...")
        self.writer.write_step(0, 0.0)

    def _init_physics(self):
        # Convert config physics to Warp struct
        p = self.config.physics
        
        params_np = np.zeros(1, dtype=ParamsStruct.numpy_dtype())
        params_np[0]['gamma'] = p.gamma
        params_np[0]['gas_constant'] = p.gas_constant
        params_np[0]['rho_floor'] = p.rho_floor
        params_np[0]['p_floor'] = p.p_floor
        params_np[0]['half'] = 0.5
        params_np[0]['one'] = 1.0
        # ... others
        params_np[0]['rho_inf'] = p.rho_inf
        params_np[0]['u_inf'] = p.u_inf
        params_np[0]['v_inf'] = p.v_inf
        params_np[0]['p_inf'] = p.p_inf
        
        self.params = wp.array(params_np, dtype=ParamsStruct, device=self.device)

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
            print(f"Applying Isentropic Vortex IC (beta={beta}, radius={radius})...")
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
                    0.0,
                    float(beta),
                    float(radius)
                ],
                device=self.device
            )
        else:
            print(f"Warning: Unknown IC '{ic_name}'. State zeroed.")

    def compute_rhs(self, t, q_in, rhs_out):
        """
        Callback for the Time Integrator.
        """
        # Zero RHS? The kernel overwrites, but accumulation might be needed if multiple physics?
        # FR Kernel overwrites.
        
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
                3, # Fixed Level for now
                self.params,
                t
            ],
            device=self.device
        )

    def run(self):
        print(f"Starting simulation: {self.config.case_name}")
        t = 0.0
        dt = 0.001 # Fixed DT for now. 
        # Ideally, calculate DT based on CFL: dt = CFL * dx / max_wave_speed
        
        log_freq = self.config.io.write_interval
        t_final = self.config.simulation.t_final
        
        while t < t_final:
            self.integrator.step_ssp_rk3(
                self.compute_rhs,
                dt,
                t,
                self.quadtree.num_blocks
            )
            t += dt
            self.state.step += 1
            
            if self.state.step % log_freq == 0:
                print(f"Step {self.state.step}, Time {t:.4f}")
                self.writer.write_step(self.state.step, t)

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        solver = PazuzuSolver(sys.argv[1])
        solver.run()
    else:
        print("Usage: python pazuzu.py <config.yaml>")
