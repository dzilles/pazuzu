import warp as wp
import os
import numpy as np
from src.kernels import common_kernels as ck

class TimeIntegrator:
    """
    Manages the explicit time integration loop for the simulation.

    This class handles the time-stepping process using a Low-Storage Strong Stability Preserving 
    Runge-Kutta 3 (SSP-RK3) scheme. It orchestrates calls to the solver's RHS computation 
    and handles I/O operations at specified intervals.

    Attributes:
        solver (BaseSolver): The physics solver instance (e.g., Euler2DSolver).
        device (str): The Warp compute device ("cpu" or "cuda").
        Q_stage1 (wp.array): Buffer for the first RK intermediate stage.
        Q_stage2 (wp.array): Buffer for the second RK intermediate stage.
    """
    def __init__(self, solver):
        """
        Initializes the TimeIntegrator.

        Args:
            solver (BaseSolver): An initialized solver instance containing the state Q and rhs buffer.
        """
        self.solver = solver
        self.device = solver.device
        
        # Allocate Runge-Kutta intermediate buffers
        # Assumes the solver has already initialized its state vector Q
        self.Q_stage1 = wp.zeros_like(solver.Q)
        self.Q_stage2 = wp.zeros_like(solver.Q)
        
    def solve(self, writer=None, max_steps=None):
        """
        Executes the main time-stepping loop.

        This method performs adaptive time-stepping based on the CFL condition provided by the solver.
        Configuration parameters (t_final, CFL, write_interval) are retrieved from the solver's config.

        Args:
            writer (HDF5Writer, optional): Writer object for saving results. If None, saves .npz files.
            max_steps (int, optional): Maximum number of time steps to run. Defaults to None (unlimited).

        Raises:
            ValueError: If simulation instability (NaNs) is detected.
        """
        # Retrieve configuration
        if not self.solver.config:
            raise ValueError("Solver configuration is missing.")
            
        t_final = self.solver.config.simulation.t_final
        CFL = self.solver.config.numerics.cfl
        write_interval = self.solver.config.io.write_interval
        
        t = 0.0
        step = 0
        last_output_time = 0.0
        output_steps_dir = "data/steps"
        
        # Ensure output directory exists if using legacy .npz output
        if writer is None and not os.path.exists(output_steps_dir):
            os.makedirs(output_steps_dir)
            
        # Use a CUDA stream for asynchronous execution if running on GPU
        use_stream = self.device == "cuda"
        if use_stream:
            stream = wp.Stream(device=self.device)
            
        # Direct references to solver's main state arrays for clarity
        Q = self.solver.Q
        rhs = self.solver.rhs
        
        print(f"Starting simulation: t_final={t_final}, CFL={CFL}, write_interval={write_interval}")

        while t < t_final:
            if max_steps is not None and step >= max_steps:
                print(f"Reached maximum steps ({max_steps}). Stopping.")
                break
                
            # Calculate adaptive time step based on current flow state
            dt = self.solver.calculate_dt(CFL)
            
            # Adjust dt to hit t_final exactly
            if t + dt > t_final:
                dt = t_final - t
            
            dtype = self.solver.dtype_warp
            
            # Cast time variables to Warp scalars for kernel compatibility
            dt_warp = dtype(dt)
            t_warp = dtype(t)
            
            # --- SSP-RK3 Time Stepping Scheme ---
            # Stage 1: Q(1) = Q_n + dt * RHS(Q_n)
            # Stage 2: Q(2) = 0.75 * Q_n + 0.25 * (Q(1) + dt * RHS(Q(1)))
            # Stage 3: Q(n+1) = 1/3 * Q_n + 2/3 * (Q(2) + dt * RHS(Q(2)))
            
            if use_stream:
                with stream:
                    # Stage 1
                    self.solver.compute_rhs(t, dt, Q, rhs)
                    wp.launch(kernel=ck.rk_stage_1, dim=Q.shape, inputs=[Q, rhs, dt_warp], outputs=[self.Q_stage1], device=self.device)
                    
                    # Stage 2
                    c1_s2 = dtype(0.75)
                    c2_s2 = dtype(0.25)
                    self.solver.compute_rhs(t+dt, dt, self.Q_stage1, rhs)
                    wp.launch(kernel=ck.rk_stage_2, dim=Q.shape, inputs=[Q, self.Q_stage1, rhs, dt_warp, self.Q_stage2, c1_s2, c2_s2], device=self.device)
                    
                    # Stage 3 (Final update)
                    c1_s3 = dtype(1.0/3.0)
                    c2_s3 = dtype(2.0/3.0)
                    self.solver.compute_rhs(t+0.5*dt, dt, self.Q_stage2, rhs)
                    wp.launch(kernel=ck.rk_stage_3, dim=Q.shape, inputs=[Q, self.Q_stage2, rhs, dt_warp, Q, c1_s3, c2_s3], device=self.device)
                
                stream.synchronize()
            else:
                # Synchronous execution (CPU)
                
                # Stage 1
                self.solver.compute_rhs(t, dt, Q, rhs)
                wp.launch(kernel=ck.rk_stage_1, dim=Q.shape, inputs=[Q, rhs, dt_warp], outputs=[self.Q_stage1], device=self.device)
                
                # Stage 2
                c1_s2 = dtype(0.75)
                c2_s2 = dtype(0.25)
                self.solver.compute_rhs(t+dt, dt, self.Q_stage1, rhs)
                wp.launch(kernel=ck.rk_stage_2, dim=Q.shape, inputs=[Q, self.Q_stage1, rhs, dt_warp, self.Q_stage2, c1_s2, c2_s2], device=self.device)
                
                # Stage 3
                c1_s3 = dtype(1.0/3.0)
                c2_s3 = dtype(2.0/3.0)
                self.solver.compute_rhs(t+0.5*dt, dt, self.Q_stage2, rhs)
                wp.launch(kernel=ck.rk_stage_3, dim=Q.shape, inputs=[Q, self.Q_stage2, rhs, dt_warp, Q, c1_s3, c2_s3], device=self.device)
                
                wp.synchronize()
            
            # Post-step hook (e.g., for filtering)
            if hasattr(self.solver, 'post_step'):
                self.solver.post_step()

            t += dt
            step += 1
            
            # --- Logging and Output ---
            # Basic logging every 10 steps or if output is due
            should_write = (t - last_output_time) >= write_interval or t >= t_final or step == 1

            if step % 10 == 0 or should_write:
                 print(f"Step: {step}, t = {t:.4f} / {t_final}, dt = {dt:.3e}")

            if should_write:
                last_output_time = t
                # Check for NaNs to detect instability early
                q_np = Q.numpy()
                if np.isnan(q_np).any():
                    print(f"!!! Simulation unstable at step {step} (t={t:.4f})")
                    break
                
                # Write output
                if writer:
                    writer.write_step(step, t, q_np)
                else:
                    # Legacy fallback
                    np.savez(f"{output_steps_dir}/step_{step:04d}.npz", q=q_np)
