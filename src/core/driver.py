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
        
    def solve(self, t_final, CFL, log_frequency=10, writer=None, max_steps=None):
        """
        Executes the main time-stepping loop until t_final is reached.

        This method performs adaptive time-stepping based on the CFL condition provided by the solver.
        It executes the 3-stage SSP-RK3 integration and handles data output.

        Args:
            t_final (float): The final simulation time to reach.
            CFL (float): The CFL number (Courant-Friedrichs-Lewy) for time step stability.
            log_frequency (int, optional): Frequency of logging/saving steps (every N steps). Defaults to 10.
            writer (HDF5Writer, optional): Writer object for saving results. If None, saves .npz files.
            max_steps (int, optional): Maximum number of time steps to run. Defaults to None (unlimited).

        Raises:
            ValueError: If simulation instability (NaNs) is detected.
        """
        t = 0.0
        step = 0
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
        
        while t < t_final:
            if max_steps is not None and step >= max_steps:
                print(f"Reached maximum steps ({max_steps}). Stopping.")
                break
                
            # Calculate adaptive time step based on current flow state
            dt = self.solver.calculate_dt(CFL)
            
            # --- SSP-RK3 Time Stepping Scheme ---
            # Stage 1: Q(1) = Q_n + dt * RHS(Q_n)
            # Stage 2: Q(2) = 0.75 * Q_n + 0.25 * (Q(1) + dt * RHS(Q(1)))
            # Stage 3: Q(n+1) = 1/3 * Q_n + 2/3 * (Q(2) + dt * RHS(Q(2)))
            
            if use_stream:
                with stream:
                    # Stage 1
                    self.solver.compute_rhs(t, dt, Q, rhs)
                    wp.launch(kernel=ck.rk_stage_1, dim=Q.shape, inputs=[Q, rhs, dt], outputs=[self.Q_stage1], device=self.device)
                    
                    # Stage 2
                    self.solver.compute_rhs(t+dt, dt, self.Q_stage1, rhs)
                    wp.launch(kernel=ck.rk_stage_2, dim=Q.shape, inputs=[Q, self.Q_stage1, rhs, dt], outputs=[self.Q_stage2], device=self.device)
                    
                    # Stage 3 (Final update)
                    self.solver.compute_rhs(t+0.5*dt, dt, self.Q_stage2, rhs)
                    wp.launch(kernel=ck.rk_stage_3, dim=Q.shape, inputs=[Q, self.Q_stage2, rhs, dt], outputs=[Q], device=self.device)
                
                stream.synchronize()
            else:
                # Synchronous execution (CPU)
                
                # Stage 1
                self.solver.compute_rhs(t, dt, Q, rhs)
                wp.launch(kernel=ck.rk_stage_1, dim=Q.shape, inputs=[Q, rhs, dt], outputs=[self.Q_stage1], device=self.device)
                
                # Stage 2
                self.solver.compute_rhs(t+dt, dt, self.Q_stage1, rhs)
                wp.launch(kernel=ck.rk_stage_2, dim=Q.shape, inputs=[Q, self.Q_stage1, rhs, dt], outputs=[self.Q_stage2], device=self.device)
                
                # Stage 3
                self.solver.compute_rhs(t+0.5*dt, dt, self.Q_stage2, rhs)
                wp.launch(kernel=ck.rk_stage_3, dim=Q.shape, inputs=[Q, self.Q_stage2, rhs, dt], outputs=[Q], device=self.device)
                
                wp.synchronize()
            
            t += dt
            step += 1
            
            # --- Logging and Output ---
            if step % log_frequency == 0:
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
                
                print(f"Step: {step}, t = {t:.4f} / {t_final}, dt = {dt:.3e}")
