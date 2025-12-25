import warp as wp
import os
import numpy as np
from src.kernels import common_kernels as ck

class TimeIntegrator:
    """
    Manages the explicit time integration loop for the simulation.

    This class handles the time-stepping process using the Classic Runge-Kutta 4 (RK4) scheme.
    It orchestrates calls to the solver's RHS computation and handles I/O operations at specified intervals.

    Attributes:
        solver (BaseSolver): The physics solver instance (e.g., Euler2DSolver).
        device (str): The Warp compute device ("cpu" or "cuda").
    """
    def __init__(self, solver):
        """
        Initializes the TimeIntegrator.

        Args:
            solver (BaseSolver): An initialized solver instance containing the state Q and rhs buffer.
        """
        self.solver = solver
        self.device = solver.device
        
        # Buffers are now managed by solver.state
        
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
        
        # Use state from solver
        state = self.solver.state
        state.t = 0.0
        state.step = 0
        
        t = state.t
        step = state.step
        
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
        Q = state.q
        rhs = state.rhs
        Q_old = state.q_old
        Q_temp = state.q_temp
        
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
            
            # --- Classic RK4 Time Stepping Scheme ---
            # yn+1 = yn + dt/6 * (k1 + 2k2 + 2k3 + k4)
            
            if use_stream:
                with stream:
                    # Start of step: Copy Q (yn) to Q_old
                    wp.copy(Q_old, Q)
                    
                    # Stage 1: k1
                    # Input: Q_old (yn)
                    # Compute k1 = f(t, Q_old)
                    # Accumulate: Q (initially yn) += dt * k1/6 -> yn + dt*k1/6
                    # Prepare next: Q_temp = Q_old + 0.5 * dt * k1
                    self.solver.compute_rhs(t, dt, Q_old, rhs)
                    weight_accum_s1 = dtype(1.0/6.0)
                    weight_next_s1 = dtype(0.5)
                    wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                              inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s1, weight_next_s1], 
                              device=self.device)
                    
                    # Stage 2: k2
                    # Input: Q_temp (yn + 0.5*dt*k1)
                    # Compute k2 = f(t + 0.5*dt, Q_temp)
                    # Accumulate: Q += dt * 2*k2/6 (k2/3) -> yn + dt*k1/6 + dt*k2/3
                    # Prepare next: Q_temp = Q_old + 0.5 * dt * k2
                    self.solver.compute_rhs(t + 0.5*dt, dt, Q_temp, rhs)
                    weight_accum_s2 = dtype(1.0/3.0)
                    weight_next_s2 = dtype(0.5)
                    wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                              inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s2, weight_next_s2], 
                              device=self.device)
                    
                    # Stage 3: k3
                    # Input: Q_temp (yn + 0.5*dt*k2)
                    # Compute k3 = f(t + 0.5*dt, Q_temp)
                    # Accumulate: Q += dt * 2*k3/6 (k3/3) -> yn + dt*k1/6 + dt*k2/3 + dt*k3/3
                    # Prepare next: Q_temp = Q_old + 1.0 * dt * k3
                    self.solver.compute_rhs(t + 0.5*dt, dt, Q_temp, rhs)
                    weight_accum_s3 = dtype(1.0/3.0)
                    weight_next_s3 = dtype(1.0)
                    wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                              inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s3, weight_next_s3], 
                              device=self.device)
                              
                    # Stage 4: k4
                    # Input: Q_temp (yn + dt*k3)
                    # Compute k4 = f(t + dt, Q_temp)
                    # Accumulate: Q += dt * k4/6 -> yn + dt/6*(k1 + 2k2 + 2k3 + k4)
                    self.solver.compute_rhs(t + dt, dt, Q_temp, rhs)
                    weight_accum_s4 = dtype(1.0/6.0)
                    wp.launch(kernel=ck.rk4_final_update, dim=Q.shape, 
                              inputs=[rhs, Q, dt_warp, weight_accum_s4], 
                              device=self.device)
                
                stream.synchronize()
            else:
                # Synchronous execution (CPU)
                
                # Start of step: Copy Q (yn) to Q_old
                wp.copy(Q_old, Q)
                
                # Stage 1
                self.solver.compute_rhs(t, dt, Q_old, rhs)
                weight_accum_s1 = dtype(1.0/6.0)
                weight_next_s1 = dtype(0.5)
                wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                          inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s1, weight_next_s1], 
                          device=self.device)
                
                # Stage 2
                self.solver.compute_rhs(t + 0.5*dt, dt, Q_temp, rhs)
                weight_accum_s2 = dtype(1.0/3.0)
                weight_next_s2 = dtype(0.5)
                wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                          inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s2, weight_next_s2], 
                          device=self.device)
                
                # Stage 3
                self.solver.compute_rhs(t + 0.5*dt, dt, Q_temp, rhs)
                weight_accum_s3 = dtype(1.0/3.0)
                weight_next_s3 = dtype(1.0)
                wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                          inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s3, weight_next_s3], 
                          device=self.device)
                          
                # Stage 4
                self.solver.compute_rhs(t + dt, dt, Q_temp, rhs)
                weight_accum_s4 = dtype(1.0/6.0)
                wp.launch(kernel=ck.rk4_final_update, dim=Q.shape, 
                          inputs=[rhs, Q, dt_warp, weight_accum_s4], 
                          device=self.device)
                
                wp.synchronize()
            
            # Post-step hook (e.g., for filtering)
            if hasattr(self.solver, 'post_step'):
                self.solver.post_step()

            t += dt
            step += 1
            
            # Update state object
            state.t = t
            state.step = step
            
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