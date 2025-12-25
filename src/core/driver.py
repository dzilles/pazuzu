import warp as wp
import os
import numpy as np

class TimeIntegrator:
    """
    Manages the explicit time integration loop for the simulation.

    This class handles the time-stepping process using a provided TimeStepper strategy.
    It orchestrates calls to the solver's RHS computation and handles I/O operations at specified intervals.

    Attributes:
        solver (BaseSolver): The physics solver instance (e.g., Euler2DSolver).
        stepper (TimeStepper): The time-stepping algorithm strategy.
        device (str): The Warp compute device ("cpu" or "cuda").
    """
    def __init__(self, solver, stepper):
        """
        Initializes the TimeIntegrator.

        Args:
            solver (BaseSolver): An initialized solver instance containing the state Q and rhs buffer.
            stepper (TimeStepper): A time-stepper instance (e.g., RK4Stepper).
        """
        self.solver = solver
        self.stepper = stepper
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
        min_dt = self.solver.config.numerics.min_dt
        
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
        
        print(f"Starting simulation: t_final={t_final}, CFL={CFL}, write_interval={write_interval}, min_dt={min_dt:.1e}")

        while t < t_final:
            if max_steps is not None and step >= max_steps:
                print(f"Reached maximum steps ({max_steps}). Stopping.")
                break
                
            # Calculate adaptive time step based on current flow state
            dt = self.solver.calculate_dt(CFL)

            # --- Safety Check: Check if dt is too small ---
            if dt < min_dt:
                raise RuntimeError(f"Simulation aborted: Time step dt={dt:.3e} is below threshold min_dt={min_dt:.3e}. The simulation is likely unstable.")
            
            # Adjust dt to hit t_final exactly
            if t + dt > t_final:
                dt = t_final - t
            
            # --- Time Stepping Strategy ---
            self.stepper.step(self.solver, t, dt)
            
            # Post-step hook (e.g., for filtering)
            if hasattr(self.solver, 'post_step'):
                self.solver.post_step()

            t += dt
            step += 1
            
            # Update state object
            state.t = t
            state.step = step
            
            # --- Logging and Output ---
            # Output every write_interval steps or when simulation finishes
            should_write = (step % write_interval == 0) or (t >= t_final)

            if should_write:
                print(f"Step: {step}, t = {t:.4f} / {t_final}, dt = {dt:.3e}")
                
                # Check for NaNs to detect instability early
                q_np = Q.numpy()
                if np.isnan(q_np).any():
                    raise RuntimeError(f"Simulation became unstable: NaN values detected at step {step} (t={t:.4f})")
                
                # Write output
                if writer:
                    writer.write_step(step, t, q_np)
                else:
                    # Legacy fallback
                    np.savez(f"{output_steps_dir}/step_{step:04d}.npz", q=q_np)