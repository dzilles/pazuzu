import warp as wp
import os
import numpy as np
from src.kernels import common_kernels as ck

class TimeIntegrator:
    def __init__(self, solver):
        self.solver = solver
        self.device = solver.device
        
        # Allocate RK buffers
        # We assume solver has Q initialized
        self.Q_stage1 = wp.zeros_like(solver.Q)
        self.Q_stage2 = wp.zeros_like(solver.Q)
        
    def solve(self, t_final, CFL, log_frequency=10, writer=None):
        t = 0.0
        step = 0
        output_steps_dir = "data/steps"
        
        # Only create if writing to npz (no writer provided)
        if writer is None and not os.path.exists(output_steps_dir):
            os.makedirs(output_steps_dir)
            
        use_stream = self.device == "cuda"
        if use_stream:
            stream = wp.Stream(device=self.device)
            
        Q = self.solver.Q
        rhs = self.solver.rhs
        
        while t < t_final:
            dt = self.solver.calculate_dt(CFL)
            
            # --- RK3 Time Stepping ---
            
            if use_stream:
                with stream:
                    # Stage 1
                    self.solver.compute_rhs(t, dt, Q, rhs)
                    wp.launch(kernel=ck.rk_stage_1, dim=Q.shape, inputs=[Q, rhs, dt], outputs=[self.Q_stage1], device=self.device)
                    
                    # Stage 2
                    self.solver.compute_rhs(t+dt, dt, self.Q_stage1, rhs)
                    wp.launch(kernel=ck.rk_stage_2, dim=Q.shape, inputs=[Q, self.Q_stage1, rhs, dt], outputs=[self.Q_stage2], device=self.device)
                    
                    # Stage 3
                    self.solver.compute_rhs(t+0.5*dt, dt, self.Q_stage2, rhs)
                    wp.launch(kernel=ck.rk_stage_3, dim=Q.shape, inputs=[Q, self.Q_stage2, rhs, dt], outputs=[Q], device=self.device)
                
                stream.synchronize()
            else:
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
            
            if step % log_frequency == 0:
                q_np = Q.numpy()
                if np.isnan(q_np).any():
                    print(f"!!! Simulation unstable at step {step}")
                    break
                
                if writer:
                    writer.write_step(step, t, q_np)
                else:
                    np.savez(f"{output_steps_dir}/step_{step:04d}.npz", q=q_np)
                
                print(f"Step: {step}, t = {t:.4f} / {t_final}, dt = {dt:.3e}")
