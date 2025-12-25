from abc import ABC, abstractmethod
import warp as wp
from src.kernels import common_kernels as ck

class TimeStepper(ABC):
    """
    Abstract base class for time-stepping algorithms.
    """
    def __init__(self, device: str = "cpu"):
        self.device = device
        self.use_stream = device == "cuda"
        if self.use_stream:
            self.stream = wp.Stream(device=self.device)

    @abstractmethod
    def step(self, solver, t, dt):
        """
        Performs a single time step.

        Args:
            solver (BaseSolver): The physics solver instance.
            t (float): Current simulation time.
            dt (float): Time step size.
        """
        pass

class RK4Stepper(TimeStepper):
    """
    Classic Runge-Kutta 4 (RK4) time-stepping scheme.
    """
    def step(self, solver, t, dt):
        state = solver.state
        dtype = solver.dtype_warp
        dt_warp = dtype(dt)
        
        Q = state.q
        rhs = state.rhs
        Q_old = state.q_old
        Q_temp = state.q_temp
        
        if self.use_stream:
            with self.stream:
                self._rk4_logic(solver, t, dt, dt_warp, dtype, Q, rhs, Q_old, Q_temp)
            self.stream.synchronize()
        else:
            self._rk4_logic(solver, t, dt, dt_warp, dtype, Q, rhs, Q_old, Q_temp)
            wp.synchronize()

    def _rk4_logic(self, solver, t, dt, dt_warp, dtype, Q, rhs, Q_old, Q_temp):
        # Start of step: Copy Q (yn) to Q_old
        wp.copy(Q_old, Q)
        
        # Stage 1: k1
        solver.compute_rhs(t, dt, Q_old, rhs)
        weight_accum_s1 = dtype(1.0/6.0)
        weight_next_s1 = dtype(0.5)
        wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                  inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s1, weight_next_s1], 
                  device=self.device)
        
        # Stage 2: k2
        solver.compute_rhs(t + 0.5*dt, dt, Q_temp, rhs)
        weight_accum_s2 = dtype(1.0/3.0)
        weight_next_s2 = dtype(0.5)
        wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                  inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s2, weight_next_s2], 
                  device=self.device)
        
        # Stage 3: k3
        solver.compute_rhs(t + 0.5*dt, dt, Q_temp, rhs)
        weight_accum_s3 = dtype(1.0/3.0)
        weight_next_s3 = dtype(1.0)
        wp.launch(kernel=ck.rk4_stage_update, dim=Q.shape, 
                  inputs=[Q_old, rhs, Q, Q_temp, dt_warp, weight_accum_s3, weight_next_s3], 
                  device=self.device)
                  
        # Stage 4: k4
        solver.compute_rhs(t + dt, dt, Q_temp, rhs)
        weight_accum_s4 = dtype(1.0/6.0)
        wp.launch(kernel=ck.rk4_final_update, dim=Q.shape, 
                  inputs=[rhs, Q, dt_warp, weight_accum_s4], 
                  device=self.device)
