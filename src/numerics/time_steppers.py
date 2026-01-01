import warp as wp
from typing import Optional
from src.core.simulation_state import SimulationState
from src.kernels.common_kernels import rk_stage_1, rk_stage_2, rk_stage_3, rk4_stage_update, rk4_final_update

class TimeIntegrator:
    """
    Time Integrators:
    1. Low-Storage SSP-RK3 (3rd Order)
    2. Classic RK4 (4th Order)
    """
    def __init__(self, state: SimulationState):
        self.state = state
        self.dt = 0.0

    def step(self, rhs_function, dt: float, time: float):
        pass

    def step_rk4(self, rhs_function, dt: float, time: float, num_active: int):
        """
        Classic Runge-Kutta 4th Order.
        """
        scalar_dtype = self.state.scalar_dtype
        
        # 1. Clear Accumulator (q_accum)
        # We use the pre-allocated buffer in SimulationState.
        self.state.q_accum.zero_()
        
        # --- Stage 1 ---
        # k1 = f(q_n)
        rhs_function(time, self.state.q, self.state.rhs)
        
        # Update Accumulator: q_accum += dt/6 * k1
        # Update Next Stage Input: q_temp = q_n + 0.5*dt * k1
        wp.launch(
            kernel=rk4_stage_update,
            dim=(num_active, self.state.Np),
            inputs=[
                self.state.q,       # q_old
                self.state.rhs,     # rhs (k1)
                self.state.q_accum, # q_accum
                self.state.q_temp,  # q_next
                scalar_dtype(dt),
                scalar_dtype(1.0/6.0), # weight_accum
                scalar_dtype(0.5)      # weight_next
            ],
            device=self.state.device
        )
        
        # --- Stage 2 ---
        # k2 = f(q_temp)
        rhs_function(time + 0.5*dt, self.state.q_temp, self.state.rhs)
        
        # Update Accumulator: q_accum += dt/3 * k2  (2/6 = 1/3)
        # Update Next Stage Input: q_temp = q_n + 0.5*dt * k2
        wp.launch(
            kernel=rk4_stage_update,
            dim=(num_active, self.state.Np),
            inputs=[
                self.state.q,
                self.state.rhs,     # k2
                self.state.q_accum,
                self.state.q_temp,
                scalar_dtype(dt),
                scalar_dtype(1.0/3.0),
                scalar_dtype(0.5)
            ],
            device=self.state.device
        )
        
        # --- Stage 3 ---
        # k3 = f(q_temp)
        rhs_function(time + 0.5*dt, self.state.q_temp, self.state.rhs)
        
        # Update Accumulator: q_accum += dt/3 * k3
        # Update Next Stage Input: q_temp = q_n + 1.0*dt * k3
        wp.launch(
            kernel=rk4_stage_update,
            dim=(num_active, self.state.Np),
            inputs=[
                self.state.q,
                self.state.rhs,     # k3
                self.state.q_accum,
                self.state.q_temp,
                scalar_dtype(dt),
                scalar_dtype(1.0/3.0),
                scalar_dtype(1.0)
            ],
            device=self.state.device
        )
        
        # --- Stage 4 ---
        # k4 = f(q_temp)
        rhs_function(time + dt, self.state.q_temp, self.state.rhs)
        
        # Final Accumulate: q_accum += dt/6 * k4
        # And construct final result: q_n+1 = q_n + q_accum
        # Wait, q_accum stores the *delta*.
        # So q_n+1 = q_n + q_accum.
        # But we can just add q_accum to q_n directly now.
        
        # Actually, let's just finish accumulation.
        wp.launch(
            kernel=rk4_final_update,
            dim=(num_active, self.state.Np),
            inputs=[
                self.state.rhs,     # k4
                self.state.q_accum,
                scalar_dtype(dt),
                scalar_dtype(1.0/6.0)
            ],
            device=self.state.device
        )
        
        # q_n = q_n + q_accum
        # We can reuse rk_stage_1 logic (q = q + 1.0 * accum) or specialized add.
        # Let's use rk_stage_1 as a generic adder: q_out = q + 1.0 * rhs
        # where q=q_n, rhs=q_accum, dt=1.0.
        wp.launch(
            kernel=rk_stage_1,
            dim=(num_active, self.state.Np),
            inputs=[
                self.state.q,
                self.state.q_accum,
                scalar_dtype(1.0),
                self.state.q
            ],
            device=self.state.device
        )

    def step_ssp_rk3(self, rhs_function, dt: float, time: float, num_active: int):
        scalar_dtype = self.state.scalar_dtype
        
        # --- Stage 1 ---
        # Compute RHS(Q_n) -> stored in self.state.rhs
        rhs_function(time, self.state.q, self.state.rhs)
        
        # Q(1) = Q_n + dt * RHS
        # Store Q(1) in q_temp
        wp.launch(
            kernel=rk_stage_1,
            dim=(num_active, self.state.Np),
            inputs=[
                self.state.q,
                self.state.rhs,
                scalar_dtype(dt),
                self.state.q_temp # Output to Temp
            ],
            device=self.state.device
        )
        
        # --- Stage 2 ---
        # Compute RHS(Q(1))
        rhs_function(time + dt, self.state.q_temp, self.state.rhs)
        
        # Q(2) = 0.75*Q_n + 0.25*(Q(1) + dt * RHS)
        # We can overwrite q_temp with Q(2)
        wp.launch(
            kernel=rk_stage_2,
            dim=(num_active, self.state.Np),
            inputs=[
                self.state.q,       # Q_n
                self.state.q_temp,  # Q(1)
                self.state.rhs,     # RHS(Q(1))
                scalar_dtype(dt),
                self.state.q_temp,  # Output Q(2)
                scalar_dtype(0.75),
                scalar_dtype(0.25)
            ],
            device=self.state.device
        )

        # --- Stage 3 ---
        # Compute RHS(Q(2))
        # Time for stage 3? SSP-RK3 stages are effectively at t, t+dt, t+dt/2.
        # But standard is just t.
        rhs_function(time + 0.5*dt, self.state.q_temp, self.state.rhs)
        
        # Q_n+1 = 1/3*Q_n + 2/3*(Q(2) + dt * RHS)
        # Output to self.state.q (Update state)
        wp.launch(
            kernel=rk_stage_3,
            dim=(num_active, self.state.Np),
            inputs=[
                self.state.q,      # Q_n
                self.state.q_temp, # Q(2)
                self.state.rhs,    # RHS(Q(2))
                scalar_dtype(dt),
                self.state.q,      # Output Q_n+1
                scalar_dtype(1.0/3.0),
                scalar_dtype(2.0/3.0)
            ],
            device=self.state.device
        )