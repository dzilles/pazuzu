import warp as wp
from typing import Optional
from src.core.simulation_state import SimulationState
from src.kernels.common_kernels import rk_stage_1, rk_stage_2, rk_stage_3

class TimeIntegrator:
    """
    Low-Storage Strong Stability Preserving Runge-Kutta (SSP-RK3) Time Integrator.
    
    Update rule:
    Stage 1: Q(1) = Q_n + dt * RHS(Q_n)
    Stage 2: Q(2) = 0.75*Q_n + 0.25*(Q(1) + dt * RHS(Q(1)))
    Stage 3: Q_n+1 = 1/3*Q_n + 2/3*(Q(2) + dt * RHS(Q(2)))
    """
    def __init__(self, state: SimulationState):
        self.state = state
        self.dt = 0.0

    def step(self, rhs_function, dt: float, time: float):
        """
        Performs one time step.
        
        Args:
            rhs_function (callable): Function that computes RHS. Signature: (t, q_in, rhs_out)
            dt (float): Time step size.
            time (float): Current simulation time.
        """
        self.dt = dt
        
        num_active = self.state.active_block_indices.shape[0] # Assuming packed for now, or need num_active count?
        # SimulationState keeps track of active indices, but we need the COUNT of active blocks to launch kernels efficiently.
        # Ideally, `state` should track `num_active`.
        # For Phase 1/2 (Uniform), num_active is tracked by Quadtree.
        # Let's assume passed in or stored in state.
        # UPDATE: Quadtree tracks num_blocks. But TimeIntegrator doesn't have Quadtree.
        # We should add `num_active` to SimulationState or pass it.
        # For now, let's launch over MAX_BLOCKS but use the active check in kernel?
        # Or better, the RHS function usually handles the launch dimensions.
        # The RK kernels are element-wise (block-wise).
        
        # NOTE: rk_stage kernels currently iterate (e, i).
        # We should ideally launch over `num_active * Np`.
        # But we need `num_active`.
        # Let's inspect SimulationState again. It doesn't store num_active explicitly yet (just active_block_indices).
        # I will assume `rhs_function` returns `num_active` or we just rely on MAX_BLOCKS for now with mask?
        # Actually, let's assume `state` has a property or we pass it.
        # Let's pass `num_active` to step.
        pass

    def step_ssp_rk3(self, rhs_function, dt: float, time: float, num_active: int):
        
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
                dt,
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
                dt,
                self.state.q_temp,  # Output Q(2)
                0.75,
                0.25
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
                dt,
                self.state.q,      # Output Q_n+1
                1.0/3.0,
                2.0/3.0
            ],
            device=self.state.device
        )
        
        self.state.t += dt
        self.state.step += 1
