import warp as wp
from typing import Any
from src.core.config import PazuzuConfig
from src.core.simulation_state import SimulationState
from src.core.basis import Basis
from src.kernels import fr_kernels, fv_kernels, indicator_kernels

class HybridRHS:
    """
    Orchestrates the computation of the Right-Hand Side (RHS) for the hybrid FR/FV solver.
    """
    def __init__(self, config: PazuzuConfig, state: SimulationState, basis: Basis):
        self.config = config
        self.state = state
        self.basis = basis
        
        # Ensure filter matrix is computed if shock capturing is enabled
        if self.config.numerics.shock_capturing.enabled:
            if self.basis.filter_matrix is None:
                self.basis.compute_filter_matrix(
                    alpha=self.config.numerics.shock_capturing.indicator_alpha,
                    order=self.config.numerics.shock_capturing.indicator_order
                )

    def compute_rhs(self, t: float, q: Any, rhs: Any):
        """
        Computes the RHS of the semi-discrete equation: dq/dt = RHS(q).
        
        Args:
            t: Current simulation time.
            q: Current state vector (Warp array).
            rhs: Output RHS vector (Warp array).
        """
        # 1. Clear RHS
        self.state.zero_rhs()
        
        # 2. Shock Detection (if enabled)
        if self.config.numerics.shock_capturing.enabled:
            self._detect_troubled_cells()
        
        # 3. Compute FR Update (High-Order)
        # This kernel will skip blocks where solver_mode == 1
        self._compute_fr_fluxes(t, q, rhs)
        
        # 4. Compute FV Update (Low-Order)
        if self.config.numerics.shock_capturing.enabled:
            self._compute_fv_fluxes(t, q, rhs)

    def _detect_troubled_cells(self):
        """ Runs the Persson-Peraire indicator and updates the solver mode. """
        # Compute Indicator (Se)
        wp.launch(
            kernel=indicator_kernels.compute_persson_peraire,
            dim=self.state.num_active_blocks,
            inputs=[
                self.state.q,
                self.state.active_block_indices,
                self.basis.filter_matrix,
                self.state.element_indicator,
                self.state.num_active_blocks,
                0 # Monitor Density (Component 0)
            ],
            device=self.state.device
        )
        
        # Mark Cells (Update solver_mode)
        wp.launch(
            kernel=indicator_kernels.mark_troubled_cells,
            dim=self.state.num_active_blocks,
            inputs=[
                self.state.element_indicator,
                self.state.active_block_indices,
                self.state.solver_mode,
                self.config.numerics.shock_capturing.indicator_threshold,
                self.state.num_active_blocks
            ],
            device=self.state.device
        )

    def _compute_fr_fluxes(self, t: float, q: Any, rhs: Any):
        """ Launches the FR update kernel. """
        # Dimensions for 1D launch over all nodes
        dim = self.state.num_active_blocks * self.basis.Np
        
        wp.launch(
            kernel=fr_kernels.compute_fr_update,
            dim=dim,
            inputs=[
                q,
                self.state.active_block_indices,
                self.state.neighbors,
                self.state.solver_mode,
                self.state.num_active_blocks,
                rhs,
                # Geometry
                self.basis.nodes_1d,
                self.basis.D1D,
                self.basis.dg_L,
                self.basis.dg_R,
                self.state.root_bounds,
                self.state.block_levels,
                # Physics
                self.config.physics,
                t
            ],
            device=self.state.device
        )

    def _compute_fv_fluxes(self, t: float, q: Any, rhs: Any):
        """ Launches the FV update kernel. """
        dim = self.state.num_active_blocks * self.basis.Np
        
        wp.launch(
            kernel=fv_kernels.compute_fv_update,
            dim=dim,
            inputs=[
                q,
                rhs,
                self.state.active_block_indices,
                self.state.neighbors,
                self.state.solver_mode,
                self.basis.weights_1d,
                self.state.root_bounds, # Need to ensure this exists
                self.state.block_levels,
                self.config.physics
            ],
            device=self.state.device
        )
