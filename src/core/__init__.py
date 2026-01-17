"""
Core module for the Pazuzu solver.

This module contains the central data structures and configurations used throughout the application, including:
- SimulationState: Manages the memory pool and dynamic simulation data.
- Config: Defines the configuration schema using Pydantic.
- Basis: Handles the nodal DG basis functions and operators.
- BoundaryConditionManager: Manages boundary condition logic and mapping.
"""

from .simulation_state import SimulationState as SimulationState
from .config import PazuzuConfig as PazuzuConfig, SolverType as SolverType
from .basis import Basis as Basis
from .boundary_condition_manager import BoundaryConditionManager as BoundaryConditionManager
