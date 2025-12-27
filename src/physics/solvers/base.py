from abc import ABC, abstractmethod
import warp as wp
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from src.core.config import PazuzuConfig

class BaseSolver(ABC):
    """
    Abstract base class for all physics solvers.

    This class defines the interface that any specific physics implementation (e.g., Euler, Navier-Stokes)
    must adhere to. It handles the shared state initialization and provides signatures for the core
    computational steps required by the TimeIntegrator.

    Attributes:
        mesh (Mesh): The computational mesh.
        basis (Basis): The polynomial basis and operators.
        device (str): The compute device ("cpu" or "cuda").
        config (PazuzuConfig, optional): Configuration object.
        state (SimulationState): The simulation state container (initialized by concrete classes).
    """
    def __init__(self, mesh, basis, config: Optional["PazuzuConfig"] = None):
        """
        Initializes the base solver.

        Args:
            mesh (Mesh): The computational mesh.
            basis (Basis): The DG basis functions.
            config (PazuzuConfig, optional): Configuration object. Defaults to None.
        """
        self.mesh = mesh
        self.basis = basis
        self.device = mesh.device
        self.config = config
        
        # State variables (to be initialized by concrete class)
        self.state = None
    
    @abstractmethod
    def initialize(self, initial_condition_func):
        """
        Initializes the state vector Q based on a provided function.

        Args:
            initial_condition_func (callable): A function `f(x, y)` returning primitive variables.
        """
        pass

    @abstractmethod
    def compute_rhs(self, t, dt, q, rhs):
        """
        Computes the right-hand side (RHS) of the semi-discrete equation dQ/dt = RHS(Q).

        This method encapsulates the spatial discretization (volume integrals, surface fluxes).

        Args:
            t (float): Current simulation time.
            dt (float): Current time step size.
            q (wp.array): Input state vector.
            rhs (wp.array): Output buffer where the computed RHS will be stored.
        """
        pass

    @abstractmethod
    def calculate_dt(self, CFL):
        """
        Calculates a stable time step size based on the CFL condition.

        Args:
            CFL (float): The Courant-Friedrichs-Lewy number.

        Returns:
            float: The calculated time step size (dt).
        """
        pass
