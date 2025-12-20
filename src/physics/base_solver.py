from abc import ABC, abstractmethod
import warp as wp

class BaseSolver(ABC):
    def __init__(self, mesh, basis, config=None):
        self.mesh = mesh
        self.basis = basis
        self.device = mesh.device
        self.config = config
        
        # State variables (to be initialized by concrete class)
        self.Q = None
        self.rhs = None
    
    @abstractmethod
    def initialize(self, initial_condition_func):
        """Initialize the state Q."""
        pass

    @abstractmethod
    def compute_rhs(self, t, dt, q, rhs):
        """Compute RHS = dQ/dt."""
        pass

    @abstractmethod
    def calculate_dt(self, CFL):
        """Calculate stable time step."""
        pass
