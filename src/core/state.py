import warp as wp
import numpy as np

class SimulationState:
    """
    Encapsulates the dynamic state of the simulation.

    Holds the primary state vector, right-hand side buffer, and intermediate buffers
    required for time integration (e.g., RK4 stages). managing allocation and device placement.

    Attributes:
        q (wp.array): Current state vector (Conservative variables).
        rhs (wp.array): Right-hand side (time derivative).
        q_old (wp.array): State at the beginning of the time step.
        q_temp (wp.array): Intermediate state for multi-stage integration.
        filter_buffer (wp.array): Optional buffer for solution filtering.
        t (float): Current simulation time.
        step (int): Current time step index.
    """
    def __init__(self, shape, dtype, device, use_filtering=False):
        """
        Allocates simulation buffers.

        Args:
            shape (tuple): Shape of the state array (NumElements, Np, 4).
            dtype (wp.dtype): Warp data type (e.g., wp.vec4).
            device (str): Compute device ("cpu" or "cuda").
            use_filtering (bool): Whether to allocate a buffer for filtering.
        """
        self.device = device
        self.shape = shape
        self.dtype = dtype
        
        # Primary State
        self.q = wp.zeros(shape, dtype=dtype, device=device)
        self.rhs = wp.zeros(shape, dtype=dtype, device=device)
        
        # Time Integration Buffers
        self.q_old = wp.zeros(shape, dtype=dtype, device=device)
        self.q_temp = wp.zeros(shape, dtype=dtype, device=device)
        
        # Optional Filter Buffer
        self.filter_buffer = None
        if use_filtering:
            self.filter_buffer = wp.zeros(shape, dtype=dtype, device=device)
            
        # Scalar State
        self.t = 0.0
        self.step = 0

    def zero_rhs(self):
        """Clears the RHS buffer."""
        self.rhs.zero_()

    def numpy(self):
        """Returns the current state Q as a numpy array."""
        return self.q.numpy()
