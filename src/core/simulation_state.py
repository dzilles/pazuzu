import warp as wp
import numpy as np

MAX_BLOCKS = 10000

class SimulationState:
    """
    Encapsulates the dynamic state of the simulation using a Block-AMR Memory Pool.

    Holds the primary state vector, right-hand side buffer, and intermediate buffers
    allocated for the maximum possible number of blocks (MAX_BLOCKS).

    Attributes:
        q (wp.array): Current state vector (MAX_BLOCKS, Np, 4).
        rhs (wp.array): Right-hand side (MAX_BLOCKS, Np, 4).
        active_block_indices (wp.array): Indices of active blocks.
        t (float): Current simulation time.
        step (int): Current time step index.
        x (wp.array): Physical X coordinates of solution points (MAX_BLOCKS, Np).
        y (wp.array): Physical Y coordinates of solution points (MAX_BLOCKS, Np).
    """
    def __init__(self, Np, dtype, device, scalar_dtype=wp.float32, max_blocks=MAX_BLOCKS, use_filtering=False, basis=None):
        """
        Allocates simulation buffers.

        Args:
            Np (int): Number of solution points per element.
            dtype (wp.dtype): Warp data type for State (e.g., wp.vec4).
            device (str): Compute device ("cpu" or "cuda").
            scalar_dtype (wp.dtype): Warp data type for Scalars (e.g., wp.float32).
            max_blocks (int): Size of the memory pool.
            use_filtering (bool): Whether to allocate a buffer for filtering.
            basis (Basis, optional): The basis object for over-integration.
        """
        self.device = device
        self.max_blocks = max_blocks
        self.Np = Np
        self.dtype = dtype
        self.scalar_dtype = scalar_dtype
        
        self.pool_shape = (max_blocks, Np)

        # Primary State
        self.q = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        self.rhs = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        
        # Grid Coordinates (Physical)
        self.x = wp.zeros(self.pool_shape, dtype=scalar_dtype, device=device)
        self.y = wp.zeros(self.pool_shape, dtype=scalar_dtype, device=device)

        # Block Management
        self.active_block_indices = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        # Neighbors: (MAX_BLOCKS, 4). Indices: 0:Left, 1:Right, 2:Bottom, 3:Top
        # Stores pool index of the neighbor. -1 if no neighbor (boundary).
        self.neighbors = wp.full((max_blocks, 4), -1, dtype=wp.int32, device=device)
        
        # Time Integration Buffers
        self.q_old = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        self.q_temp = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        
        # Optional Filter Buffer
        self.filter_buffer = None
        if use_filtering:
            self.filter_buffer = wp.zeros(self.pool_shape, dtype=dtype, device=device)
            
        # Over-Integration / Decoupled Flux Buffers
        self.f_x_n = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        self.f_y_n = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        
        self.q_q = None
        self.div_q = None
        if basis is not None and hasattr(basis, 'Nq'):
            quad_shape = (max_blocks, basis.Nq)
            self.q_q = wp.zeros(quad_shape, dtype=dtype, device=device)
            self.div_q = wp.zeros(quad_shape, dtype=dtype, device=device)

        # Limiter Buffers (Per-Block)
        self.q_avg = wp.zeros(max_blocks, dtype=dtype, device=device)
        self.q_min = wp.zeros(max_blocks, dtype=dtype, device=device)
        self.q_max = wp.zeros(max_blocks, dtype=dtype, device=device)
        self.grad_x = wp.zeros(max_blocks, dtype=dtype, device=device)
        self.grad_y = wp.zeros(max_blocks, dtype=dtype, device=device)

        # Scalar State
        self.t = 0.0
        self.step = 0

    def zero_rhs(self):
        """Clears the RHS buffer."""
        self.rhs.zero_()

    def numpy(self):
        """Returns the current state q as a numpy array."""
        return self.q.numpy()
