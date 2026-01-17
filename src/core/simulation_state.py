import warp as wp
from typing import Optional, Any
from src.kernels.structs import BoundaryState32, BoundaryState64

class SimulationState:
    """
    Encapsulates the dynamic state of the simulation using a Block-AMR Memory Pool.

    Holds the primary state vector, right-hand side buffer, and intermediate buffers
    allocated for the maximum possible number of blocks (max_blocks).

    Attributes:
        q (wp.array): Current state vector (max_blocks, Np, 4).
        rhs (wp.array): Right-hand side (max_blocks, Np, 4).
        active_block_indices (wp.array): Indices of active blocks. A compact list of indices pointing to valid slots in the memory pool.
        num_active_blocks (int): The current number of active blocks in the simulation.
        t (float): Current simulation time.
        step (int): Current time step index.
        x (wp.array): Physical X coordinates of solution points (max_blocks, Np).
        y (wp.array): Physical Y coordinates of solution points (max_blocks, Np).
        bc_mask (wp.array): Boundary Condition ID for each face (max_blocks, 4). -1 if internal or default.
        bc_data (wp.array): Array of BoundaryState structs containing parameters.
    """
    def __init__(self, Np: int, dtype: Any, device: str, max_blocks: int, scalar_dtype: Any=wp.float32, use_filtering: bool=False, basis: Optional[Any]=None):
        """
        Allocates simulation buffers.

        Args:
            Np (int): Number of solution points per element.
            dtype (wp.dtype): Warp data type for State (e.g., wp.vec4).
            device (str): Compute device ("cpu" or "cuda").
            max_blocks (int): Size of the memory pool.
            scalar_dtype (wp.dtype): Warp data type for Scalars (e.g., wp.float32).
            use_filtering (bool): Whether to allocate a buffer for filtering.
            basis (Basis, optional): The basis object for over-integration.
        """
        self.device = device
        self.max_blocks = max_blocks
        self.num_active_blocks = 0
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

        # --- [NEW] Immersed Boundary Method State ---
        # Signed Distance Field (phi): >0 Fluid, <0 Solid, =0 Interface
        self.phi = wp.zeros(self.pool_shape, dtype=scalar_dtype, device=device)
        # --------------------------------------------

        # Block Management
        self.active_block_indices = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        # Neighbors: (MAX_BLOCKS, 4). Indices: 0:Left, 1:Right, 2:Bottom, 3:Top
        # Stores pool index of the neighbor. -1 if no neighbor (boundary).
        self.neighbors = wp.full((max_blocks, 4), -1, dtype=wp.int32, device=device)
        
        # Boundary Conditions
        self.bc_mask = wp.full((max_blocks, 4), -1, dtype=wp.int32, device=device)
        
        # Initial dummy BC data to allow kernel type inference
        bc_struct = BoundaryState64 if dtype == wp.vec4d else BoundaryState32
        self.bc_data = wp.zeros(1, dtype=bc_struct, device=device)

        # Time Integration Buffers
        self.q_old = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        self.q_temp = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        self.q_accum = wp.zeros(self.pool_shape, dtype=dtype, device=device)
        
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

        # Shock Capturing / Limiting
        self.element_indicator = wp.zeros(max_blocks, dtype=scalar_dtype, device=device)
        self.solver_mode = wp.zeros(max_blocks, dtype=wp.int32, device=device) # 0: FR, 1: FV

        # Geometry Meta-Data (Managed by Quadtree, stored here for Kernels)
        self.block_levels = wp.zeros(max_blocks, dtype=wp.int32, device=device)
        self.root_bounds = wp.zeros(4, dtype=scalar_dtype, device=device)

        # Scalar State
        self.t = 0.0
        self.step = 0

    def zero_rhs(self):
        """Clears the RHS buffer."""
        self.rhs.zero_()

    def numpy(self):
        """Returns the current state q as a numpy array."""
        return self.q.numpy()
