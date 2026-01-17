# Developer Guidelines

## Memory Management
*   **No Per-Element Allocations:** Never allocate arrays per-element. This kills performance and fragments memory.
*   **Global Memory Pool:** Always allocate major arrays (fluxes, state variables) based on `MAX_BLOCKS` at startup.
*   **Active Blocks:** Use the `simulation_state` to manage which blocks are active and map them to the underlying memory pool.

## Kernel Design
*   **Block-Based Iteration:** Kernels should launch over `num_active_blocks`.
*   **Indirection:** Use the `active_block_indices` map to access the global memory pool from within kernels.

## Testing
*   **Frequency:** Run `pytest` frequently.
*   **Independence:** Each phase of the project is designed to be independently verifiable. Ensure your changes don't regress previous phases.
*   **Environment:** Run tests within the configured virtual environment.

