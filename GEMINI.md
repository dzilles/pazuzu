# Gemini Context: Warp-based 2D Euler CFD Solver

## Project Overview
This project is a high-performance **2D Discontinuous Galerkin (DG) Euler Solver** written in Python, utilizing **Nvidia Warp** (`warp-lang`) for hardware acceleration (CPU/CUDA). It solves the compressible Euler equations on both structured (Cartesian) and unstructured meshes.

## Key Technologies
*   **Language:** Python 3.13+
*   **Compute Backend:** [Nvidia Warp](https://github.com/NVIDIA/warp) (CPU/GPU acceleration)
*   **Numerical Methods:** 
    *   Discontinuous Galerkin (DG) Spatial Discretization
    *   SSP-RK3 Time Integration (Low-Storage Strong Stability Preserving Runge-Kutta 3)
    *   Lax-Friedrichs / Rusanov Numerical Flux
*   **I/O:** HDF5 (via `h5py`) and XMF for visualization (Paraview compatible).
*   **Meshing:** Gmsh (`gmsh` Python API) support.

## Project Structure

```text
/
├── run_simulation.py       # Main entry point for running simulations
├── config/                 # Configuration files
│   └── default.yaml        # Default simulation parameters
├── src/
│   ├── core/               # Core infrastructure
│   │   ├── driver.py       # Time integration (TimeIntegrator class)
│   │   └── basis.py        # Polynomial basis definitions
│   ├── geometry/           # Mesh handling
│   │   └── mesh.py         # Mesh class (Cartesian & Unstructured)
│   ├── physics/            # Physics solvers and equations
│   │   ├── euler_2d.py     # Main Euler2DSolver class
│   │   ├── equations.py    # Flux functions, state conversions
│   │   └── initial_conditions.py
│   ├── kernels/            # Warp kernels (compiled code)
│   │   ├── common_kernels.py # RK stages, basic math
│   │   └── euler_kernels.py  # Flux loops, boundary conditions
│   └── io/                 # Input/Output
│       └── data_writer.py  # HDF5/XMF writer
└── tests/                  # Verification and testing
    ├── run_tests.py        # Test runner
    └── verification/       # Standard CFD test cases (Vortex, Channel, Cylinder)
```

## Setup & Environment
The project relies on a pre-configured virtual environment located in `.venv`. You should use the python interpreter within this environment to ensure all dependencies are available.

**Dependencies (installed in `.venv`):**
*   `warp-lang`
*   `numpy`
*   `h5py`
*   `gmsh`
*   `pyyaml`
*   `matplotlib` (for visualization scripts)
*   `pytest` (for testing)

## Usage

### Running a Simulation
To run a simulation, use the python interpreter from the virtual environment.

```bash
# Run with default configuration
./.venv/bin/python run_simulation.py

# Run with a specific configuration
./.venv/bin/python run_simulation.py config/my_simulation.yaml
```

### Configuration (`.yaml`)
Configuration files control all aspects of the simulation. Key sections:
*   **`simulation`**: Time stepping (`t_final`, `cfl`), output settings (`output_dir`, `log_frequency`), and device selection (`cpu` or `cuda`).
*   **`mesh`**: 
    *   `type: cartesian` (needs `nx`, `ny`, bounds)
    *   `type: unstructured` (needs `filename` to a `.msh` file)
*   **`numerical`**: `polynomial_degree` for the DG scheme.
*   **`initial_condition`**: `name` (e.g., "vortex", "uniform").

### Visualization
Results are saved as `.h5` files with accompanying `.xmf` descriptors in the `output_dir`. These can be opened directly in **Paraview**.

There is also a helper script:
```bash
./.venv/bin/python scripts/visualize.py --file data/results.h5
```

## Development & Testing

### Running Tests
The project includes a test suite in the `tests/` directory. You should use `pytest` from the virtual environment.

```bash
# Run all tests (unit + verification)
./.venv/bin/python -m pytest

# Run a specific verification test
./.venv/bin/python tests/verification/channel_flow/verify.py
```

### Verification Cases
Specific physics verification cases (like the Isentropic Vortex or Cylinder Flow) are located in `tests/verification/`. These usually have their own generation and verification scripts (e.g., `generate_mesh.py`, `verify.py`).

## Coding Conventions
*   **Warp Usage:** Compute-heavy loops should be implemented as Warp kernels in `src/kernels/` and launched from the solver classes.
*   **Type Safety:** Use Warp's type system (`wp.array`, `wp.vec2`, etc.) within kernels.
*   **I/O:** Heavy data should be written to HDF5. Text output should be minimal (logging).
