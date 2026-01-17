<p align="center">
  <img src="assets/logo.png" alt="Pazuzu Logo" width="200">
</p>

# Pazuzu: Warp-based 2D DG CFD Solver

Pazuzu is a high-performance **2D Discontinuous Galerkin (DG) Solver** for compressible fluid dynamics, written in Python and accelerated by **Nvidia Warp** (`warp-lang`) for both CPU and CUDA backends.

It supports both the **Euler equations** and the **Navier-Stokes equations** on structured (Cartesian) meshes.

## Key Features

*   **Compute Backend:** [Nvidia Warp](https://github.com/NVIDIA/warp) for JIT-compiled CPU/GPU acceleration.
*   **Physics Solvers:**
    *   2D Compressible Euler Equations.
    *   2D Compressible Navier-Stokes Equations.
*   **Numerical Methods:**
    *   **Spatial Discretization:** Discontinuous Galerkin (DG) with arbitrary polynomial degree.
    *   **Time Integration:** SSP-RK3 (Strong Stability Preserving Runge-Kutta 3).
    *   **Numerical Fluxes:** Lax-Friedrichs / Rusanov.
*   **Meshing:**
    *   Built-in Cartesian mesh generator.
*   **I/O:**
    *   HDF5 data storage (via `h5py`).
    *   XMF descriptors for seamless visualization in **ParaView**.

## Project Structure

```text
├── solver.py               # Main entry point
├── config/                 # YAML configuration files
├── src/
│   ├── core/               # Driver, state, and basis functions
│   ├── geometry/           # Mesh handling (Cartesian)
│   ├── physics/            # Euler and Navier-Stokes solver logic
│   ├── kernels/            # Warp kernels (fluxes, BCs, RK stages)
│   ├── numerics/           # Time stepping algorithms
│   └── io/                 # HDF5/XMF data writers
└── tests/
    ├── unit/               # Unit tests for kernels and core components
    └── verification/       # Standard CFD test cases (Vortex, Channel, Cylinder, etc.)
```

## Setup & Installation

### Prerequisites

*   Python 3.13+
*   Nvidia GPU (optional, for CUDA acceleration)

### Environment Setup

It is recommended to use the provided virtual environment:

```bash
# Create and activate virtual environment (if not already present)
python -m venv .venv
source .venv/bin/bin/activate  # On Linux/macOS

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Running a Simulation

Use the `solver.py` script with a YAML configuration file.

```bash
# Run with default configuration
python solver.py

# Run with a specific configuration
python solver.py config/my_simulation.yaml
```

### Configuration

Simulations are controlled via `.yaml` files in the `config/` directory. You can specify:
*   **Simulation type:** `euler_2d` or `navier_stokes_2d`.
*   **Mesh:** Resolution and bounds.
*   **Numerical:** Polynomial degree for DG.
*   **Device:** `cpu` or `cuda`.

### Visualization

Results are written to the `output_dir` specified in the config (default: `data/`).

1.  Open the `.xmf` file in **ParaView**.

## Development & Testing

### Running Tests

We use `pytest` for unit and verification testing.

```bash
# Run all tests
python -m pytest

# Run specific verification cases
python tests/verification/vortex_2D/verify.py
```

### Verification Cases

The `tests/verification/` directory contains several standard cases:
*   **Isentropic Vortex:** Order of accuracy verification.
    *   [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dzilles/pazuzu/blob/feature/quadtree-amr-fr/tests/verification/vortex_2D/vortex_analysis.ipynb)
*   **Channel Flow:** Wall boundary conditions.

*   **Acoustic Pulse:** Wave propagation.
*   **Sod Shock Tube:** Discontinuity handling.

