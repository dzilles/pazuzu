# AGENT PROTOCOL & CONTEXT

> **SYSTEM INSTRUCTION:** This file contains the absolute truth about the project architecture, coding rules, and current status. Read this before modifying any code.

## 1. PROJECT OVERVIEW & STATUS

**Goal:** An explicit Discontinuous Galerkin (DG) solver for 2D Euler equations using Nvidia Warp (GPU).

**Current Architecture:**
- **Method:** Nodal DG (Strong Form).
- **Time Integration:** Low-Storage SSP-RK3.
- **Flux:** Rusanov (Local Lax-Friedrichs).
- **Mesh:** Currently Cartesian (logic being refactored for unstructured/Gmsh).

**Roadmap (The "North Star"):**
1.  **Modularization:** Split monolithic scripts into Topology, Physics, State, Operator.
2.  **I/O:** Move from `.npz` to **HDF5/XDMF** or **VTK**.
3.  **Meshing:** Implement **Gmsh** import (.msh) for complex geometries.
4.  **Physics:** Implement reactive flow (Combustion) and Moving Mesh (Rotors).

## 2. CODING RULES & CONSTRAINTS (CRITICAL)

### A. Nvidia Warp Specifics
1.  **Kernel Purity:** Code inside `@wp.kernel` or `@wp.func` MUST NOT use Python lists, dicts, or dynamic allocation.
2.  **Types:** Strictly use `wp.vec2`, `wp.vec4`, `wp.float32`. **NEVER** use `np.float64` inside kernels.
3.  **Loops:** Use `for i in range(Start, End):`. No `enumerate` or `zip`.
4.  **Atomic Ops:** Use `wp.atomic_add()` for flux accumulation.

### B. Architecture & Style
1.  **Host/Device Separation:** Python host code and Warp device kernels must be strictly separated.
2.  **Memory:** Heavy data resides in `wp.array` on `device="cuda"`. Only copy to CPU (`.numpy()`) for I/O.
3.  **Interfaces:** Use `PhysicsModel` abstraction for new equations. No hardcoded `if/else` in solver loops.
4.  **Docstrings:** Every new function requires a docstring explaining inputs/outputs.

## 3.

## 4. PROJECT STRUCTURE (Reference)

Pazuzu/
├── config/                 # YAML configs (Planned)
├── data/                   # Output (.npz, .h5)
├── meshes/                 # Gmsh files (.msh)
├── scripts/                # helper scripts
├── src/
│   ├── core/               # Solver logic (solver.py)
│   ├── geometry/           # Mesh handling (mesh.py)
│   ├── kernels/            # GPU Kernels (warp_kernels.py)
│   └── physics/            # Equations & Fluxes (equations.py)
├── run_simulation.py       # Entry point
└── Agents.md               # THIS FILE