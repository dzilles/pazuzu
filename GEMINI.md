# Pazuzu: Quadtree-AMR Flux Reconstruction Solver

**Branch:** `feature/quadtree-amr-fr`

This document outlines the architectural rewrite of Pazuzu from a generic unstructured DG solver to a specialized **Cartesian Quadtree** solver using **Flux Reconstruction (FR)** and **Block-Adaptive Mesh Refinement (AMR)**.

---

## 1. Core Architecture

### **A. Mesh Topology: Cartesian Quadtree**
* **Structure:** Hierarchical Quadtree using **Morton Encoding (Z-ordering)** for fast spatial indexing on GPU.
* **Connectivity:** Implicit connectivity (parent/child relationships) replaces explicit neighbor lists.
* **AMR Strategy:** Block-based refinement. A "Memory Pool" of fixed-size blocks (e.g., $8 \times 8$ or $16 \times 16$ cells) is pre-allocated on the GPU. Active leaf nodes map to slots in this pool.

### **B. Numerical Method: Flux Reconstruction (FR)**
* **Formulation:** Differential Form ($\frac{\partial f}{\partial x} + g_{corr}$). 
* **Basis:** Gauss-Lobatto-Legendre (GLL) nodal basis (Tensor Product).
* **Correction Functions:** Radau or Legendre polynomials (e.g., Huynh's $g_2$) to recover high-order accuracy.
* **Shock Capturing:** A hybrid **Sub-Cell Finite Volume** method. "Troubled" blocks switch from high-order FR to a robust FV update on the sub-grid nodes.

### **C. Geometry: Immersed Boundary Method (IBM)**
* **Representation:** Signed Distance Fields (SDF).
* **Boundary Handling:** Ghost-cell forcing or cut-cell integration (depending on implementation phase). No body-fitted meshing required.

---

## 2. Project Structure

The directory structure has been refactored to support the new paradigm.

```text
├── config/                 # YAML configs (now including AMR/IBM settings)
├── src/
│   ├── core/
│   │   ├── basis.py        # GLL nodes, weights, & Correction Polynomials
│   │   ├── simulation_state.py # Memory Pool (MAX_BLOCKS allocation)
│   │   └── config.py       # Config dataclasses
│   ├── geometry/
│   │   ├── quadtree.py     # Morton encoding, Refinement/Coarsening logic
│   │   └── ibm.py          # Level-set / SDF management
│   ├── physics/
│   │   └── laws/           # Pure physics (Euler/NS Flux functions) - UNCHANGED
│   ├── numerics/
│   │   ├── flux_reconstruction.py # 1D Correction function derivatives
│   │   └── time_steppers.py       # SSP-RK3 / RK4
│   ├── kernels/
│   │   ├── fr_kernels.py   # Differential form FR update
│   │   ├── amr_kernels.py  # Inter-block interpolation (Prolongation/Restriction)
│   │   ├── fv_kernels.py   # Sub-cell Finite Volume update
│   │   └── boundary_conditions.py # Standard Riemann/Slip BCs
│   └── io/                 # HDF5 writers (Updated for Block-Structured data)
└── tests/                  # Unit and verification tests
```

## 3. Implementation Phases

### Phase 1: The Foundation (Quadtree & Memory)
*   **Goal:** Establish grid data structures without physics.
*   **Key Tasks:**
    *   [x] Implement Quadtree class (Morton codes, uniform refinement).
    *   [x] Implement SimulationState with Memory Pool allocation (MAX_BLOCKS).
    *   [x] Create mapping logic: Leaf Node -> Pool Index.

### Phase 2: Flux Reconstruction (Static Grid)
*   **Goal:** Run physics on a uniform, non-adaptive grid.
*   **Key Tasks:**
    *   [x] Implement 1D Correction Polynomial derivatives in basis.py.
    *   [x] Write fr_kernels.py (Differential form: compute flux gradients + add correction).
    *   [x] Verify order of accuracy (Isentropic Vortex).

### Phase 3: Adaptive Mesh Refinement (AMR)
*   **Goal:** Enable dynamic hanging-node refinement.
*   **Key Tasks:**
    *   [x] Implement P-Multigrid kernels (Interpolation between Coarse/Fine blocks).
    *   [x] Handle Mortar Interfaces (2:1 balance) in flux kernels.
    *   [x] Implement the Mark -> Refine -> Balance loop.

### Phase 4: Robustness (Sub-Cell FV)
*   **Goal:** Shock capturing.
*   **Key Tasks:**
    *   [x] Implement "Troubled Cell" detectors (Persson-Peraire).
    *   [x] Write fv_kernels.py (First-order Godunov on sub-grid).
    *   [x] Test on Sod Shock Tube and Double Mach Reflection.

### Phase 5: Immersed Boundary Method (IBM)
*   **Goal:** Complex geometry support.
*   **Key Tasks:**
    *   [ ] Load/Generate Signed Distance Fields (SDF).
    *   [ ] Implement Ghost-Node forcing kernels.

## 4. Developer Guidelines
*   **Environment:** Always use the virtual environment located at `.venv`.
    *   **Activate (Windows):** `.venv\Scripts\activate`
    *   **Run Tests:** `.venv\Scripts\pytest` or `python -m pytest` (after activation).
    *   **Run Linting:** `.venv\Scripts\ruff check src`
    *   **Run Type Checking:** `.venv\Scripts\mypy src`
*   **Memory Management:** Never allocate arrays per-element. Always allocate (MAX_BLOCKS, ...) at startup.
*   **Kernel Design:** Kernels should launch over `num_active_blocks`. Use the `active_block_indices` map to access the global memory pool.
*   **Testing:** Run pytest frequently. Each phase is designed to be independently verifiable.