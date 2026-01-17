# Karman Vortex Street (Euler)

This test case simulates inviscid flow over a circular cylinder.

## Configuration
- **Domain**: `[-5, 15] x [-5, 5]`
- **Cylinder**: Radius `0.5` at `(0, 0)`
- **Mach**: ~0.85 (u=1, c=sqrt(1.4*1/1) ~ 1.18) -> Subsonic/Transonic.
- **AMR**: Geometry-aware refinement near cylinder.

## Running
To run the full simulation:
```bash
python solver.py tests/verification/karman_vortex_street/karman_euler.yaml
```

## Notes
- Since this is an Euler simulation (inviscid) with slip walls, physical vortex shedding (Karman street) theoretically should not occur due to d'Alembert's paradox (potential flow solution).
- However, numerical viscosity from the Flux Reconstruction scheme (Rusanov/HLLC) and the "staircase" effect of the Immersed Boundary Method (even with ghost cells) usually introduce enough perturbation and separation to trigger shedding.
- Boundary Conditions are currently Transmissive (Zero Gradient) at domain edges. Ensure the domain is large enough to avoid reflection artifacts.
