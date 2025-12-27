from src.core.config import PazuzuConfig, SolverType
from src.physics.solvers.euler_2d import Euler2DSolver
from src.physics.solvers.navier_stokes_2d import NavierStokes2DSolver

class SolverBuilder:
    """
    Builder class for creating solver instances.
    Decouples the main driver from specific physics implementations.
    """
    def __init__(self, mesh, basis, config: PazuzuConfig):
        """
        Initializes the SolverBuilder.

        Args:
            mesh (Mesh): The computational mesh.
            basis (Basis): The DG basis.
            config (PazuzuConfig): Configuration object.
        """
        self.mesh = mesh
        self.basis = basis
        self.config = config

    def build(self):
        """
        Builds and returns the appropriate solver instance based on the configuration.

        Returns:
            BaseSolver: An instance of a class inheriting from BaseSolver.
        """
        if self.config.solver_type == SolverType.EULER_2D:
            return Euler2DSolver(self.mesh, self.basis, config=self.config)
        elif self.config.solver_type == SolverType.NAVIER_STOKES_2D:
            return NavierStokes2DSolver(self.mesh, self.basis, config=self.config)
        else:
            raise ValueError(f"Unknown solver type: {self.config.solver_type}")
