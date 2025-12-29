from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, List, Union, Literal
import yaml
from pydantic import BaseModel, Field, ValidationError, ConfigDict

class SolverType(str, Enum):
    EULER_2D = "euler_2d"
    NAVIER_STOKES_2D = "navier_stokes_2d"

class PhysicsConfig(BaseModel):
    gamma: float = Field(default=1.4, gt=1.0, description="Adiabatic index (heat capacity ratio)")
    gas_constant: float = Field(default=287.0, gt=0.0, description="Specific gas constant (R)")
    # Freestream / Background state
    rho_inf: float = Field(default=1.0, gt=0.0, description="Freestream density")
    u_inf: float = Field(default=0.0, description="Freestream x-velocity")
    v_inf: float = Field(default=0.0, description="Freestream y-velocity")
    p_inf: float = Field(default=1.0, gt=0.0, description="Freestream pressure")
    rho_floor: float = Field(default=1.0e-5, gt=0.0, description="Minimum allowed density")
    p_floor: float = Field(default=1.0e-5, gt=0.0, description="Minimum allowed pressure")

class NumericsConfig(BaseModel):
    cfl: float = Field(default=0.4, gt=0, description="CFL number (must be > 0)")
    polynomial_order: int = Field(default=1, ge=0, description="Polynomial degree for DG")
    precision: Literal["single", "double"] = Field(default="single", description="Floating point precision (single or double)")

class IOConfig(BaseModel):
    output_dir: str = Field(default="output", description="Directory for output files")
    write_interval: int = Field(..., gt=0, description="Number of simulation steps between outputs")

class SimulationConfig(BaseModel):
    t_final: float = Field(default=1.0, gt=0, description="Final simulation time")
    device: str = Field(default="cpu", description="Compute device (cpu or cuda)")

class InitialConditionConfig(BaseModel):
    name: str = Field(..., description="Name of the initial condition")
    params: Dict[str, Any] = Field(default_factory=dict, description="Parameters for the initial condition")

class AmrConfig(BaseModel):
    max_blocks: int = Field(default=10000, gt=0, description="Maximum number of blocks in the memory pool")
    initial_depth: int = Field(default=3, ge=0, description="Initial refinement depth (uniform)")

class MeshConfig(BaseModel):
    x_min: float = Field(default=-1.0, description="Domain x-min")
    x_max: float = Field(default=1.0, description="Domain x-max")
    y_min: float = Field(default=-1.0, description="Domain y-min")
    y_max: float = Field(default=1.0, description="Domain y-max")
    periodic_x: bool = Field(default=False, description="Periodic boundary in X")
    periodic_y: bool = Field(default=False, description="Periodic boundary in Y")

class PazuzuConfig(BaseModel):
    case_name: str = Field(default="simulation", description="Name of the simulation case")
    solver_type: SolverType = Field(default=SolverType.EULER_2D, description="Type of solver to use")
    mesh: MeshConfig = Field(default_factory=MeshConfig, description="Mesh/Domain configuration")
    amr: AmrConfig = Field(default_factory=AmrConfig, description="AMR configuration")
    initial_condition: Union[str, InitialConditionConfig] = Field(default="vortex", description="Initial condition configuration")
    physics: PhysicsConfig = Field(default_factory=PhysicsConfig)
    numerics: NumericsConfig = Field(default_factory=NumericsConfig)
    io: IOConfig = Field(..., description="I/O configuration")
    simulation: SimulationConfig = Field(default_factory=SimulationConfig)
    boundaries: Dict[str, Dict[str, Any]] = Field(default_factory=dict, description="Boundary conditions map")

    model_config = ConfigDict(frozen=True)

    @classmethod
    def from_yaml(cls, path: str) -> "PazuzuConfig":
        path_obj = Path(path)
        if not path_obj.exists():
            raise FileNotFoundError(f"Configuration file not found: {path}")
        
        with open(path_obj, 'r') as f:
            try:
                data = yaml.safe_load(f)
            except yaml.YAMLError as e:
                raise ValueError(f"Error parsing YAML file: {e}")
        
        try:
            return cls(**data)
        except ValidationError as e:
            raise ValueError(f"Configuration validation failed:\n{e}")