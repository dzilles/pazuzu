from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, List, Union, Literal
import yaml
from pydantic import BaseModel, Field, ValidationError

class FluxType(str, Enum):
    RUSANOV = "rusanov"
    HLLC = "hllc"
    ROE = "roe"

class SolverType(str, Enum):
    EULER_2D = "euler_2d"
    NAVIER_STOKES_2D = "navier_stokes_2d"

class LimiterType(str, Enum):
    NONE = "none"
    BARTH_JESPERSEN = "barth_jespersen"
    MINMOD = "minmod"

class PhysicsConfig(BaseModel):
    gamma: float = Field(default=1.4, description="Adiabatic index (heat capacity ratio)")
    gas_constant: float = Field(default=287.0, description="Specific gas constant (R)")
    # Freestream / Background state
    rho_inf: float = Field(default=1.0, description="Freestream density")
    u_inf: float = Field(default=0.0, description="Freestream x-velocity")
    v_inf: float = Field(default=0.0, description="Freestream y-velocity")
    p_inf: float = Field(default=1.0, description="Freestream pressure")

class NumericsConfig(BaseModel):
    cfl: float = Field(default=0.4, gt=0, description="CFL number (must be > 0)")
    polynomial_order: int = Field(default=1, ge=0, description="Polynomial degree for DG")
    precision: Literal["single", "double"] = Field(default="single", description="Floating point precision (single or double)")
    flux_type: FluxType = Field(default=FluxType.HLLC, description="Numerical flux scheme")
    limiter: LimiterType = Field(default=LimiterType.NONE, description="Limiter for shock capturing")
    use_filtering: bool = Field(default=False, description="Enable exponential filtering against aliasing")
    filter_alpha: float = Field(default=36.0, description="Filter strength (alpha)")
    filter_order: int = Field(default=16, description="Filter order")
    min_dt: float = Field(default=1e-15, gt=0, description="Minimum allowable time step before aborting")

class IOConfig(BaseModel):
    output_dir: str = Field(default="output", description="Directory for output files")
    write_interval: float = Field(..., gt=0, description="Simulation time interval between outputs")

class SimulationConfig(BaseModel):
    t_final: float = Field(default=1.0, gt=0, description="Final simulation time")
    ramp_time: float = Field(default=1.0, description="Ramp time for boundary conditions")
    device: str = Field(default="cpu", description="Compute device (cpu or cuda)")
    max_steps: Optional[int] = Field(default=None, description="Maximum number of time steps")

class InitialConditionConfig(BaseModel):
    name: str = Field(..., description="Name of the initial condition")
    params: Dict[str, Any] = Field(default_factory=dict, description="Parameters for the initial condition")

class PazuzuConfig(BaseModel):
    case_name: str = Field(default="simulation", description="Name of the simulation case")
    solver_type: SolverType = Field(default=SolverType.EULER_2D, description="Type of solver to use")
    mesh_file: str = Field(..., description="Path to the mesh file (.msh)")
    initial_condition: Union[str, InitialConditionConfig] = Field(default="vortex", description="Initial condition configuration")
    periodic_pairs: List[List[Union[str, int]]] = Field(default_factory=list, description="Periodic boundary pairs [Tag1, Tag2, Axis]")
    physics: PhysicsConfig = Field(default_factory=PhysicsConfig)
    numerics: NumericsConfig = Field(default_factory=NumericsConfig)
    io: IOConfig = Field(..., description="I/O configuration")
    simulation: SimulationConfig = Field(default_factory=SimulationConfig)
    boundaries: Dict[str, Dict[str, Any]] = Field(default_factory=dict, description="Boundary conditions map")

    class Config:
        frozen = True

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
