from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, List, Union, Literal
import yaml
from pydantic import BaseModel, Field, ValidationError

class FluxType(str, Enum):
    RUSANOV = "rusanov"
    HLLC = "hllc"
    ROE = "roe"

class LimiterType(str, Enum):
    NONE = "none"
    BARTH_JESPERSEN = "barth_jespersen"
    MINMOD = "minmod"

class PhysicsConfig(BaseModel):
    gamma: float = Field(default=1.4, description="Adiabatic index (heat capacity ratio)")
    gas_constant: float = Field(default=287.0, description="Specific gas constant (R)")

class NumericsConfig(BaseModel):
    cfl: float = Field(default=0.4, gt=0, description="CFL number (must be > 0)")
    polynomial_order: int = Field(default=1, ge=0, description="Polynomial degree for DG")
    precision: Literal["single", "double"] = Field(default="single", description="Floating point precision (single or double)")
    flux_type: FluxType = Field(default=FluxType.HLLC, description="Numerical flux scheme")
    limiter: LimiterType = Field(default=LimiterType.NONE, description="Limiter for shock capturing")
    use_filtering: bool = Field(default=False, description="Enable exponential filtering against aliasing")
    filter_alpha: float = Field(default=36.0, description="Filter strength (alpha)")
    filter_order: int = Field(default=16, description="Filter order")

class IOConfig(BaseModel):
    output_dir: str = Field(default="output", description="Directory for output files")
    write_interval: float = Field(..., gt=0, description="Simulation time interval between outputs")

class SimulationConfig(BaseModel):
    t_final: float = Field(default=1.0, gt=0, description="Final simulation time")
    ramp_time: float = Field(default=1.0, description="Ramp time for boundary conditions")
    device: str = Field(default="cpu", description="Compute device (cpu or cuda)")
    max_steps: Optional[int] = Field(default=None, description="Maximum number of time steps")

class PazuzuConfig(BaseModel):
    case_name: str = Field(default="simulation", description="Name of the simulation case")
    mesh_file: str = Field(..., description="Path to the mesh file (.msh)")
    initial_condition: str = Field(default="vortex", description="Initial condition name")
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
