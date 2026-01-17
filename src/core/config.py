import os
from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, Literal
import yaml
from pydantic import BaseModel, Field, ValidationError, ConfigDict, model_validator

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

class ShockCapturingConfig(BaseModel):
    enabled: bool = Field(default=False, description="Enable hybrid shock capturing (FR/FV)")
    indicator_threshold: float = Field(default=0.001, gt=0.0, description="Energy threshold (Se) for marking troubled cells")
    indicator_alpha: float = Field(default=36.0, gt=0.0, description="Spectral filter strength (alpha)")
    indicator_order: int = Field(default=4, ge=1, description="Spectral filter order")

class NumericsConfig(BaseModel):
    cfl: float = Field(default=0.4, gt=0, description="CFL number (must be > 0)")
    polynomial_order: int = Field(default=1, ge=0, description="Polynomial degree for DG")
    precision: Literal["single", "double"] = Field(default="single", description="Floating point precision (single or double)")
    flux: Literal["rusanov", "hllc"] = Field(default="rusanov", description="Numerical flux function")
    hllc_fallback: bool = Field(default=True, description="Fallback to Rusanov flux if HLLC fails")
    shock_capturing: ShockCapturingConfig = Field(default_factory=ShockCapturingConfig, description="Shock capturing settings")
    time_integrator: Literal["ssp_rk3", "rk4"] = Field(default="ssp_rk3", description="Time integration scheme")
    dt_min: float = Field(default=1e-9, gt=0, description="Minimum allowed timestep")
    dt_init: float = Field(default=0.001, gt=0, description="Initial timestep")
    dt_static: Optional[float] = Field(default=None, gt=0, description="Optional fixed timestep to override CFL-based DT")

class IOConfig(BaseModel):
    output_dir: Optional[str] = Field(default=None, description="Directory for output files")
    write_interval: int = Field(..., gt=0, description="Number of simulation steps between outputs")

class SimulationConfig(BaseModel):
    t_final: float = Field(default=1.0, gt=0, description="Final simulation time")
    max_steps: int = Field(default=1000000, ge=1, description="Maximum number of steps")
    ramp_up_time: float = Field(default=0.0, ge=0.0, description="Time period to ramp up inlet velocities")
    device: str = Field(default="automatic", description="Compute device (cpu, cuda, or automatic)")

class InitialConditionConfig(BaseModel):
    name: str = Field(..., description="Name of the initial condition")
    params: Dict[str, Any] = Field(default_factory=dict, description="Parameters for the initial condition")

class AmrConfig(BaseModel):
    enabled: bool = Field(default=False, description="Enable Adaptive Mesh Refinement")
    max_blocks: int = Field(default=10000, gt=0, description="Maximum number of blocks in the memory pool")
    max_depth: int = Field(default=10, ge=0, description="Maximum allowed refinement depth")
    refinement_threshold: Optional[float] = Field(default=None, description="Gradient threshold for adaptive refinement")
    coarsening_threshold: Optional[float] = Field(default=None, description="Gradient threshold for adaptive coarsening")
    refine_interval: int = Field(default=0, ge=0, description="Steps between AMR updates (0 to disable dynamic AMR)")

class MeshConfig(BaseModel):
    initial_depth: int = Field(default=3, ge=0, description="Initial refinement depth (uniform)")
    x_min: float = Field(default=-1.0, description="Domain x-min")
    x_max: float = Field(default=1.0, description="Domain x-max")
    y_min: float = Field(default=-1.0, description="Domain y-min")
    y_max: float = Field(default=1.0, description="Domain y-max")
    periodic_x: bool = Field(default=False, description="Periodic boundary in X")
    periodic_y: bool = Field(default=False, description="Periodic boundary in Y")

    @model_validator(mode='after')
    def check_bounds(self):
        if self.x_max <= self.x_min:
            raise ValueError(f"x_max ({self.x_max}) must be greater than x_min ({self.x_min})")
        if self.y_max <= self.y_min:
            raise ValueError(f"y_max ({self.y_max}) must be greater than y_min ({self.y_min})")
        return self

class IBMConfig(BaseModel):
    enabled: bool = Field(default=False, description="Enable Immersed Boundary Method")
    mode: Literal["analytical", "stl_file"] = Field(default="analytical", description="Geometry source mode")
    boundary_type: Literal["slip", "no_slip"] = Field(default="slip", description="Boundary condition type for immersed surfaces")
    stl_path: Optional[str] = Field(default=None, description="Path to the .stl file (if mode is stl_file)")
    invert_inside_outside: bool = Field(default=False, description="Invert the signed distance field (useful if STL normals are flipped)")
    geometric_params: Dict[str, Any] = Field(default_factory=dict, description="Parameters for analytical shapes (e.g., radius, center)")

class PazuzuConfig(BaseModel):
    case_name: str = Field(default="simulation", description="Name of the simulation case")
    solver_type: SolverType = Field(default=SolverType.EULER_2D, description="Type of solver to use")
    mesh: MeshConfig = Field(default_factory=MeshConfig, description="Mesh/Domain configuration")
    amr: AmrConfig = Field(default_factory=AmrConfig, description="AMR configuration")
    ibm: IBMConfig = Field(default_factory=IBMConfig, description="Immersed Boundary Method configuration")
    initial_condition: InitialConditionConfig = Field(default_factory=lambda: InitialConditionConfig(name="vortex"), description="Initial condition configuration")
    physics: PhysicsConfig = Field(default_factory=PhysicsConfig)
    numerics: NumericsConfig = Field(default_factory=NumericsConfig)
    io: IOConfig = Field(..., description="I/O configuration")
    simulation: SimulationConfig = Field(default_factory=SimulationConfig)
    boundaries: Dict[str, Dict[str, Any]] = Field(default_factory=dict, description="Boundary conditions map")

    model_config = ConfigDict(frozen=True)

    @model_validator(mode='after')
    def check_amr_depths(self):
        if self.amr.enabled and self.amr.max_depth < self.mesh.initial_depth:
            raise ValueError(f"amr.max_depth ({self.amr.max_depth}) must be >= mesh.initial_depth ({self.mesh.initial_depth})")
        return self

    @model_validator(mode='before')
    @classmethod
    def normalize_ic(cls, data: Any) -> Any:
        if isinstance(data, dict):
             if 'initial_condition' in data:
                 ic = data['initial_condition']
                 if isinstance(ic, str):
                     data['initial_condition'] = {'name': ic, 'params': {}}
        return data

    def override(self, overrides: Dict[str, Any]) -> "PazuzuConfig":
        """
        Returns a new PazuzuConfig instance with updated values.
        Support dot-notation for nested fields, e.g., {'numerics.cfl': 0.1}.
        Raises KeyError if a key does not exist in the original config.
        """
        data = self.model_dump()
        
        for path, value in overrides.items():
            keys = path.split('.')
            d = data
            for key in keys[:-1]:
                if key not in d:
                    raise KeyError(f"Invalid config path: '{path}'. Key '{key}' not found.")
                d = d[key]
            
            if keys[-1] not in d:
                raise KeyError(f"Invalid config path: '{path}'. Key '{keys[-1]}' not found.")
                
            d[keys[-1]] = value
            
        return PazuzuConfig(**data)

    @classmethod
    def from_yaml(cls, path: str, project_root: Optional[str] = None) -> "PazuzuConfig":
        path_obj = Path(path)
        if not path_obj.exists():
            raise FileNotFoundError(f"Configuration file not found: {path}")
        
        with open(path_obj, 'r') as f:
            try:
                data = yaml.safe_load(f)
            except yaml.YAMLError as e:
                raise ValueError(f"Error parsing YAML file: {e}")
        
        # --- [NEW] Centralized Output Logic ---
        # If 'output_dir' is NOT specified in YAML, enforce the project-relative structure.
        if project_root:
            io_data = data.get('io', {})
            if 'output_dir' not in io_data:
                # Default: pazuzu/output/<case_name>
                case_name = data.get('case_name', 'simulation')
                default_dir = os.path.join(project_root, "output", case_name)
                
                # Inject back into data dict before validation
                if 'io' not in data:
                    data['io'] = {}
                data['io']['output_dir'] = default_dir

            # IBM STL Path Resolution
            ibm_data = data.get('ibm', {})
            if ibm_data.get('mode') == 'stl_file' and 'stl_path' in ibm_data:
                stl_path = ibm_data['stl_path']
                if stl_path and not os.path.isabs(stl_path):
                    # Resolve relative to project root
                    if 'ibm' not in data: # Should exist if we are here
                        data['ibm'] = {}
                    data['ibm']['stl_path'] = os.path.join(project_root, stl_path)
        # --------------------------------------

        try:
            # If output_dir was still None (no project_root and no yaml key), 
            # fall back to "output" in current dir
            if 'io' in data and data['io'].get('output_dir') is None:
                data['io']['output_dir'] = "output"
                
            return cls(**data)
        except ValidationError as e:
            raise ValueError(f"Configuration validation failed:\n{e}")