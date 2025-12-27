import pytest
import numpy as np
from src.core.config import PhysicsConfig, PazuzuConfig
from src.core.boundary_condition_manager import BoundaryConditionManager
from src.kernels import boundary_conditions as bc
from pydantic import ValidationError

def test_physics_config_validation():
    # Valid config
    PhysicsConfig(gamma=1.4, gas_constant=287.0, cp=1004.5, rho_inf=1.0, p_inf=1.0)
    
    # Invalid gamma
    with pytest.raises(ValidationError):
        PhysicsConfig(gamma=0.9)
    
    # Invalid rho_inf
    with pytest.raises(ValidationError):
        PhysicsConfig(rho_inf=-1.0)
        
    # Invalid p_inf
    with pytest.raises(ValidationError):
        PhysicsConfig(p_inf=0.0)

def test_bc_manager_validation():
    class MockMesh:
        def __init__(self):
            self.num_elements = 1
            self.boundary_tags_host = np.zeros((1, 4), dtype=np.int32)
            self.physical_groups = {"Inlet": 1}
            
    class MockConfig:
        def __init__(self, boundaries):
            self.boundaries = boundaries
            self.physics = PhysicsConfig()
            
    mesh = MockMesh()
    
    # Valid inlet
    boundaries = {"Inlet": {"type": "inlet", "params": {"rho": 1.0, "u": 1.0, "v": 0.0, "p": 1.0}}}
    manager = BoundaryConditionManager(mesh, MockConfig(boundaries))
    manager.setup_boundary_conditions()
    
    # Invalid rho in inlet
    boundaries = {"Inlet": {"type": "inlet", "params": {"rho": -1.0, "u": 1.0, "v": 0.0, "p": 1.0}}}
    manager = BoundaryConditionManager(mesh, MockConfig(boundaries))
    with pytest.raises(ValueError, match="invalid rho"):
        manager.setup_boundary_conditions()
        
    # Missing parameter
    boundaries = {"Inlet": {"type": "inlet", "params": {"u": 1.0, "v": 0.0, "p": 1.0}}}
    manager = BoundaryConditionManager(mesh, MockConfig(boundaries))
    with pytest.raises(ValueError, match="missing required parameter"):
        manager.setup_boundary_conditions()

if __name__ == "__main__":
    pytest.main([__file__])
