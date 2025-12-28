import pytest
from src.core.config import PazuzuConfig, AmrConfig

def test_amr_config_defaults():
    """Verify default values for AMR configuration."""
    conf = AmrConfig()
    assert conf.max_blocks == 10000
    assert conf.max_depth == 5
    assert conf.refinement_threshold == 0.1

def test_amr_config_validation():
    """Verify validation logic for AMR configuration."""
    # max_blocks must be > 0
    with pytest.raises(Exception): # Pydantic raises ValidationError
        AmrConfig(max_blocks=0)
    
    # max_depth must be >= 0
    with pytest.raises(Exception):
        AmrConfig(max_depth=-1)

def test_pazuzu_config_integration():
    """Verify PazuzuConfig correctly incorporates AmrConfig."""
    # Create a minimal valid config dictionary
    config_dict = {
        "mesh_file": "dummy.msh", # Kept in config definition for now, though unused by logic
        "io": {"write_interval": 10},
        "initial_condition": {"name": "test"},
        "amr": {
            "max_blocks": 500,
            "max_depth": 3
        }
    }
    
    # We might have removed mesh_file in the previous turns? 
    # Let's check the replace call in history. 
    # Yes, I removed mesh_file in src/core/config.py.
    
    config_dict_no_mesh = {
        "io": {"write_interval": 10},
        "initial_condition": {"name": "test"},
        "amr": {
            "max_blocks": 500,
            "max_depth": 3
        }
    }
    
    conf = PazuzuConfig(**config_dict_no_mesh)
    assert conf.amr.max_blocks == 500
    assert conf.amr.max_depth == 3
