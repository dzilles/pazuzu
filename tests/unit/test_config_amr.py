import pytest
from src.core.config import PazuzuConfig, AmrConfig

def test_amr_config_defaults():
    """Verify default values for AMR configuration."""
    conf = AmrConfig()
    assert conf.max_blocks == 10000
    assert conf.initial_depth == 3

def test_amr_config_validation():
    """Verify validation logic for AMR configuration."""
    # max_blocks must be > 0
    with pytest.raises(Exception): # Pydantic raises ValidationError
        AmrConfig(max_blocks=0)

def test_pazuzu_config_integration():
    """Verify PazuzuConfig correctly incorporates AmrConfig."""
    # Create a minimal valid config dictionary
    config_dict = {
        "io": {"write_interval": 10},
        "initial_condition": {"name": "test"},
        "amr": {
            "max_blocks": 500,
            "initial_depth": 5
        }
    }
    
    conf = PazuzuConfig(**config_dict)
    assert conf.amr.max_blocks == 500
    assert conf.amr.initial_depth == 5