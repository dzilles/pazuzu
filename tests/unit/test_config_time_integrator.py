
import pytest
from src.core.config import PazuzuConfig, NumericsConfig

def test_default_time_integrator():
    config_dict = {
        "case_name": "test",
        "solver_type": "euler_2d",
        "io": {"write_interval": 10}
    }
    config = PazuzuConfig(**config_dict)
    assert config.numerics.time_integrator == "ssp_rk3"

def test_rk4_time_integrator():
    config_dict = {
        "case_name": "test",
        "solver_type": "euler_2d",
        "numerics": {"time_integrator": "rk4"},
        "io": {"write_interval": 10}
    }
    config = PazuzuConfig(**config_dict)
    assert config.numerics.time_integrator == "rk4"

def test_invalid_time_integrator():
    config_dict = {
        "case_name": "test",
        "solver_type": "euler_2d",
        "numerics": {"time_integrator": "invalid_scheme"},
        "io": {"write_interval": 10}
    }
    with pytest.raises(Exception): # Pydantic validation error
        PazuzuConfig(**config_dict)
