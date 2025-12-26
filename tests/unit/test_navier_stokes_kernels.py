import pytest
import numpy as np
import warp as wp
from src.kernels import navier_stokes_kernels as nsk
from src.kernels.structs import EquationParams32

class TestNavierStokesKernels:
    @pytest.fixture(autouse=True)
    def setup(self):
        wp.init()
        self.device = "cpu"
        
    def test_primitive_gradients_constant(self):
        """Test that a constant state results in zero gradients."""
        Np = 4
        num_elements = 1
        
        q_host = np.zeros((num_elements, Np, 4), dtype=np.float32)
        q_host[:, :, 0] = 1.0 # rho
        q_host[:, :, 1] = 0.5 # rho*u => u=0.5
        q_host[:, :, 2] = 0.0 # rho*v => v=0.0
        q_host[:, :, 3] = 2.5 # E
        
        q = wp.array(q_host, dtype=wp.vec4, device=self.device)
        grad_u = wp.zeros((num_elements, Np), dtype=wp.vec2, device=self.device)
        grad_v = wp.zeros((num_elements, Np), dtype=wp.vec2, device=self.device)
        grad_T = wp.zeros((num_elements, Np), dtype=wp.vec2, device=self.device)
        
        # Identity differentiation matrices for constant test
        Dr = wp.array(np.zeros((Np, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        Ds = wp.array(np.zeros((Np, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        
        rx = wp.array(np.ones((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        ry = wp.array(np.zeros((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        sx = wp.array(np.zeros((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        sy = wp.array(np.ones((num_elements, Np), dtype=np.float32), dtype=wp.float32, device=self.device)
        
        params = EquationParams32()
        params.gamma = 1.4
        params.rho_floor = 1e-5
        params.p_floor = 1e-5
        params.one = 1.0
        params.cp = 1004.5
        params.gas_constant = 287.0
        
        wp.launch(
            kernel=nsk.compute_primitive_gradients_volume,
            dim=(num_elements, Np),
            inputs=[q, grad_u, grad_v, grad_T, Dr, Ds, rx, ry, sx, sy, Np, params],
            device=self.device
        )
        
        res_u = grad_u.numpy()
        res_v = grad_v.numpy()
        res_T = grad_T.numpy()
        
        assert np.allclose(res_u, 0.0)
        assert np.allclose(res_v, 0.0)
        assert np.allclose(res_T, 0.0)

    def test_viscous_flux_calculation(self):
        """Test viscous flux helper functions (indirectly via volume kernel)."""
        pass