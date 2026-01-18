import pytest
import numpy as np
import warp as wp
import sys
import os

# Add project root (parent of src) to path to allow 'src.' imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.kernels import common_kernels
from src.kernels.structs import EquationParams32

class TestCommonKernels:
    
        def test_rk_stage_1(self, device):
    
            wp.init()
    
            num_elems = 2
    
            Np = 4
    
            q = np.ones((num_elems, Np, 4), dtype=np.float32)
    
            rhs = np.full((num_elems, Np, 4), 0.5, dtype=np.float32)
    
            dt = 0.1
    
            active_indices = np.arange(num_elems, dtype=np.int32)
            params = EquationParams32()
    
            q_wp = wp.array(q, dtype=wp.vec4, device=device)
    
            rhs_wp = wp.array(rhs, dtype=wp.vec4, device=device)
    
            q_out_wp = wp.zeros_like(q_wp)
    
            active_wp = wp.array(active_indices, dtype=wp.int32, device=device)
    
            
    
            wp.launch(
    
                kernel=common_kernels.rk_stage_1,
    
                dim=(num_elems, Np),
    
                inputs=[q_wp, rhs_wp, dt, q_out_wp, active_wp, params],
    
                device=device
    
            )
    
            
    
            expected = q + dt * rhs
    
            np.testing.assert_allclose(q_out_wp.numpy(), expected, atol=1e-6)
    
    
    
        def test_rk_stage_2(self, device):
    
            wp.init()
    
            num_elems = 1
    
            Np = 1
    
            q = np.array([[[1.0, 1.0, 1.0, 1.0]]], dtype=np.float32)
    
            q_1 = np.array([[[1.1, 1.1, 1.1, 1.1]]], dtype=np.float32)
    
            rhs = np.array([[[0.5, 0.5, 0.5, 0.5]]], dtype=np.float32)
    
            dt = 0.1
    
            c1 = 0.75
    
            c2 = 0.25
    
            active_indices = np.arange(num_elems, dtype=np.int32)
            params = EquationParams32()
    
            
    
            q_wp = wp.array(q, dtype=wp.vec4, device=device)
    
            q_1_wp = wp.array(q_1, dtype=wp.vec4, device=device)
    
            rhs_wp = wp.array(rhs, dtype=wp.vec4, device=device)
    
            q_out_wp = wp.zeros_like(q_wp)
    
            active_wp = wp.array(active_indices, dtype=wp.int32, device=device)
    
            
    
            wp.launch(
    
                kernel=common_kernels.rk_stage_2,
    
                dim=(num_elems, Np),
    
                inputs=[q_wp, q_1_wp, rhs_wp, dt, q_out_wp, c1, c2, active_wp, params],
    
                device=device
    
            )
    
            
    
            # 0.75 * Q_n + 0.25 * (Q_1 + dt * RHS)
    
            expected = 0.75 * q + 0.25 * (q_1 + dt * rhs)
    
            np.testing.assert_allclose(q_out_wp.numpy(), expected, atol=1e-6)
    
    
    
        def test_apply_filter_matrix(self, device):
    
            wp.init()
    
            num_elems = 1
    
            Np = 2
    
            q = np.array([[[1.0, 1.0, 1.0, 1.0], [2.0, 2.0, 2.0, 2.0]]], dtype=np.float32)
    
            # Identity filter
    
            F = np.eye(Np, dtype=np.float32)
    
            
    
            q_wp = wp.array(q, dtype=wp.vec4, device=device)
    
            F_wp = wp.array(F, dtype=wp.float32, device=device)
    
            q_out_wp = wp.zeros_like(q_wp)
    
            
    
            wp.launch(
    
                kernel=common_kernels.apply_filter_matrix,
    
                dim=(num_elems, Np),
    
                inputs=[q_wp, F_wp, q_out_wp],
    
                device=device
    
            )
    
            
    
            np.testing.assert_allclose(q_out_wp.numpy(), q, atol=1e-6)
    
    