import pytest
import numpy as np
import os
import sys

# Ensure we can import the study script
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from study_convergence import run_study

def test_vortex_convergence_regression():
    """
    Runs the full isentropic vortex convergence study (P=1, 2, 3 on Depths 5, 6, 7)
    and asserts that the L2 errors and convergence rates match the expected baselines
    within a 5% tolerance.
    """
    
    # Expected results from baseline run (CUDA, Double Precision)
    # Format: P -> {'expected_slope': float, 'finest_error': float}
    # finest_error is the error at the highest depth (Depth 7)
    
    expected_metrics = {
        1: {
            'expected_slope': 1.761,
            'finest_error': 1.690555e-04
        },
        2: {
            'expected_slope': 2.549,
            'finest_error': 7.594914e-06
        },
        3: {
            'expected_slope': 3.805,
            'finest_error': 8.921087e-08
        }
    }
    
    # 5% Tolerance
    tolerance = 0.05
    
    print("\nStarting Regression Test: Vortex Convergence Study...")
    results = run_study()
    
    for P, metrics in expected_metrics.items():
        assert P in results, f"Missing results for Polynomial Order P={P}"
        
        data = results[P]
        errors = data['errors']
        h_values = data['h']
        
        # 1. Calculate actual slope (Global Order of Convergence)
        log_h = np.log(h_values)
        log_err = np.log(errors)
        coeffs = np.polyfit(log_h, log_err, 1)
        slope = coeffs[0]
        
        # 2. Get finest mesh error (last element)
        finest_error = errors[-1]
        
        expected_slope = metrics['expected_slope']
        expected_finest_error = metrics['finest_error']
        
        print(f"\n--- Checking P={P} ---")
        print(f"Slope: Actual={slope:.4f}, Expected={expected_slope:.4f}")
        print(f"Error (Finest): Actual={finest_error:.4e}, Expected={expected_finest_error:.4e}")
        
        # Check Slope: Should not be significantly worse (lower) than expected
        # We allow it to be *better* (higher), but flag if it drops by > 5%
        assert slope >= expected_slope * (1 - tolerance), \
            f"P={P} Convergence rate degraded! Actual: {slope:.4f}, Expected: {expected_slope:.4f}"
            
        # Check Error: Should not be significantly higher than expected
        # We allow it to be *smaller* (better), but flag if it increases by > 5%
        assert finest_error <= expected_finest_error * (1 + tolerance), \
            f"P={P} Error increased! Actual: {finest_error:.4e}, Expected: {expected_finest_error:.4e}"

    print("\nRegression Test Passed: Convergence rates and errors are within tolerance.")
