#!/usr/bin/env python3
"""
Test script to compare different pressure solver methods.
This demonstrates how to use the new SciPy-based pressure solvers.
"""

import time
import numpy as np
from src.solvers import pressure_solver

def test_pressure_solver_methods():
    """Test and compare different pressure solver methods."""
    
    # Test parameters
    methods = ['brentq', 'root_scalar', 'minimize_scalar', 'secant']
    Vtarget = 1e-6  # Target volume
    vol_tol = 1e-8  # Volume tolerance
    
    print("Testing Pressure Solver Methods")
    print("=" * 50)
    
    results = {}
    
    for method in methods:
        print(f"\nTesting method: {method}")
        
        # Mock objects for testing (you would use real FEniCS objects in practice)
        class MockPressure:
            def __init__(self):
                self.value = 100.0
            def assign(self, p):
                self.value = p
        class MockDisplacement:
            def solve(self):
                pass
            def get(self):
                return None
        class MockPhase:
            def get_old(self):
                return None
        class MockHistory:
            def update(self, u):
                pass
        
        # Create mock objects
        pressure = MockPressure()
        displacement = MockDisplacement()
        phase = MockPhase()
        history = MockHistory()
        
        # Time the solver
        start_time = time.time()
        
        try:
            iterations, final_pressure = pressure_solver(
                Vtarget=Vtarget,
                phase=phase,
                displacement=displacement,
                history=history,
                pressure=pressure,
                vol_tol=vol_tol,
                method=method
            )
            
            elapsed_time = time.time() - start_time
            
            results[method] = {
                'iterations': iterations,
                'final_pressure': final_pressure,
                'elapsed_time': elapsed_time,
                'success': iterations > 0
            }
            
            print(f"  Success: {iterations > 0}")
            print(f"  Iterations: {iterations}")
            print(f"  Final pressure: {final_pressure:.6f}")
            print(f"  Time: {elapsed_time:.4f} seconds")
            
        except Exception as e:
            results[method] = {
                'iterations': -1,
                'final_pressure': pressure.value,
                'elapsed_time': time.time() - start_time,
                'success': False,
                'error': str(e)
            }
            print(f"  Failed: {e}")
    
    # Summary
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    
    successful_methods = [m for m, r in results.items() if r['success']]
    
    if successful_methods:
        print(f"Successful methods: {', '.join(successful_methods)}")
        
        # Find fastest method
        fastest = min(successful_methods, key=lambda m: results[m]['elapsed_time'])
        print(f"Fastest method: {fastest} ({results[fastest]['elapsed_time']:.4f}s)")
        
        # Find most efficient method (fewest iterations)
        most_efficient = min(successful_methods, key=lambda m: results[m]['iterations'])
        print(f"Most efficient method: {most_efficient} ({results[most_efficient]['iterations']} iterations)")
    
    else:
        print("No methods succeeded!")
    
    return results

if __name__ == "__main__":
    test_pressure_solver_methods() 