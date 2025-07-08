from dolfin import lhs, rhs, LinearVariationalProblem, LinearVariationalSolver
import logging
from scipy.optimize import root_scalar, minimize_scalar, brentq
import numpy as np
from .material_model import compute_fracture_volume

logger = logging.getLogger(__name__)

def setup_solvers(E_du, E_phi, unew, pnew, bc_u, bc_phi):
    p_disp = LinearVariationalProblem(lhs(E_du), rhs(E_du), unew, bc_u)
    p_phi = LinearVariationalProblem(lhs(E_phi), rhs(E_phi), pnew, bc_phi)

    solver_disp = LinearVariationalSolver(p_disp)
    solver_phi = LinearVariationalSolver(p_phi)

    solver_phi.parameters["linear_solver"] = "gmres"
    solver_phi.parameters["preconditioner"] = "ilu"

    return solver_disp, solver_phi


def pressure_solver(Vtarget, phase, displacement, history, pressure, vol_tol, method='brentq'):
    """
    Solve for pressure using scipy optimization methods.
    
    Args:
        method: 'brentq', 'root_scalar', 'minimize_scalar', or 'secant' (original method)
    """
    
    def objective_function(pn):
        """Objective function: difference between target and computed volume"""
        pressure.assign(pn)
        displacement.solve()
        unew = displacement.get()
        history.update(unew)
        
        VK = compute_fracture_volume(phase.get_old(), unew)
        return Vtarget - VK
    
    def volume_function(pn):
        """Volume function for root finding"""
        pressure.assign(pn)
        displacement.solve()
        unew = displacement.get()
        history.update(unew)
        
        return compute_fracture_volume(phase.get_old(), unew)
    
    # Initial pressure estimate
    pn_initial = float(pressure)
    
    try:
        if method == 'brentq':
            # Brent's method - very robust for root finding
            result = brentq(
                lambda p: volume_function(p) - Vtarget,
                a=0.0,  # Lower bound
                b=pn_initial * 10000,  # Upper bound (adjust as needed)
                xtol=vol_tol,
                maxiter=50
            )
            pn = result
            iterations = 1  # brentq doesn't return iteration count easily
            
        elif method == 'root_scalar':
            # General root finding with different methods
            result = root_scalar(
                lambda p: volume_function(p) - Vtarget,
                x0=pn_initial,
                x1=pn_initial * 1.1,  # Second guess for secant method
                method='secant',
                xtol=vol_tol,
                maxiter=50
            )
            pn = result.root
            iterations = result.iterations
            
        elif method == 'minimize_scalar':
            # Minimize the absolute difference
            result = minimize_scalar(
                lambda p: abs(objective_function(p)),
                bounds=(0.0, pn_initial * 10),
                method='bounded',
                options={'xatol': vol_tol, 'maxiter': 50}
            )
            pn = result.x
            iterations = result.nit
            
        elif method == 'secant':
            # Original secant method for comparison
            return _original_secant_method(Vtarget, phase, displacement, history, pressure, vol_tol)
            
        else:
            raise ValueError(f"Unknown method: {method}. Use 'brentq', 'root_scalar', 'minimize_scalar', or 'secant'")
        
        # Final check
        final_error = abs(objective_function(pn)) / (abs(Vtarget) if abs(Vtarget) > 1e-15 else 1.0)

        if final_error > vol_tol:
            logger.warning(f"Method {method} converged with error {final_error:.2e} > {vol_tol:.2e}")
            return -1, pn
        
        logger.debug(f"Pressure solver ({method}) converged in {iterations} iterations")
        return iterations, pn
        
    except Exception as e:
        logger.error(f"Pressure solver ({method}) failed: {e}")
        return -1, pn_initial

def _original_secant_method(Vtarget, phase, displacement, history, pressure, vol_tol):
    """Original secant method for comparison"""
    check_val = Vtarget if abs(Vtarget) > 1e-15 else 1.0
    
    pn = float(pressure)
    prevs = []
    ite = 0
    max_ite = 20
    pold = phase.get_old()
    
    while True:
        ite += 1
        pressure.assign(pn)
        displacement.solve()
        unew = displacement.get()
        history.update(unew)
        
        VK = compute_fracture_volume(pold, unew)
        errV = Vtarget - VK
        
        prevs.append((pn, VK))
        if len(prevs) > 2: 
            prevs.pop(0)
        
        if abs(errV) / abs(check_val) < vol_tol: 
            break

        if len(prevs) == 2:
            print("LLegando", prevs[1][1] - prevs[0][1])

        if len(prevs) == 2 and abs(prevs[1][1] - prevs[0][1]) > 1e-14:
            pn = prevs[1][0] + (Vtarget - prevs[1][1]) * (prevs[1][0] - prevs[0][0]) / (prevs[1][1] - prevs[0][1])
        else:
            pn *= 1.01
        
        pn = max(pn, 0.0)
        if ite > max_ite:
            logger.warning(f"Original secant method reached max iterations ({ite}) with error {abs(errV)/abs(check_val):.2e} > {vol_tol:.2e}")
            return -1, pn
    
    return ite, pn
