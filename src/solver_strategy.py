import logging
import json
import time
from mpi4py import MPI
from dolfin import project

from .problem import Problem
from .config import Config
from .output_utils import create_output_files, write_output
from .material_model import compute_fracture_volume
from .utils import fracture_length, compute_opening_overtime
from .solvers import pressure_solver

logger = logging.getLogger(__name__)

class StaggeredSolver:
    """
    Manages the simulation run, including the time-stepping loop,
    staggered solution scheme, pressure adjustments, and output.
    """
    def __init__(self, problem: Problem, config: Config):
        self.problem = problem
        self.config = config
        self.case_dir = config.case_dir

        # Unpack problem components for easier access
        self.displacement = problem.displacement
        self.phase = problem.phase
        self.history = problem.history
        self.pressure = problem.pressure
        self.sigt = problem.sigt
        self.Vsig = problem.Vsig

        # Simulation control variables
        self.t = 0.0
        self.step = 0
        p_init = config.get("p_init", 100)
        self.pn = p_init
        self.pn_old = p_init
        self.pressure.assign(p_init)

        # Adaptive time step parameters
        solver_params = config.get("solver_parameters", {})
        self.dt = solver_params.get("dt")
        self.dt_min = solver_params.get("dt_min", 1e-5)
        self.dt_max = solver_params.get("dt_max", 1e-2)
        self.dt_growth = solver_params.get("dt_growth", 1.2)
        self.dt_shrink = solver_params.get("dt_shrink", 0.5)
        self.t_max = solver_params.get("t_max", 1.0)

        # Tolerances
        self.vol_tol = solver_params.get("tolerances", {}).get("volume", 1e-5)
        self.phi_tol = solver_params.get("tolerances", {}).get("phi", 0.5e-3)

        self._setup_output()

    def _setup_output(self):
        """Initializes output files and saves initial state."""
        self.case_dir.mkdir(parents=True, exist_ok=True)
        self.out_xml, self.u_ts, self.phi_ts = create_output_files(self.case_dir)
        self.fname = open(self.case_dir / "output.csv", 'w')
        self.fname.write("time,pressure,volume,length,wplus,wminus\n")

        # Save parameters used
        with open(self.case_dir / "parameters_used.json", 'w') as f:
            json.dump(self.config.params, f, indent=4)
        
        logger.info("Performing initial solve...")
        self.phase.solve()
        logger.info("Setup complete.")

    def _adjust_pressure(self, V0):
        """Use SciPy-based pressure solver instead of manual secant method."""
        Q0 = self.config["material_parameters"]["Qo"]
        Vtarget = V0 + self.dt * Q0
        
        # Get the optimization method from config, default to 'brentq'
        method = self.config.get("pressure_solver_method", "brentq")
        
        # Call the new pressure solver
        ite_p, pn = pressure_solver(
            Vtarget=Vtarget,
            phase=self.phase,
            displacement=self.displacement,
            history=self.history,
            pressure=self.pressure,
            vol_tol=self.vol_tol,
            method=method
        )
        
        return ite_p, pn

    def _solve_step(self):
        # (Code from Simulation.solve_step)
        err_phi = 1.0
        outer_ite = 0
        V0 = compute_fracture_volume(self.phase.get_old(), self.displacement.get())

        while err_phi > self.phi_tol:
            outer_ite += 1
            ite_p, pn = self._adjust_pressure(V0)
            if ite_p < 0:
                raise RuntimeError("Pressure adjustment failed to converge.")
            err_phi = self.phase.solve()
            if outer_ite > 15:
                raise RuntimeError(f"Outer staggered loop failed to converge (err_phi={err_phi:.2e})")

        self.displacement.update()
        self.phase.update()
        return pn

    def run(self):
        # (Code from Simulation.run)
        logger.info("--- Starting Simulation ---")
        saved_vtus = 0
        final_length = self.config.get("l_max", 0.0)
        progress = 0
        size = MPI.COMM_WORLD.Get_size()

        try:
            while self.t <= self.t_max and progress < 0.95:
                start_time_step = time.time()
                self.step += 1
                self.t += self.dt

                self.pn = self._solve_step()
                
                pnew = self.phase.get()
                unew = self.displacement.get()

                vol_frac = compute_fracture_volume(pnew, unew)
                w_plus, w_minus = compute_opening_overtime(unew, pnew, self.config["meshing_parameters"]["h"]/20)

                stress_expr = (1-pnew)**2 * self.problem.sigma(unew, self.problem.E_expr, self.config.params["material_parameters"].get("nu", 0.3))
                self.sigt.assign(project(stress_expr, self.Vsig))

                if size == 1:
                    fracture_length_value = fracture_length(pnew, x2=final_length/2, x1=-final_length/2)
                else:
                    fracture_length_value = 1e-2 # Placeholder for MPI
                
                self.fname.write(f"{self.t},{self.pn},{vol_frac},{fracture_length_value},{w_plus},{w_minus}\\n")

                if self.step % self.config.get("output_frequency", 10) == 0:
                    write_output(self.out_xml, unew, pnew, self.sigt, self.t)
                    saved_vtus += 1
                
                if self.step % self.config.get("store_frequency", 1) == 0:
                    self.fname.flush()

                elapsed_time_step = time.time() - start_time_step
                if MPI.COMM_WORLD.rank == 0:
                    if final_length > 0:
                        progress = min(fracture_length_value / final_length, 1.0)
                        logger.info(
                            f"Progress: {progress*100:.1f}% | Time: {self.t:.2e}s | Step: {self.step} | dt: {self.dt:.2e} | "
                            f"Step time: {elapsed_time_step:.2f}s"
                        )
                    else:
                        logger.info(f"Time: {self.t:.2e}s | Step: {self.step} | Length: {fracture_length_value:.4f}")

        except Exception as e:
            logger.error(f"Simulation stopped due to error: {e}", exc_info=True)
        finally:
            logger.info("Simulation Finished.")
            if self.fname:
                self.fname.close()
                logger.info("Closed output CSV file.") 