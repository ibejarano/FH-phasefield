import logging
from dolfin import TensorFunctionSpace, Function, project
import json
import time
import datetime
from sys import stdout
import functools

from mesh_setup import setup_gmsh, setup_rect_mesh
from material_model import epsilon, select_sigma, select_psi, compute_fracture_volume, get_E_expression
from variational_forms import define_variational_forms
from output_utils import create_output_files, write_output, store_time_series
from boundary_conditions import setup_boundary_conditions
from fields.history import HistoryField
from fields.phase import PhaseField
from fields.displacement import DisplacementField

logger = logging.getLogger(__name__)

class Simulation:
    def __init__(self, data):
        self.data = data
        self.caseDir = data.get("caseDir", "results")
        self.p_init = data.get("p_init", 100)
        self.setup()

    def setup(self):
        """Initializes mesh, function spaces, variables, solvers, and output."""
        logger.info("Setting up simulation...")
        ## MESHING ##
        if self.data["mesh_data"]["type"] == "rectangle":
            self.mesh, self.boundary_markers = setup_rect_mesh(self.data)
        elif self.data["mesh_data"]["type"] == "gmsh":
            self.mesh, self.boundary_markers = setup_gmsh(self.caseDir, self.data)
        else:
            logger.error("config mesh data not recognized")
            raise RuntimeError("config mesh data not recognized")

        self.Vsig = TensorFunctionSpace(self.mesh, "DG", 0)

        # Material models
        self.E_expr = get_E_expression(self.data)
        nu = self.data["nu"]
        self.psi = select_psi("linear")
        self.sigma = select_sigma("linear")
        logger.info("Using energy model: Linear")

        # Fields definitions
        self.history = HistoryField(self.mesh, self.psi, self.E_expr, nu, self.data)
        self.displacement = DisplacementField(self.mesh)
        self.phase = PhaseField(self.mesh)

        # Boundary Conditions
        self.bc_u, self.bc_phi = setup_boundary_conditions(self.phase, self.displacement, self.boundary_markers, self.data)

        # Functions
        self.sigt = Function(self.Vsig, name="stress")

        # Variational Forms and Solvers
        E_du, E_phi, self.pressure = define_variational_forms(
            epsilon, self.sigma, self.history.get(), self.phase, self.displacement,
            self.data, self.boundary_markers, self.E_expr, nu
        )

        self.displacement.setup_solver(E_du, self.bc_u)
        self.phase.setup_solver(E_phi, self.bc_phi)

        # Output setup
        self.out_xml, self.u_ts, self.phi_ts = create_output_files(self.caseDir)
        self.fname = open(f"./{self.caseDir}/output.csv", 'w')
        self.fname.write("time,pressure,volume\n")

        # Simulation control variables
        self.t = 0.0
        self.step = 0
        self.pn = self.p_init
        self.pressure.assign(self.pn)

        # Initial solve/state setup
        logger.info("Performing initial solve...")
        self.phase.solve()

        # Save parameters used
        outfile = open(f"./{self.caseDir}/parameters_used.json", 'w')
        outfile.write(json.dumps(self.data, indent=4))
        outfile.close()
        logger.info("Setup complete.")

    def adjust_pressure(self, V0):
        """Iteratively adjusts pressure to match target volume inflow."""
        errV = 1.0
        errV1 = None
        errV2 = None
        Q0 = self.data["Qo"]
        DT = self.data["dt"]
        Vtarget = V0 + DT * Q0
        pn = float(self.pressure)
        pn1 = pn
        ite = 0
        vol_tol = self.data["tolerances"]["volume"]
        check_val = Vtarget if abs(Vtarget) > 1e-15 else 1.0

        unew = self.displacement.get()
        pold = self.phase.get_old()

        while abs(errV) / abs(check_val) > vol_tol:
            ite += 1
            if errV1 is not None and errV2 is not None and abs(errV1 - errV2) > 1e-12:
                pn = pn1 - errV1 * (pn1 - pn2) / (errV1 - errV2)
            else:
                pn *= 1.01

            self.pressure.assign(pn)
            self.displacement.solve()
            self.history.update(unew)

            VK = compute_fracture_volume(pold, unew)
            errV = Vtarget - VK

            pn2 = pn1
            errV2 = errV1
            pn1 = pn
            errV1 = errV

            if ite > 20:
                logger.warning(f"Adjust pressure reached max iterations ({ite}) with error {abs(errV)/abs(check_val):.2e} > {vol_tol:.2e}")
                return -1, pn

        # logger.info(f"Adjust pressure converged in {ite} iterations. P = {pn:.4e}, V_err = {errV:.4e}")
        return ite, pn

    def solve_step(self, V0):
        """Solves one coupled time step using staggered approach."""
        err_phi = 1.0
        phi_tol = self.data["tolerances"]["phi"]
        outer_ite = 0

        while err_phi > phi_tol:
            outer_ite += 1
            ite_p, pn = self.adjust_pressure(V0)

            if ite_p < 0:
                logger.error("Pressure adjustment failed to converge. Stopping simulation.")
                raise RuntimeError("Pressure adjustment failed.")

            err_phi = self.phase.solve()

            if outer_ite > 15:
                logger.error(f"Outer loop reached max iterations ({outer_ite}) with phi error {err_phi:.2e} > {phi_tol:.2e}")
                raise RuntimeError("Outer staggered loop failed to converge.")

        self.displacement.update()
        self.phase.update()
        return pn

    def run(self):
        """Runs the main simulation loop."""
        logger.info("--- Starting Simulation ---")
        total_steps = int(self.data["t_max"] / self.data["dt"])
        start_time = time.time()
        logger.info(f"Total steps: {total_steps} | Total time: {self.data['t_max']:.4f}")
        saved_vtus = 0
        total_vtus = int(total_steps / self.data.get("output_frequency", 10))
        try:
            while self.t <= self.data["t_max"]:
                self.step += 1
                self.t += self.data["dt"]

                pnew = self.phase.get()
                unew = self.displacement.get()

                V0 = compute_fracture_volume(pnew, unew)
                self.pn = self.solve_step(V0)
                vol_frac = compute_fracture_volume(pnew, unew)

                stress_expr = (1-pnew)**2 * self.sigma(unew, self.E_expr, self.data.get("nu", 0.3))
                self.sigt.assign(project(stress_expr, self.Vsig))

                self.fname.write(f"{self.t},{self.pn},{vol_frac}\n")

                elapsed = time.time() - start_time
                avg_time = elapsed / self.step
                remaining_steps = total_steps - self.step
                remaining_time = remaining_steps * avg_time
                end_time = datetime.datetime.now() + datetime.timedelta(seconds=remaining_time)

                if self.step % self.data.get("output_frequency", 10) == 0:
                    write_output(self.out_xml, unew, pnew, self.sigt, self.t)
                    self.fname.flush()
                    saved_vtus += 1
                msg = f"Estimated end: {end_time.strftime('%Y-%m-%d %H:%M:%S')} | Vtks saved: {saved_vtus}/{total_vtus} | Step: {self.step}/{total_steps}"
                logger.info(msg)

        except Exception as e:
            logger.error(f"Simulation stopped due to error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            logger.info("Simulation Finished.")
            if hasattr(self, 'fname') and self.fname:
                self.fname.close()
                logger.info("Closed output CSV file.")