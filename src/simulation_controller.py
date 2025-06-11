import logging
from dolfin import TensorFunctionSpace, Function, project
import json
import time
from mpi4py import MPI

from mesh_setup import setup_gmsh, setup_rect_mesh
from material_model import epsilon, select_sigma, select_psi, compute_fracture_volume, get_E_expression
from variational_forms import define_variational_forms
from output_utils import create_output_files, write_output
from boundary_conditions import setup_shallow_bc, create_markers
from fields.history import HistoryField
from fields.phase import PhaseField
from fields.displacement import DisplacementField
from utils import fracture_length, compute_opening_overtime

logger = logging.getLogger(__name__)

class Simulation:
    def __init__(self, data):
        self.data = data
        self.caseDir = data.get("caseDir", "results")
        self.p_init = data.get("p_init", 100)
        self.pn_old = self.p_init
        self.setup()

    def setup(self):
        """Initializes mesh, function spaces, variables, solvers, and output."""
        logger.info("Setting up simulation...")
        ## MESHING ##
        if self.data["mesh_data"]["type"] == "rectangle":
            self.mesh = setup_rect_mesh(self.data)
        elif self.data["mesh_data"]["type"] == "gmsh":
            self.mesh = setup_gmsh(self.caseDir, self.data)
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
        self.bc_u, self.bc_phi = setup_shallow_bc(self.phase, self.displacement, self.data)
        markers = create_markers(self.mesh)
        # Functions
        self.sigt = Function(self.Vsig, name="stress")

        # Variational Forms and Solvers
        E_du, E_phi, self.pressure = define_variational_forms(
            epsilon, self.sigma, self.history.get(), self.phase, self.displacement,
            self.data, self.E_expr, nu, markers
        )

        self.displacement.setup_solver(E_du, self.bc_u)
        self.phase.setup_solver(E_phi, self.bc_phi)

        # Output setup
        self.out_xml, self.u_ts, self.phi_ts = create_output_files(self.caseDir)
        self.fname = open(f"./{self.caseDir}/output.csv", 'w')
        self.fname.write("time,pressure,volume,length,wplus,wminus\n")

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

        # Adaptive time step parameters
        self.dt = self.data["dt"]
        self.dt_min = self.data.get("dt_min", 1e-5)
        self.dt_max = self.data.get("dt_max", 1e-2)
        self.dt_growth = self.data.get("dt_growth", 1.2)
        self.dt_shrink = self.data.get("dt_shrink", 0.5)

    def adjust_pressure(self, V0):
        """Iteratively adjusts pressure to match target volume inflow."""
        Q0 = self.data["Qo"]
        Vtarget = V0 + self.dt * Q0
        vol_tol = self.data["tolerances"]["volume"]
        check_val = Vtarget if abs(Vtarget) > 1e-15 else 1.0

        pn = float(self.pressure)
        prevs = []  # Guarda (pn, VK) para método secante
        ite = 0
        max_ite = 20

        pold = self.phase.get_old()
        errV = 1.0

        while True:
            ite += 1

            self.pressure.assign(pn)
            self.displacement.solve()
            unew = self.displacement.get()
            self.history.update(unew)

            VK = compute_fracture_volume(pold, unew)
            errV = Vtarget - VK
            #logger.info(f"Iteration {ite}: Pressure={pn:.4f}, Volume={VK:.6e}, Target={Vtarget:.6e}, Error={(errV/Vtarget):.2e}")

            prevs.append((pn, VK))
            if len(prevs) > 2:
                prevs.pop(0)

            # Chequeo de convergencia
            if abs(errV) / abs(check_val) < vol_tol:
                break

            # Método secante si hay dos valores previos y diferencia de volumen suficiente
            if len(prevs) == 2 and abs(prevs[1][1] - prevs[0][1]) > 1e-12:
                # Interpolación lineal para encontrar presión que da Vtarget
                pn = prevs[1][0] + (Vtarget - prevs[1][1]) * (prevs[1][0] - prevs[0][0]) / (prevs[1][1] - prevs[0][1])
            else:
                pn *= 1.01  # Incremento simple si no hay suficiente info

            # Opcional: evita presión negativa
            pn = max(pn, 0.0)

            if ite > max_ite:
                logger.warning(f"Adjust pressure reached max iterations ({ite}) with error {abs(errV)/abs(check_val):.2e} > {vol_tol:.2e}")
                return -1, pn

        return ite, pn

    def solve_step(self):
        """Solves one coupled time step using staggered approach."""
        err_phi = 1.0
        phi_tol = self.data["tolerances"]["phi"]
        outer_ite = 0
        V0 = compute_fracture_volume(self.phase.get_old(), self.displacement.get())

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
        start_time = time.time()
        saved_vtus = 0
        final_length = self.data.get("l_max", 0.0)
        progress = 0
        size = MPI.COMM_WORLD.Get_size()
        try:
            while self.t <= self.data["t_max"] and progress < 0.95:
                start_time_step = time.time()
                self.step += 1
                self.t += self.dt

                # Guardar phi anterior para dt adaptativo
                self.pn = self.solve_step()
                
                pnew = self.phase.get()
                unew = self.displacement.get()

                vol_frac = compute_fracture_volume(pnew, unew)
                w_plus, w_minus = compute_opening_overtime(unew, pnew, self.data["h"]/20)

                stress_expr = (1-pnew)**2 * self.sigma(unew, self.E_expr, self.data.get("nu", 0.3))
                self.sigt.assign(project(stress_expr, self.Vsig))
                if size == 1:
                    fracture_length_value = fracture_length(pnew, x2=final_length/2, x1=-final_length/2)
                else:
                    fracture_length_value = 1e-2
                self.fname.write(f"{self.t},{self.pn},{vol_frac},{fracture_length_value},{w_plus},{w_minus}\n")

                pn_new = self.pn
                # Cambio relativo de presión
                delta_p = abs(pn_new - self.pn_old) / max(abs(self.pn_old), 1e-8)

                # Adaptar dt según cambio de presión
                if delta_p > self.data.get("pressure_tol_adapt", 0.01) and self.dt > self.dt_min:
                    self.dt = max(self.dt * self.dt_shrink, self.dt_min)
                    logger.info(f"Decreasing dt to {self.dt:.2e} (pressure drop: {delta_p:.2e})")
                elif delta_p < self.data.get("pressure_tol_adapt", 0.005) and self.dt < self.dt_max:
                    self.dt = min(self.dt * self.dt_growth, self.dt_max)
                    logger.info(f"Increasing dt to {self.dt:.2e} (pressure stable: {delta_p:.2e})")

                self.pn_old = pn_new

                if self.step % self.data.get("output_frequency", 10) == 0:
                    write_output(self.out_xml, unew, pnew, self.sigt, self.t)
                    saved_vtus += 1

                if self.step % self.data.get("store_frequency", 5) == 0:
                    self.fname.flush()

                elapsed_time_step = time.time() - start_time_step
                if MPI.COMM_WORLD.rank == 0:
                    if final_length > 0:
                        progress = min(fracture_length_value / final_length, 1.0)
                        if progress > 0:
                            msg = (
                                f"Fracture progress: {progress*100:.1f}% | "
                                f"Time: {self.t:.2e} s | "
                                f"Step/VTUS: {self.step}/{saved_vtus} | "
                                f" | completed in: {elapsed_time_step:.2f} s | "
                                f" fr length {fracture_length_value:.4f} | "
                            )
                            logger.info(msg)
                        else:
                            logger.info("Waiting for fracture to start growing to estimate end time.")
                    else:
                        logger.info(f"Fracture length: {fracture_length_value:.4f}")

        except Exception as e:
            logger.error(f"Simulation stopped due to error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            logger.info("Simulation Finished.")
            if hasattr(self, 'fname') and self.fname:
                self.fname.close()
                logger.info("Closed output CSV file.")