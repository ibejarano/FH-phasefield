from dolfin import TensorFunctionSpace, Function, TrialFunction, TestFunction, project, FunctionSpace
from dolfin import errornorm
import json
import time
import datetime
from sys import argv
from sys import stdout
import os
import functools

from mesh_setup import setup_gmsh, setup_rect_mesh, set_function_spaces
from material_model import epsilon, select_sigma, select_psi, compute_fracture_volume, get_E_expression
from variational_forms import define_variational_forms
from solvers import setup_solvers
from output_utils import create_output_files, write_output, store_time_series
from boundary_conditions import setup_boundary_conditions
from history_field import HistoryField
print = functools.partial(print, flush=True)

class Simulation:
    def __init__(self, data):
        self.data = data
        self.caseDir = data.get("caseDir", "results")
        self.p_init = data.get("p_init", 100) # Get initial pressure from data or default
        self.setup()

    def setup(self):
        """Initializes mesh, function spaces, variables, solvers, and output."""
        print("Setting up simulation...")
        ## MESHING ##

        if self.data["mesh_data"]["type"] == "rectangle":
            self.mesh, self.boundary_markers = setup_rect_mesh(self.data)
        elif self.data["mesh_data"]["type"] == "gmsh":
            self.mesh, self.boundary_markers = setup_gmsh(self.caseDir, self.data)
        else:
            RuntimeError("config mesh data not recognized")

        self.V, self.W = set_function_spaces(self.mesh)
        self.WW = FunctionSpace(self.mesh, "DG", 0)
        self.Vsig = TensorFunctionSpace(self.mesh, "DG", 0)
        # Trial and Test functions (used in variational forms)
        self.p, self.q = TrialFunction(self.V), TestFunction(self.V)
        self.u, self.v = TrialFunction(self.W), TestFunction(self.W)

        # Material models
        self.E_expr = get_E_expression(self.data)
        nu = self.data["nu"]
        self.psi = select_psi("linear")
        self.sigma = select_sigma("linear") # Or select based on data if needed
        
        self.history = HistoryField(self.WW, self.psi, self.E_expr, nu, self.data)
        print(f"Using energy model: Linear")

        # Boundary Conditions
        self.bc_u, self.bc_phi = setup_boundary_conditions(self.V, self.W, self.boundary_markers, self.data)

        # Functions
        self.unew, self.uold, self.ut = Function(self.W), Function(self.W), Function(self.W, name="displacement")
        self.pnew, self.pold, self.phit = Function(self.V), Function(self.V), Function(self.V, name="phi")
        self.sigt = Function(self.Vsig, name="stress")

        # Variational Forms and Solvers
        E_du, E_phi, self.pressure = define_variational_forms(
            epsilon, self.sigma, self.history.get(), self.pold, self.u,
            self.v, self.p, self.q, self.data, self.boundary_markers, self.E_expr, nu
        )

        self.solver_disp, self.solver_phi = setup_solvers(E_du, E_phi, self.unew, self.pnew, self.bc_u, self.bc_phi)
        
        self.solver_disp.solve()
        self.solver_phi.solve()

        # Output setup
        self.out_xml, self.u_ts, self.phi_ts = create_output_files(self.caseDir)
        self.fname = open(f"./{self.caseDir}/output.csv", 'w')
        self.fname.write("time,pressure,volume\n") # Add header to CSV

        # Simulation control variables
        self.t = 0.0
        self.step = 0
        self.pn = self.p_init
        self.pressure.assign(self.pn)

        # Initial solve/state setup
        print("Performing initial solve...")
        self.solver_phi.solve()
        self.pold.assign(self.pnew)

        # Save parameters used
        outfile = open(f"./{self.caseDir}/parameters_used.json", 'w')
        outfile.write(json.dumps(self.data, indent=4)) # Use indent for readability
        outfile.close()
        print("Setup complete.")

        self.E_func = project(self.E_expr, self.V)
        from dolfin import XDMFFile
        with XDMFFile(f"./{self.caseDir}/E_regions.xdmf") as xdmf:
            xdmf.write(self.E_func, 0.0)
        print(f"Saved E regions to {self.caseDir}/E_regions.xdmf")

    def adjust_pressure(self, V0):
        """Iteratively adjusts pressure to match target volume inflow."""
        errV = 1.0 # Initialize with a value > tolerance
        errV1 = None
        errV2 = None
        Q0 = self.data["Qo"]
        DT = self.data["dt"]
        Vtarget = V0 + DT * Q0
        pn = float(self.pressure) # Get current pressure value
        pn1 = pn # Initialize for secant method
        ite = 0
        vol_tol = self.data["tolerances"]["volume"]

        # Use Vtarget for relative tolerance check, avoid division by zero if Vtarget is 0
        check_val = Vtarget if abs(Vtarget) > 1e-15 else 1.0

        while abs(errV) / abs(check_val) > vol_tol:
            ite += 1
            # Secant method step
            if errV1 is not None and errV2 is not None and abs(errV1 - errV2) > 1e-12:
                pn = pn1 - errV1 * (pn1 - pn2) / (errV1 - errV2)
            else: # Fallback: small perturbation if secant fails or first steps
                pn *= 1.01

            self.pressure.assign(pn)
            self.solver_disp.solve()

            VK = compute_fracture_volume(self.pold, self.unew) # Use pold here as BCs depend on it
            errV = Vtarget - VK

            # Update state for next iteration/step
            self.uold.assign(self.unew)
            self.history.update(self.unew)

            # Store previous values for secant method
            pn2 = pn1
            errV2 = errV1
            pn1 = pn
            errV1 = errV

            if ite > 20: # Max iterations guard
                print(f"*** Warning: Adjust pressure reached max iterations ({ite}) with error {abs(errV)/abs(check_val):.2e} > {vol_tol:.2e}")
                return -1, pn # Indicate failure

        # print(f"Adjust pressure converged in {ite} iterations. P = {pn:.4e}, V_err = {errV:.4e}")
        return ite, pn

    def solve_phase_field(self):
        """Solves the phase-field equation and returns the error."""
        self.solver_phi.solve()
        err_phi = errornorm(self.pnew, self.pold, norm_type='l2', mesh=self.mesh)
        self.pold.assign(self.pnew)
        return err_phi

    def solve_step(self, V0):
        """Solves one coupled time step using staggered approach."""
        err_phi = 1.0 # Initialize error > tolerance
        phi_tol = self.data["tolerances"]["phi"]
        outer_ite = 0

        while err_phi > phi_tol:
            outer_ite += 1
            # Adjust pressure and solve displacement
            ite_p, pn = self.adjust_pressure(V0)

            if ite_p < 0:
                print("*** Error: Pressure adjustment failed to converge. Stopping simulation. ***")
                raise RuntimeError("Pressure adjustment failed.") # Raise exception instead of exit()

            # Solve phase field
            err_phi = self.solve_phase_field()

            # Max iterations guard for the outer loop
            if outer_ite > 15:
                 print(f"*** Warning: Outer loop reached max iterations ({outer_ite}) with phi error {err_phi:.2e} > {phi_tol:.2e}")
                 # Decide whether to continue or stop
                 # break # Option: continue with the current state
                 raise RuntimeError("Outer staggered loop failed to converge.") # Option: stop

        # Update time-dependent functions after step converges
        self.ut.assign(self.unew)
        self.phit.assign(self.pnew)
        return pn # Return the converged pressure

    def run(self):
        """Runs the main simulation loop."""
        print("\n--- Starting Simulation ---")
        total_steps = int(self.data["t_max"] / self.data["dt"])
        start_time = time.time()
        print(f"Total steps: {total_steps} | Total time: {self.data['t_max']:.4f}")
        saved_vtus = 0
        total_vtus = int(total_steps / self.data.get("output_frequency", 10))
        try:
            while self.t <= self.data["t_max"]:
                self.step += 1
                self.t += self.data["dt"]

                # Calculate volume at the beginning of the step
                V0 = compute_fracture_volume(self.phit, self.ut)
                
                # Solve the coupled step
                self.pn = self.solve_step(V0)

                # Calculate volume at the end of the step
                vol_frac = compute_fracture_volume(self.phit, self.ut)

                # Compute and project stress before writing output
                stress_expr = (1-self.phit)**2 * self.sigma(self.ut, self.E_expr, self.data.get("nu", 0.3))
                self.sigt.assign(project(stress_expr, self.Vsig))

                # Write data to CSV
                self.fname.write(f"{self.t},{self.pn},{vol_frac}\n")


                elapsed = time.time() - start_time
                avg_time = elapsed / self.step
                remaining_steps = total_steps - self.step
                remaining_time = remaining_steps * avg_time
                end_time = datetime.datetime.now() + datetime.timedelta(seconds=remaining_time)
                

                # Save output files periodically
                if self.step % self.data.get("output_frequency", 10) == 0: # Use frequency from config or default
                    write_output(self.out_xml, self.ut, self.phit, self.sigt, self.t)
                    # store_time_series(self.u_ts, self.phi_ts, self.ut, self.phit, self.t)
                    self.fname.flush() # Ensure CSV data is written to disk
                    saved_vtus +=1
                msg = f"Estimated end: {end_time.strftime('%Y-%m-%d %H:%M:%S')} | Vtks saved: {saved_vtus}/{total_vtus} | Step: {self.step}/{total_steps}"
                print('\r' + msg + ' '*20, end='')
                stdout.flush()
                # Optional: Add convergence checks or other break conditions here

        except Exception as e:
            print(f"\n--- Simulation stopped due to error: {e} ---")
            import traceback
            traceback.print_exc()
        finally:
            # Cleanup: Close files
            print("\n--- Simulation Finished ---")
            if hasattr(self, 'fname') and self.fname:
                self.fname.close()
                print("Closed output CSV file.")