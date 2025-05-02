from dolfin import *
import json
from sys import argv
import os
from utils import read_data
from mesh_setup import setup_gmsh, setup_rect_mesh, set_function_spaces
from material_model import epsilon, select_sigma, select_psi, H, compute_fracture_volume
from variational_forms import define_variational_forms
from solvers import setup_solvers
from output_utils import create_output_files, write_output, store_time_series
from boundary_conditions import setup_boundary_conditions

class Simulation:
    def __init__(self, data):
        self.data = data
        self.caseDir = os.path.join("./results/", argv[1])
        self.p_init = data.get("p_init", 100) # Get initial pressure from data or default
        self.setup()

    def setup(self):
        """Initializes mesh, function spaces, variables, solvers, and output."""
        print("Setting up simulation...")
        ## MESHING ##
        self.mesh, self.boundary_markers = setup_gmsh(self.caseDir, self.data)
        self.V, self.W, self.WW = set_function_spaces(self.mesh)

        # Trial and Test functions (used in variational forms)
        self.p, self.q = TrialFunction(self.V), TestFunction(self.V)
        self.u, self.v = TrialFunction(self.W), TestFunction(self.W)

        # Material models
        psi_model = self.data.get("psi_model", "linear")
        self.psi = select_psi(psi_model)
        self.sigma = select_sigma("linear") # Or select based on data if needed
        print(f"Using energy model: {psi_model}")

        # Boundary Conditions
        self.bc_u, self.bc_phi = setup_boundary_conditions(self.V, self.W, self.boundary_markers, self.data)

        # Functions
        self.unew, self.uold, self.ut = Function(self.W), Function(self.W), Function(self.W, name="displacement")
        self.pnew, self.pold, self.Hold, self.phit = Function(self.V), Function(self.V), Function(self.V), Function(self.V, name="phi")

        # Variational Forms and Solvers
        E_du, E_phi, self.pressure = define_variational_forms(
            self.W, self.V, epsilon, self.sigma, H, self.psi, self.pold,
            self.u, self.v, self.p, self.q, self.unew, self.Hold, self.data, self.boundary_markers
        )
        self.solver_disp, self.solver_phi = setup_solvers(E_du, E_phi, self.unew, self.pnew, self.bc_u, self.bc_phi)
        
        self.solver_disp.solve()
        self.solver_phi.solve()

         # Output setup
        self.out_xml, self.u_ts, self.phi_ts = create_output_files(self.mesh, self.caseDir)
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
            self.Hold.assign(project(self.psi(self.unew, self.data), self.WW))

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

                # Write data to CSV
                self.fname.write(f"{self.t},{self.pn},{vol_frac}\n")

                print(f"Step: {self.step} | Time: {self.t:.4f} | Pressure: {self.pn:.4e} | Volume: {vol_frac:.4e}")

                # Save output files periodically
                if self.step % self.data.get("output_frequency", 10) == 0: # Use frequency from config or default
                    print(f"  Saving output at step {self.step}...")
                    write_output(self.out_xml, self.ut, self.phit, self.t, self.step)
                    store_time_series(self.u_ts, self.phi_ts, self.ut, self.phit, self.t)
                    self.fname.flush() # Ensure CSV data is written to disk

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