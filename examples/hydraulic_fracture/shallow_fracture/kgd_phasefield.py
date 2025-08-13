from dolfin import Mesh
from mpi4py import MPI
from dolfin import set_log_level, LogLevel
import time
import sys

from variational_forms.phase_field import phase_field_problem, solve_step_staggered
from boundary_conditions import setup_bc, create_markers
from fields.history import HistoryField
from fields.phase import PhaseField
from fields.displacement import DisplacementField
from fields.stress import StressField
from output_utils import write_output, create_xml_output
from utils import compute_opening_overtime, setup_logging


logger = setup_logging("Shallow KGD Fracture")
set_log_level(LogLevel.ERROR)

if len(sys.argv) >= 2:
    CASE_DIR = sys.argv[1]
    MESH_NAME = sys.argv[2]
else:
    CASE_DIR = "output"
    MESH_NAME = "mesh"

logger.info("Setting up the simulation problem...")

RELACION_HL = 3 # el parametro de transicion l_c cubre minimamente 3 elementos de la malla
## MESHING ##
comm = MPI.COMM_WORLD
mesh = Mesh(comm, f"{CASE_DIR}/{MESH_NAME}.xml")

# Material models
E = 2e8
nu = 0.3
Gc = 2.7
Q0 = 1e-3

# Computar parametros de Lame
mu = E / (2 * (1 + nu))
lmbda = E*nu / ((1 + nu)*(1 - 2*nu))

# Fields definitions
history = HistoryField(mesh, lmbda, mu)
displacement = DisplacementField(mesh)
phase = PhaseField(mesh)
stress = StressField(mesh, lmbda, mu)

# Phasefield params
h_elem = mesh.hmin()
l_c = h_elem * RELACION_HL
l_init = l_c * 2.5

# Solver params
p_init = 1000 # Presion inicial

# Boundary Conditions
# en upper_face_free estoy definiendo la superficie libre para el caso shallow
bcs_u, bcs_phi = setup_bc(phase, displacement, l_init, h_elem, upper_face_free=True, symmetric=True)

markers = create_markers(mesh)

# Variational Forms
E_du, E_phi, pressure = phase_field_problem(
    phase, displacement, history,
    lmbda, mu, Gc,
    p_init, l_c
)

# Setup solvers
displacement.setup_solver(E_du, bcs_u)
phase.setup_solver(E_phi, bcs_phi)
phase.solve()

logger.info("Problem setup complete.")

logger.info("--- Starting Simulation ---")

# Solve to phase

saved_vtus = 0
progress = 0
size = MPI.COMM_WORLD.Get_size()

# Parametros de simulacion
t_max = 0.6
dt = 2.5e-4
fname = open(f"{CASE_DIR}/output.csv", 'w')
fname.write("time,pressure,volume,wplus,wminus\n")
store_freq = 50 # escribir vtk cada n-step
output_freq = 50 # Flush output cada n-step

out_xml = create_xml_output(CASE_DIR)

t = 0
step = 0
while t <= t_max:
    start_time_step = time.time()
    step += 1
    t += dt
    dV = dt*Q0

    try:
        pn, vol_frac = solve_step_staggered(displacement, phase, history, pressure, dV)

    except RuntimeError:
        pn, vol_frac = 1, 1
        write_output(out_xml, displacement.get(), phase.get(), stress.get(), t)
        logger.error("La presión no converge")
        break

    except Exception as e:
        logger.error(f"Simulation stopped due to error: {e}", exc_info=True)
        break

    finally:

        if MPI.COMM_WORLD.size == 1:
            w_plus, w_minus = compute_opening_overtime(displacement.get(), phase.get(), h_elem/20)
            
        else:
            w_plus, w_minus = 0, 0

        stress.update(phase.get(), displacement.get())
        
        fname.write(f"{t},{pn},{vol_frac},{w_plus},{w_minus}\n")

        if step % output_freq == 0:
            write_output(out_xml, displacement.get(), phase.get(), stress.get(), t)
            # export_phi_to_csv(pnew, mesh, CASE_DIR)
            saved_vtus += 1
            print("Saved VTU num:", saved_vtus)
        
        if step % store_freq == 0:
            fname.flush()
            # export_phi_to_csv(pnew, mesh, CASE_DIR)
        elapsed_time_step = time.time() - start_time_step
        if MPI.COMM_WORLD.rank == 0:
            logger.info(f"Time: {t:.2e}s | Step: {step}")

if fname:
    fname.close()
    logger.info("Closed output CSV file.") 