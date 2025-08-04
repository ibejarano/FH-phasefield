import logging
from dolfin import TensorFunctionSpace, Function, Mesh
from mpi4py import MPI
import time

from variational_forms.phase_field import phase_field_problem, sigma, compute_fracture_volume
from boundary_conditions import setup_bc, create_markers
from fields.history import HistoryField
from fields.phase import PhaseField
from fields.displacement import DisplacementField
from fields.stress import StressField
from output_utils import write_output
from utils import compute_opening_overtime, export_phi_to_csv

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.StreamHandler()
        ]
    )

setup_logging()

CASE_DIR = "output"

logger = logging.getLogger(__name__)
logger.info("Setting up the simulation problem...")

RELACION_HL = 3 # el parametro de transicion l_c cubre minimamente 3 elementos de la malla
## MESHING ##
comm = MPI.COMM_WORLD
mesh = Mesh(comm, "deep_fh.xml")


Vsig = TensorFunctionSpace(mesh, "DG", 0)

# Material models
E = 2e8
nu = 0.3
Gc = 2.7

# Computar parametros de Lame
mu = E / (2 * (1 + nu))
lmbda = E*nu / ((1 + nu)*(1 - 2*nu))

# Fields definitions
history = HistoryField(mesh, lmbda, mu)
displacement = DisplacementField(mesh)
phase = PhaseField(mesh)
stress = StressField(mesh)

# Phasefield params
l_init = 0.1
h_elem = 1e-3
l_c = h_elem * RELACION_HL

# Solver params
p_init = 100 # Presion inicial

# Boundary Conditions
bcs_u, bcs_phi = setup_bc(phase, displacement, l_init, h_elem)

markers = create_markers(mesh)

# Functions
sigt = Function(Vsig, name="stress")

# Variational Forms
E_du, E_phi, p_vec = phase_field_problem(
    phase, displacement, history,
    lmbda, mu, Gc,
    p_init, l_c
)

# Setup solvers
displacement.setup_solver(E_du, bcs_u)
phase.setup_solver(E_phi, bcs_phi)

logger.info("Problem setup complete.")


logger.info("--- Starting Simulation ---")
saved_vtus = 0
progress = 0
size = MPI.COMM_WORLD.Get_size()

# Parametros de simulacion
t_max = 0.5
dt = 2.5e-5
fname = open(f"{CASE_DIR}/output.csv", 'w')
fname.write("time,pressure,volume,wplus,wminus\n")
store_freq = 1
output_freq = 1

try:
    while t <= t_max and progress < 0.95:
        start_time_step = time.time()
        step += 1
        t += dt

        pn = self._solve_step() # Solver Staggered
        
        pnew = phase.get()
        unew = displacement.get()

        vol_frac = compute_fracture_volume(pnew, unew)
        w_plus, w_minus = compute_opening_overtime(unew, pnew, h_elem/20)

        stress.update(pnew, unew)

        
        fname.write(f"{t},{pn},{vol_frac},{w_plus},{w_minus}\n")

        if step % output_freq == 0:
            write_output(out_xml, unew, pnew, stress.get(), t)
            saved_vtus += 1
        
        if step % store_freq == 0:
            fname.flush()
            export_phi_to_csv(pnew, mesh, CASE_DIR)

        elapsed_time_step = time.time() - start_time_step
        if MPI.COMM_WORLD.rank == 0:
            logger.info(f"Time: {t:.2e}s | Step: {step}")

except Exception as e:
    logger.error(f"Simulation stopped due to error: {e}", exc_info=True)
finally:
    logger.info("Simulation Finished.")
    if fname:
        fname.close()
        logger.info("Closed output CSV file.") 