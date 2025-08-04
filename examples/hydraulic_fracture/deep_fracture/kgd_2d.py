import logging
from dolfin import TensorFunctionSpace, Function, Mesh
from mpi4py import MPI

from variational_forms.phase_field import phase_field_problem
from boundary_conditions import setup_bc, create_markers
from fields.history import HistoryField
from fields.phase import PhaseField
from fields.displacement import DisplacementField

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.StreamHandler()
        ]
    )

setup_logging()

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

# Setup solvers within the field objects
displacement.setup_solver(E_du, bcs_u)
phase.setup_solver(E_phi, bcs_phi)

logger.info("Problem setup complete.") 