from dolfin import Mesh, Point, RectangleMesh, MeshFunction, refine, cells
from mpi4py import MPI
from dolfin import set_log_level, LogLevel
import time
import math

from variational_forms.phase_field import phase_field_problem, solve_step_staggered
from boundary_conditions import setup_bc, create_markers
from fields.history import HistoryField
from fields.phase import PhaseField
from fields.displacement import DisplacementField
from fields.stress import StressField
from output_utils import write_output, create_xml_output
from utils import compute_opening_overtime, setup_logging


logger = setup_logging("Shallow Fracture")
set_log_level(LogLevel.ERROR)

CASE_DIR = "refined_output"

logger.info("Setting up the simulation problem...")

RELACION_HL = 4 # el parametro de transicion l_c cubre minimamente 3 elementos de la malla
## MESHING ##

nx, ny = 20, 20
Lx, Ly = 10, 5
Htop = 0.5
Hbottom = 2*Ly - Htop

comm = MPI.COMM_WORLD
mesh = RectangleMesh(comm, Point(0, -Hbottom), Point(Lx, Htop), nx, ny)

# --- 2) trayectoria como polilínea (lista de puntos)

path = [Point(0,0), Point(0.336,0.0313), Point(0.453,0.0625), Point(0.531,0.09375) , Point(0.6875,0.187) , Point(0.84375,0.281), Point(0.9375,0.375)]

# utilidades geométricas
def dist_point_to_segment(p, a, b):
    # distancia euclídea de p al segmento [a,b]
    ax, ay = a.x(), a.y()
    bx, by = b.x(), b.y()
    px, py = p.x(), p.y()
    abx, aby = bx-ax, by-ay
    apx, apy = px-ax, py-ay
    ab2 = abx*abx + aby*aby
    if ab2 == 0.0:
        # a==b
        dx, dy = px-ax, py-ay
        return math.hypot(dx, dy)
    t = (apx*abx + apy*aby)/ab2
    t = max(0.0, min(1.0, t))
    cx, cy = ax + t*abx, ay + t*aby
    return math.hypot(px-cx, py-cy)

def dist_point_to_polyline(p, poly):
    dmin = float("inf")
    for i in range(len(poly)-1):
        dmin = min(dmin, dist_point_to_segment(p, poly[i], poly[i+1]))
    return dmin

# --- 3) marcadores de refinamiento por proximidad a la trayectoria
def refine_near_path(mesh, path_points, radius):
    markers = MeshFunction("bool", mesh, mesh.topology().dim(), False)
    for cell in cells(mesh):
        if dist_point_to_polyline(cell.midpoint(), path_points) < radius:
            markers[cell] = True
    return refine(mesh, markers)  # devuelve nueva malla

# --- 4) una sola capa de refinamiento (tubo de radio r)
r = 0.8
mesh = refine_near_path(mesh, path, r)
for r in [0.60, 0.40, 0.25, 0.1, 0.1]:
    mesh = refine_near_path(mesh, path, r)

# mesh = Mesh(comm, f"{CASE_DIR}/curving.xml")

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
t_max = 3e-1
dt = 2.5e-4
fname = open(f"{CASE_DIR}/output.csv", 'w')
fname.write("time,pressure,volume,wplus,wminus\n")
store_freq = 1 # escribir vtk cada n-step
output_freq = 1 # Flush output cada n-step

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
        
        if step % store_freq == 0:
            fname.flush()
            # export_phi_to_csv(pnew, mesh, CASE_DIR)
        elapsed_time_step = time.time() - start_time_step
        if MPI.COMM_WORLD.rank == 0:
            logger.info(f"Time: {t:.2e}s | Step: {step}")

if fname:
    fname.close()
    logger.info("Closed output CSV file.") 