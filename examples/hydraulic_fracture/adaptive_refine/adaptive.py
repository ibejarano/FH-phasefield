from dolfin import *
from mpi4py import MPI
import time
import sys
from matplotlib import tri as mtri
import matplotlib.pyplot as plt
import numpy as np, math

from variational_forms.phase_field import phase_field_problem, solve_step_staggered
from boundary_conditions import setup_bc, create_markers
from fields.history import HistoryField
from fields.phase import PhaseField
from fields.displacement import DisplacementField
from fields.stress import StressField
from output_utils import write_output, create_xml_output
from utils import compute_opening_overtime, setup_logging


def extract_polyline_from_isocontour(mesh, f, level=0.8):
    # --- valores nodales (P1)
    V = FunctionSpace(mesh, "CG", 1)
    fn = interpolate(f, V)
    z = fn.compute_vertex_values(mesh)  # vector de escalares en vértices

    # --- triangulación (requiere malla triangular)
    xy = mesh.coordinates()
    tri_conn = mesh.cells()
    triang = mtri.Triangulation(xy[:,0], xy[:,1], tri_conn)

    # --- contorno
    fig = plt.figure()
    CS = plt.tricontour(triang, z, levels=[level])
    plt.close(fig)

    # --- obtener segmentos de forma compatible entre versiones
    segs = []
    try:
        # lista por nivel; tomamos el único nivel pedido
        segs = CS.allsegs[0]  # [array(Ni,2), array(Nj,2), ...]
    except Exception:
        pass
    if not segs:
        # fallback para versiones que sí traen .collections
        try:
            col = CS.collections[0]
            for path in col.get_paths():
                segs.append(path.vertices)  # array(N,2)
        except Exception:
            return None

    if not segs:
        return None

    # --- elegir la polilínea más larga
    lengths = [np.sum(np.linalg.norm(s[1:]-s[:-1], axis=1)) if len(s) > 1 else 0.0 for s in segs]
    s = segs[int(np.argmax(lengths))]
    return [Point(float(x), float(y)) for (x, y) in s]

def write_geo_with_path(domain, poly_xy, lc, lc_path, rin, rout, fname="remesh.geo"):
    # domain = (xmin, ymin, xmax, ymax) o podés escribir tu CAD aquí
    xmin, ymin, xmax, ymax = domain
    def pt(i, x, y): return f"Point({i}) = {{{x}, {y}, 0, {lc}}};\n"
    # cabecera y rectángulo
    geo = 'SetFactory("OpenCASCADE");\n'
    geo += f"Point(1) = {{{xmin}, {ymin}, 0, {lc}}};\n"
    geo += f"Point(2) = {{{xmax}, {ymin}, 0, {lc}}};\n"
    geo += f"Point(3) = {{{xmax}, {ymax}, 0, {lc}}};\n"
    geo += f"Point(4) = {{{xmin}, {ymax}, 0, {lc}}};\n"
    geo += "Line(11) = {1,2}; Line(12) = {2,3}; Line(13) = {3,4}; Line(14) = {4,1};\n"
    geo += "Curve Loop(21) = {11,12,13,14};\nPlane Surface(31) = {21};\n"

    # puntos de la polilínea (IDs a partir de 1001)
    ids = []
    for k, point in enumerate(poly_xy):
        x = point[0]
        y = point[1]
        geo += pt(1001+k, x, y)
        ids.append(1001+k)
    geo += f"Spline(2001) = {{{', '.join(map(str, ids))}}};\n"

    # campo de refinamiento
    geo += "Field[1] = Distance; Field[1].CurvesList = {2001}; Field[1].NumPointsPerCurve = 200;\n"
    geo += f"Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = {lc_path}; "
    geo += f"Field[2].SizeMax = {lc}; Field[2].DistMin = {rin}; Field[2].DistMax = {rout};\n"
    geo += f"Background Field = 2; Mesh.MeshSizeMin = {lc_path}; Mesh.MeshSizeMax = {lc};\n"
    with open(fname, "w") as f: f.write(geo)

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
t_max = 0.3
dt = 2.5e-4
fname = open(f"{CASE_DIR}/output.csv", 'w')
fname.write("time,pressure,volume,wplus,wminus\n")
store_freq = 10 # escribir vtk cada n-step
output_freq = 10 # Flush output cada n-step

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


f = phase.get()
poly = extract_polyline_from_isocontour(mesh, f)
box_domain = (0, -10, 10, 0.5)
write_geo_with_path(box_domain, poly, lc=0.4, lc_path=5e-3, rin=l_c/2, rout=2 - 3*l_c, fname=f"{MESH_NAME}_new.geo")

if fname:
    fname.close()
    logger.info("Closed output CSV file.") 