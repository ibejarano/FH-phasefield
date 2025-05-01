from dolfin import *
import json
from sys import argv
import os
from src.utils import read_data
from src.mesh_setup import setup_gmsh, setup_rect_mesh, set_function_spaces
from src.material_model import epsilon, select_sigma, select_psi, H, compute_fracture_volume
from src.variational_forms import define_variational_forms
from src.solvers import setup_solvers
from src.output_utils import create_output_files, write_output, store_time_series


def run_simulation(data):
    caseDir = os.path.join("./results/", argv[1])

    # Mallado
    h_elem = data["h"] # TODO: Cambiar con el tamaño de la malla en zona de fractura
    aspect_hl = data["aspect_hl"] # aspect_hl = e = l/h
    l = aspect_hl*h_elem # Longitud de transición

    # Condiciones Iniciales
    l0 = data["linit"]
    w0 = h_elem
    p_init = 100

    caseDir = os.path.join("./results", argv[1])

    ## MESHING ##
    mesh, boundary_markers = setup_gmsh(caseDir, data)
    V, W, WW = set_function_spaces(mesh)

    # 1 . MESH TIENE QUE SER DOMAIN
    # --- mesh esta importado, tiene que ser domain
    p, q = TrialFunction(V), TestFunction(V)
    u, v = TrialFunction(W), TestFunction(W)

    # Parametros de Lame (material isotropo)
    psi_model = data.get("psi_model", "linear")
    psi = select_psi(psi_model)
    sigma = select_sigma("linear")  # o "hyperelastic"
    print(f"Usando modelo de energía: {psi_model}")

    bcbottom  = DirichletBC(W, (0.0, 0.0), boundary_markers, 20)
    bc_u = [bcbottom]

    # Condicion de borde de la fractura, se marca la entalla inicial con el valor de 1
    class CrackDomain(SubDomain):
        def inside(self, x, on_boundary):
            center = [0, 0.0]
            return abs(x[0] - center[0]) <= l0 and abs(x[1] - center[1]) <= w0

    crack = CrackDomain()
    bc_phi = [DirichletBC(V, Constant(1.0), crack)]

    unew, uold, ut = Function(W), Function(W), Function(W, name="displacement")
    pnew, pold, Hold, phit = Function(V), Function(V), Function(V), Function(V, name="phi")

    # Funcional del desplazamiento eq (13)
    E_du, E_phi, pressure = define_variational_forms(W, V, epsilon, sigma, H, psi, pold, u, v, p, q, unew, Hold, data, boundary_markers)

    solver_disp, solver_phi = setup_solvers(E_du, E_phi, unew, pnew, bc_u, bc_phi)

    solver_disp.solve()
    solver_phi.solve()

    out_xml, u_ts, phi_ts = create_output_files(mesh, caseDir)
    fname = open(f"./{caseDir}/output.csv", 'w')

    t = 0
    pn = p_init
    pressure.assign(pn)

    outfile = open(f"./{caseDir}/parameters_used.json", 'w')
    outfile.write(json.dumps(data))
    outfile.close()

    step = 0
    solver_phi.solve()
    pold.assign(pnew)


    def adjust_pressure(solver_disp, pressure, unew, uold, pold, V0, WW, psi, data):
        errV = 1
        errV1 = None
        errV2 = None
        Q0 = data["Qo"]
        DT = data["dt"]
        Vtarget = V0 + DT * Q0
        pn = list(pressure.values())[0]
        ite = 0

        while abs(errV)/Vtarget > data["tolerances"]["volume"]:
            ite += 1
            try:
                pn = pn1 - errV1 * (pn1 - pn2)/(errV1 - errV2)
            except:
                pn *= 1.001

            pressure.assign(pn)
            solver_disp.solve()

            VK = compute_fracture_volume(pold, unew)
            errV = Vtarget - VK
            uold.assign(unew)
            Hold.assign(project(psi(unew, data), WW))

            try:
                pn2 = pn1
                errV2 = errV1
                pn1 = pn
                errV1 = errV
            except:
                pn1 = pn
                errV1 = errV

            if ite > 20:
                ite = -1
                break

        return ite, pn

    def solve_phase_field(solver_phi, pnew, pold, mesh):
        solver_phi.solve()
        err_phi = errornorm(pnew, pold, norm_type='l2', mesh=mesh)
        pold.assign(pnew)
        return err_phi

    while t <= data["t_max"]:
        step += 1
        V0 = compute_fracture_volume(phit, ut)
        err_phi = 1
        while err_phi > data["tolerances"]["phi"]:
            ite, pn = adjust_pressure(solver_disp, pressure, unew, uold, pold, V0, WW, psi, data)
            if ite < 0:
                print("*** Warning: Pressure adjust error ***")
                exit()
            err_phi = solve_phase_field(solver_phi, pnew, pold, mesh)

        ut.assign(unew)
        phit.assign(pnew)

        vol_frac = compute_fracture_volume(phit, ut)
        t += data["dt"]

        fname.write(str(t) + ",")
        fname.write(str(pn) + ",")
        fname.write(str(vol_frac) + "\n")

        print(f"Step: {step} - Converge t: {t:.4f} --- Ites: {ite}")
        # Save files
        if step % 2 == 0:
            write_output(out_xml, ut, phit, t, step)
            store_time_series(u_ts, phi_ts, ut, phit, t)

    fname.close()
    

set_log_level(50)
data = read_data("lab_gmsh")
run_simulation(data)