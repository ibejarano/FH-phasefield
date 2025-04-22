from dolfin import *
import json
from sys import argv
import os
from utils import read_data
from mesh_setup import setup_mesh_and_spaces
from material_model import epsilon, select_sigma, select_psi, H
from variational_forms import define_variational_forms
from solvers import setup_solvers
from output_utils import create_output_files, write_output, store_time_series

set_log_level(50)


def run_simulation(data, caseDir, meshName):
    # Mallado
    h_elem = data["h"] # TODO: Cambiar con el tamaño de la malla en zona de fractura
    aspect_hl = data["aspect_hl"] # aspect_hl = e = l/h
    l = aspect_hl*h_elem # Longitud de transición

    # Condiciones Iniciales
    l0 = data["linit"]
    w0 = h_elem
    p_init = 100

    # Control de simulacion
    TOL_PHI = 1e-3 # Tolerancia de phi
    TOL_VOL = 0.001 # 0.1% de tolerancia de volumen inyectado
    DT = data["dt"]
    T_FINAL = DT * 10000

    assert len(argv) == 3 , "Case name not found and mesh"
    caseDir = os.path.join("./results", argv[1])
    meshName = caseDir+"/"+argv[2]
    ## MESHING ##
    mesh, subdomains, boundary_markers, V, W, WW = setup_mesh_and_spaces(meshName)

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

    outfile = open(f"./{caseDir}/simulation_output.txt", 'w')
    outfile.write(" -- Algoritmo --- \n")
    outfile.write(" Algoritmo con control de volumen \n")
    outfile.write(json.dumps(data))
    outfile.close()

    step = 0
    solver_phi.solve()
    pold.assign(pnew)

    Q0 = data["Qo"]

    def adjust_pressure(solver_disp, pressure, unew, uold, pold, V0, Q0, DT, TOL_VOL, WW, psi, data):
        errDV = 1
        errDV1 = None
        errDV2 = None
        DV0 = DT * Q0
        pn = list(pressure.values())[0]
        ite = 0

        while abs(errDV) / DV0 > TOL_VOL:
            ite += 1
            DV0 = DT * Q0
            try:
                pn = pn1 - errDV1 * (pn1 - pn2)/(errDV1 - errDV2)
            except:
                pn *= 1.001

            pressure.assign(pn)
            solver_disp.solve()

            VK = assemble(inner(grad(pold), -unew) * dx)
            errDV = DV0 - (VK - V0)
            uold.assign(unew)
            Hold.assign(project(psi(unew, data), WW))

            try:
                pn2 = pn1
                errDV2 = errDV1
                pn1 = pn
                errDV1 = errDV
            except:
                pn1 = pn
                errDV1 = errDV

            if ite > 20:
                break

        return ite, errDV

    def solve_phase_field(solver_phi, pnew, pold, mesh):
        solver_phi.solve()
        err_phi = errornorm(pnew, pold, norm_type='l2', mesh=mesh)
        pold.assign(pnew)
        return err_phi

    while t <= T_FINAL:
        step += 1
        print(f"Step: {step}")
        ite = 0
        V0 = assemble( inner(grad(phit), -ut) * dx )
        errDV = 1
        errDV1 = None
        errDV2 = None
        DV0 = DT * Q0
        err_phi = 1
        while err_phi > TOL_PHI:
            ite, errDV = adjust_pressure(solver_disp, pressure, unew, uold, pold, V0, Q0, DT, TOL_VOL, WW, psi, data)
            if ite > 20:
                print("Simulation finished by iterations")
                break
            err_phi = solve_phase_field(solver_phi, pnew, pold, mesh)

        if ite > 20:
            print(" too much iterations ")
            break

        ut.assign(unew)
        phit.assign(pnew)

        vol_frac = assemble( inner(grad(phit), -ut) * dx )
        t += DT
        fname.write(str(t) + ",")
        fname.write(str(pn) + ",")
        fname.write(str(vol_frac) + "\n")

        print(f"Converge t: {t:.4f} dt: {DT:.2e} --- Iteraciones: {ite}")
        # Save files
        write_output(out_xml, ut, phit, t, step)
        store_time_series(u_ts, phi_ts, ut, phit, t)

    fname.close()
    

data = read_data("lab")
caseDir = os.path.join("./results/", argv[1])
meshName = caseDir+"/"+argv[2]
run_simulation(data, caseDir, meshName)