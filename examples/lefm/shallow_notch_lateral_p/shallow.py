from dolfin import *
import numpy as np
from math import sin, cos
import subprocess
from variational_forms.linear_static import elastic_energy_funcional


def reemplazar_H(ruta_geo: str, H_nuevo: float, beta_nuevo: float, Lcrack: float) -> None:
    # Leer todas las líneas
    with open(ruta_geo, 'r', encoding='utf-8') as f:
        lineas = f.readlines()

    # Escribir de nuevo, reemplazando solo la línea que empieza con "H ="
    with open(ruta_geo, 'w', encoding='utf-8') as f:
        for linea in lineas:
            if linea.lstrip().startswith("H ="):
                # Conservamos el mismo formato: "H = <valor>;"
                f.write(f"H = {H_nuevo};\n")

            elif linea.lstrip().startswith("beta_angle ="):
                f.write(f"beta_angle = {beta_nuevo};\n")
            elif linea.lstrip().startswith("Lcrack ="):
                f.write(f"Lcrack = {Lcrack};\n")
            else:
                f.write(linea)


def run_shallow_case(H_prof: float, 
                     p1: float, 
                     pxx: float,
                     Lcrack: float,
                     E: float,
                     nu:float,
                     beta: float = 0,
                     geo_name="shallow.geo", 
                     save_vtu=False,
                     tol= 1e-5,
                     remesh=True):
    from dolfin import ds
    np.seterr(divide='raise')

    if remesh:
        reemplazar_H(geo_name, H_prof, beta, Lcrack)
        cmd_mallado = ["gmsh", "-2", "shallow.geo", "-format", "msh2", "-o", "shallow.msh"] 
        subprocess.run(cmd_mallado, check=True)
        cmd_xml_transformer = ["dolfin-convert" ,"shallow.msh", "shallow.xml"]
        subprocess.run(cmd_xml_transformer, check=True)

    mesh = Mesh("shallow.xml")
    boundaries = MeshFunction("size_t", mesh, "shallow_facet_region.xml")

    check_H_sup = np.max(mesh.coordinates()[:, 1])

    # assert np.allclosecheck_H_sup == H_prof, "No coinciden las profundidades"

    print(f"Simulando caso: {check_H_sup:.4f} (malla) {H_prof:.4f} (dato) mts. profundidad")

    V = VectorFunctionSpace(mesh, "P", 2)

    bc_bottom = DirichletBC(V.sub(1), Constant(0), boundaries, 2)
    bcs = [bc_bottom]


    mu = E / (2 * (1 + nu))

    lmbda = E*nu / ((1 + nu)*(1 - 2*nu))
    # lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)

    ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = elastic_energy_funcional(u, v, lmbda, mu)

    n = FacetNormal(mesh)
    lateral_compression = Constant((pxx, 0.0))
    L_form  = - dot(p1*n, v)*ds(10) - dot(p1*n, v)*ds(11)
    L_form += dot(lateral_compression, v)*ds(1) - dot(lateral_compression, v)*ds(3)
    u_sol = Function(V, name="desplazamientos")
    solve(a == L_form, u_sol, bcs)

    # INICIO DEL POSTPROCESO
    kappa = 3 - 4 *nu # Plane - strain 
    # kappa = (3 - nu)/(1+nu) # Plane - stress
    factor = np.sqrt(2 * np.pi) * mu / (1+kappa)
    rel_KI = (p1 * np.sqrt(np.pi * Lcrack))

    npoints = 500

    beta_rad = np.deg2rad(beta)
    xs = np.linspace(Lcrack*0.75*np.cos(beta_rad), Lcrack*np.cos(beta_rad), npoints)
    ys = np.linspace(Lcrack*0.75*np.sin(beta_rad), Lcrack*np.sin(beta_rad), npoints)

    uplus_res = np.zeros((npoints, 2))
    uminus_res = np.zeros((npoints, 2))
    KI_calc = np.zeros(npoints)
    KII_calc = np.zeros(npoints)

    r_crack = np.zeros(npoints)
    d_offset = tol
    if save_vtu:
        file = File('shallow_lateral.pvd')
        file << u_sol
    for i, (x, y) in enumerate(zip(xs, ys)):
        U = u_sol(x-d_offset, y+d_offset)
        Un = (U[1] * cos(beta_rad) - U[0] * sin(beta_rad))
        Ut = (U[0] * cos(beta_rad) + U[1] * sin(beta_rad))

        uplus_res[i] = [Ut, Un]

        V = u_sol(x+d_offset, y-d_offset)
        Vn = (V[1] * cos(beta_rad) - V[0] * sin(beta_rad))
        Vt = (V[0] * cos(beta_rad) + V[1] * sin(beta_rad))

        uminus_res[i] = [Vt, Vn]

        try:
            r_crack[i] = Lcrack -  np.sqrt(x**2 + y**2)
            KI_calc[i] = factor * abs(Un - Vn)/np.sqrt(r_crack[i])
            KII_calc[i] = factor * abs(Ut - Vt)/np.sqrt(r_crack[i])
        except FloatingPointError:
            r_crack[i] = 0
            KI_calc[i] = 0
            KII_calc[i] = 0

    calc_KI = np.max(KI_calc)/rel_KI
    calc_KII = np.max(KII_calc)/rel_KI
    #K_calcs[j] = [calc_KI, calc_KII]
    us = u_sol(0, check_H_sup)[1]
    up = u_sol(0, d_offset)[1]
    um = u_sol(0, -d_offset)[1]


    return [calc_KI, calc_KII], [us, up, um]




# Vamos a editar el geo antes de simular
H = 0.25
E = 164.3e9
nu = 0.32
esp = 0.02
p1 = 183e6*esp # Presion interna
beta = 12
ratios_casos = [0, 0.2, 0.5, 1, 2, 5]
# casos_H = np.geomspace(0.1, 10, 30)
a_cracks = np.geomspace(H/100, H*3, 20)

for i, ratio in enumerate(ratios_casos):
    K_calcs = np.zeros((len(a_cracks), 2))
    up_arr = np.zeros(len(a_cracks))
    um_arr = np.zeros(len(a_cracks))
    us_arr = np.zeros(len(a_cracks))
    pxx = ratio*p1
    for j, a_crack in enumerate(a_cracks):
        #try:
        [KI, KII], [us, up, um] = run_shallow_case(H, p1, pxx, a_crack, E=E, nu=nu, beta=beta, save_vtu=True, remesh=True)
        K_calcs[j] = [KI, KII]
        up_arr[j] = up
        um_arr[j] = um
        us_arr[j] = us
        #except:
        #print(f"Simulacion ratio {ratio} no converge.")

    aH = a_cracks / H
    np.savetxt(f"./beta_{beta}/ratio_{ratio}.csv", np.array([aH, K_calcs[:, 0], K_calcs[:, 1], us_arr, up_arr,um_arr]).T, header="aH,KI,KII,us,up,um", delimiter=",", comments='')