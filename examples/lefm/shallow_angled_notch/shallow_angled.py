from dolfin import *
import numpy as np
from math import sin, cos
import subprocess
from variational_forms.linear_static import elastic_energy_funcional

def reemplazar_H(ruta_geo: str, H_nuevo: float, beta_nuevo) -> None:
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
            else:
                f.write(linea)

np.seterr(divide='raise')


# Vamos a editar el geo antes de simular
geo_name = "shallow.geo"

H_prof = 0.1
beta = 15
alpha = 0
Lcrack = 0.05
E = 164.3e9
nu = 0.32
tol = 1E-5
esp = 0.02
p1 = 183e6*esp # Presion interna

# casos = [0.04, 0.05, 0.1, 0.15, 0.2, 0.5, 1, 4, 10]
casos = np.geomspace(0.025, 10, 30)
betas = [0, 5, 12, 25]

K_calcs = np.zeros((len(casos), 2))
up_arr = np.zeros(len(casos))
um_arr = np.zeros(len(casos))
us_arr = np.zeros(len(casos))

for beta in betas:
    for j, H_prof in enumerate(casos):

        reemplazar_H(geo_name, H_prof, beta)

        cmd_mallado = ["gmsh", "-2", "shallow.geo", "-format", "msh2", "-o", "shallow.msh"] 
        subprocess.run(cmd_mallado, check=True)
        cmd_xml_transformer = ["dolfin-convert" ,"shallow.msh", "shallow.xml"]
        subprocess.run(cmd_xml_transformer, check=True)

        mesh = Mesh("shallow.xml")
        boundaries = MeshFunction("size_t", mesh, "shallow_facet_region.xml")

        check_H_sup = np.max(mesh.coordinates()[:, 1])

        # assert np.allclosecheck_H_sup == H_prof, "No coinciden las profundidades"

        print(f"Simulando caso: {H_prof:.2f} mts. profundidad")

        V = VectorFunctionSpace(mesh, "P", 2)

        bc_bottom = DirichletBC(V.sub(1), Constant(0), boundaries, 2)
        bc_left = DirichletBC(V.sub(0), Constant(0), boundaries, 1)
        bc_right = DirichletBC(V.sub(0), Constant(0), boundaries, 3)

        bcs = [bc_bottom, bc_left, bc_right]


        mu = E / (2 * (1 + nu))

        lmbda = E*nu / ((1 + nu)*(1 - 2*nu))
        # lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)

        ds = ds(subdomain_data=boundaries)
        u = TrialFunction(V)
        v = TestFunction(V)

        a = elastic_energy_funcional(u, v, lmbda, mu)

        upper_traction = Constant((0.0, p1))
        L_form  = dot(upper_traction, v)*ds(10) - dot(upper_traction, v)*ds(11)
        u_sol = Function(V, name="desplazamientos")
        solve(a == L_form, u_sol, bcs)


        # INICIO DEL POSTPROCESO
        kappa = 3 - 4 *nu # Plane - strain 
        # kappa = (3 - nu)/(1+nu) # Plane - stress
        factor = np.sqrt(2 * np.pi) * mu / (1+kappa)
        rel_KI = (p1 * np.sqrt(np.pi * Lcrack))


        npoints = 600

        beta_rad = np.deg2rad(beta)
        xs = np.linspace(Lcrack*0.75*np.cos(beta_rad), Lcrack*np.cos(beta_rad), npoints)
        ys = np.linspace(Lcrack*0.75*np.sin(beta_rad), Lcrack*np.sin(beta_rad), npoints)

        uplus_res = np.zeros((npoints, 2))
        uminus_res = np.zeros((npoints, 2))
        KI_calc = np.zeros(npoints)
        KII_calc = np.zeros(npoints)

        r_crack = np.zeros(npoints)
        d_offset = 1e-6

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

        # El mismo calculo pero de otra manera

        KI_teo = p1 * np.sqrt(np.pi * Lcrack) * (np.cos(beta_rad)**2 + alpha * np.sin(beta_rad)**2)
        KII_teo = p1 * np.sqrt(np.pi * Lcrack) * (1-alpha) * np.sin(beta_rad) * np.cos(beta_rad)
        calc_KI = np.max(KI_calc)/rel_KI
        calc_KII = np.max(KII_calc)/rel_KI
        K_calcs[j] = [calc_KI, calc_KII]
        us_arr[j] = u_sol(0, check_H_sup)[1]
        up_arr[j] = u_sol(0, d_offset)[1]
        um_arr[j] = u_sol(0, -d_offset)[1]

    file = File('shallow.pvd')
    file << u_sol
    aH = Lcrack / np.array(casos)
    np.savetxt(f"outputs_beta{beta}.csv", np.array([aH, K_calcs[:, 0], K_calcs[:, 1], us_arr, up_arr,um_arr]).T, header="aH,KI,KII,us,up,um", delimiter=",", comments='')