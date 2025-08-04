from dolfin import *
import numpy as np
from math import sin, cos
import matplotlib.pyplot as plt
import subprocess

def reemplazar_H(ruta_geo: str, nuevo_valor: float) -> None:
    # Leer todas las líneas
    with open(ruta_geo, 'r', encoding='utf-8') as f:
        lineas = f.readlines()

    # Escribir de nuevo, reemplazando solo la línea que empieza con "H ="
    with open(ruta_geo, 'w', encoding='utf-8') as f:
        for linea in lineas:
            if linea.lstrip().startswith("H ="):
                # Conservamos el mismo formato: "H = <valor>;"
                f.write(f"H = {nuevo_valor};\n")
            else:
                f.write(linea)

np.seterr(divide='raise')


# Vamos a editar el geo antes de simular
geo_name = "shallow.geo"

H_prof = 0.1
beta = 0
alpha = 0
Lcrack = 0.05
E = 164.3e9
nu = 0.32
tol = 1E-5
esp = 0.02
p1 = 183e6*esp # Presion interna

# casos = [0.04, 0.05, 0.1, 0.15, 0.2, 0.5, 1, 4, 10]
casos = np.geomspace(0.04, 10, 30)

K_calcs = np.zeros((len(casos), 2))
fig, axs = plt.subplots(2, 1, sharex=True)
up_arr = np.zeros(len(casos))
um_arr = np.zeros(len(casos))
us_arr = np.zeros(len(casos))

for j, H_prof in enumerate(casos):

    reemplazar_H(geo_name, H_prof)

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

    def epsilon(u):
        return sym(grad(u))

    def sigma(u):
        return lmbda*tr(epsilon(u))*Identity(2) + 2*mu*epsilon(u)


    ds = ds(subdomain_data=boundaries)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = inner(sigma(u), epsilon(v))*dx


    upper_traction = Constant((0.0, p1))
    L_form  = dot(upper_traction, v)*ds(10) - dot(upper_traction, v)*ds(11)
    u_sol = Function(V, name="desplazamientos")
    solve(a == L_form, u_sol, bcs)


    # INICIO DEL POSTPROCESO
    kappa = 3 - 4 *nu # Plane - strain 
    # kappa = (3 - nu)/(1+nu) # Plane - stress
    factor = np.sqrt(2 * np.pi) * mu / (1+kappa)
    rel_KI = (p1 * np.sqrt(np.pi * Lcrack))


    npoints = 200

    beta_rad = np.deg2rad(beta)
    xs = np.linspace(Lcrack*0.75*np.cos(beta_rad), Lcrack*np.cos(beta_rad), npoints)
    ys = np.linspace(Lcrack*0.75*np.sin(beta_rad), Lcrack*np.sin(beta_rad), npoints)

    uplus_res = np.zeros((npoints, 2))
    uminus_res = np.zeros((npoints, 2))
    KI_calc = np.zeros(npoints)
    KII_calc = np.zeros(npoints)

    r_crack = np.zeros(npoints)
    d_offset = 1e-5

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


    axs[0].semilogx(r_crack, KI_calc/rel_KI, label=f"{H_prof:.2f}")
    axs[1].semilogx(r_crack, KII_calc/rel_KI, label=f"{H_prof:.2f}")


    # El mismo calculo pero de otra manera

    KI_teo = p1 * np.sqrt(np.pi * Lcrack) * (np.cos(beta_rad)**2 + alpha * np.sin(beta_rad)**2)
    KII_teo = p1 * np.sqrt(np.pi * Lcrack) * (1-alpha) * np.sin(beta_rad) * np.cos(beta_rad)
    print("Teo KI" , KI_teo / rel_KI)
    calc_KI = np.max(KI_calc)/rel_KI
    print("Calc KI", calc_KI)
    print("-------")
    print("Teo KII" , KII_teo / rel_KI)
    calc_KII = np.max(KII_calc)/rel_KI
    print("Calc KII", calc_KII)
    K_calcs[j] = [calc_KI, calc_KII]
    us_arr[j] = u_sol(0, check_H_sup)[1]
    up_arr[j] = u_sol(0, d_offset)[1]
    um_arr[j] = u_sol(0, -d_offset)[1]

#axs[0].hlines(y = KI_teo/rel_KI, xmin=np.min(r_crack), xmax=np.max(r_crack), label="KI Teo")
#axs[1].hlines(y = KII_teo/rel_KI, xmin=np.min(r_crack), xmax=np.max(r_crack), label="KII Teo")
axs[0].set_ylabel("KI")
axs[1].set_ylabel("KII")
axs[1].set_xlabel("r (distancia al tip)")
axs[0].legend()
axs[1].legend()


plt.figure()

aH = Lcrack / np.array(casos)
plt.plot(aH, K_calcs[:, 0], marker="o", label="K I")
plt.plot(aH, K_calcs[:, 1], marker="o", label= "K II")
plt.xlabel("a / H")
plt.xscale("log")
plt.legend()
plt.show()
#file = File('shallow.pvd')
#file << u_sol
np.savetxt(f"outputs.csv", np.array([aH, K_calcs[:, 0], K_calcs[:, 1], us_arr, up_arr,um_arr]).T, header="aH,KI,KII,us,up,um", delimiter=",", comments='')