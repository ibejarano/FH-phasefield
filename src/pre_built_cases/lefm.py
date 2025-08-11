import subprocess
from math import sin, cos
from dolfin import *
import numpy as np

from variational_forms.linear_static import elastic_energy_funcional

def reemplazar_H(ruta_geo: str,  Lcrack: float, H_nuevo: float = None, beta_nuevo: float = None) -> None:
    # Leer todas las líneas
    with open(ruta_geo, 'r', encoding='utf-8') as f:
        lineas = f.readlines()

    # Escribir de nuevo, reemplazando solo la línea que empieza con "H ="
    with open(ruta_geo, 'w', encoding='utf-8') as f:
        for linea in lineas:
            if linea.lstrip().startswith("H ="):
                # Conservamos el mismo formato: "H = <valor>;"
                if H_nuevo is not None:
                    f.write(f"H = {H_nuevo};\n")

            elif linea.lstrip().startswith("beta_angle ="):
                if beta_nuevo is not None:
                    f.write(f"beta_angle = {beta_nuevo};\n")
            elif linea.lstrip().startswith("Lcrack ="):
                f.write(f"Lcrack = {Lcrack};\n")
            else:
                f.write(linea)

def run_gmsh(mesh_name: str, Lcrack=None , H_prof = None, beta = None , mesh=True):
    if mesh:
        reemplazar_H(f"{mesh_name}.geo", Lcrack, H_prof, beta)
        cmd_mallado = ["gmsh", "-2", f"{mesh_name}.geo", "-format", "msh2", "-o", f"{mesh_name}.msh"] 
        subprocess.run(cmd_mallado, check=True)
        cmd_xml_transformer = ["dolfin-convert" ,f"{mesh_name}.msh", f"{mesh_name}.xml"]
        subprocess.run(cmd_xml_transformer, check=True)
    
    mesh = Mesh(f"{mesh_name}.xml")
    boundaries = MeshFunction("size_t", mesh, f"{mesh_name}_facet_region.xml")
    return mesh, boundaries


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
                     mesh=True):
    from dolfin import ds
    np.seterr(divide='raise')

    mesh_name = geo_name.split('.geo')[0]
    mesh, boundaries = run_gmsh(mesh_name, Lcrack=Lcrack, H_prof=H_prof, beta=beta, mesh=mesh)

    check_H_sup = np.max(mesh.coordinates()[:, 1])

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

def run_shallow_case_symm(H_prof: float, 
                     p1: float, 
                     pxx: float,
                     Lcrack: float,
                     E: float,
                     nu:float,
                     beta: float = 0,
                     geo_name="shallow.geo", 
                     save_vtu=False,
                     tol= 1e-5,
                     mesh=True):
    from dolfin import ds
    np.seterr(divide='raise')

    mesh_name = geo_name.split('.geo')[0]
    mesh, boundaries = run_gmsh(mesh_name, Lcrack=Lcrack, H_prof=H_prof, beta=beta, mesh=mesh)

    check_H_sup = np.max(mesh.coordinates()[:, 1])

    print(f"Simulando caso: {check_H_sup:.4f} (malla) {H_prof:.4f} (dato) mts. profundidad")

    V = VectorFunctionSpace(mesh, "P", 2)

    bc_axis_y = DirichletBC(V.sub(0), Constant(0), boundaries, 1) # 1 es el tag de la cara izquierda
    bc_bottom = DirichletBC(V.sub(1), Constant(0), boundaries, 2)
    bcs = [bc_bottom, bc_axis_y]


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
    L_form -= dot(lateral_compression, v)*ds(3)
    u_sol = Function(V, name="desplazamientos")
    solve(a == L_form, u_sol, bcs)

    # INICIO DEL POSTPROCESO
    kappa = 3 - 4 *nu # Plane - strain 
    # kappa = (3 - nu)/(1+nu) # Plane - stress
    factor = np.sqrt(2 * np.pi) * mu / (1+kappa)
    rel_KI = (p1 * np.sqrt(np.pi * Lcrack))

    npoints = 800

    beta_rad = np.deg2rad(beta)
    xs = np.linspace(Lcrack*0.9*np.cos(beta_rad), Lcrack*np.cos(beta_rad), npoints)
    ys = np.linspace(Lcrack*0.9*np.sin(beta_rad), Lcrack*np.sin(beta_rad), npoints)

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


def run_left_notch(
    p1: float, 
    Lcrack: float,
    E: float,
    nu:float,
    geo_name="caso_2.geo", 
    save_vtu=False,
    tol= 1e-5,
    mesh=True   
):
    

    mesh_name = geo_name.split('.geo')[0]
    mesh, boundaries = run_gmsh(mesh_name, Lcrack, mesh=mesh)

    V = VectorFunctionSpace(mesh, "P", 1)


    def midpoint(x, on_boundary):
        return near(x[0], 0.0) and (x[0] > Lcrack) 

    bc_bottom = DirichletBC(V.sub(0), Constant(0), boundaries, 2)
    bc_top = DirichletBC(V.sub(0), Constant(0), boundaries, 4)

    bcs = [bc_bottom, bc_top]

    mu = E / (2 * (1 + nu))

    lmbda = E*nu / ((1 + nu)*(1 - 2*nu))
    #lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)

    ds = Measure("ds", subdomain_data=boundaries)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = elastic_energy_funcional(u, v, lmbda, mu)


    upper_traction = Constant((0.0, p1))
    L_form  = dot(upper_traction, v)*ds(4) - dot(upper_traction, v)*ds(2)
    u_sol = Function(V, name="desplazamientos")
    solve(a == L_form, u_sol, bcs)


    npoints = 200

    xs = np.linspace(Lcrack*0.6, Lcrack*0.999, npoints)
    uplus_res = np.zeros((npoints, 2))
    uminus_res = np.zeros((npoints, 2))
    d_offset = tol

    for i, x in enumerate(xs):
        uplus_res[i] = u_sol(x, d_offset)
        uminus_res[i] = u_sol(x, -d_offset)

    dU = np.abs(uplus_res[:, 1] - uminus_res[:, 1]) # No TOCAR
    dV = np.abs(uplus_res[:, 0] - uminus_res[:, 0])

    kappa = 3 - 4 *nu # Plane - strain 
    kappa = (3 - nu)/(1+nu) # Plane - stress
    factor = np.sqrt(2 * np.pi) * mu / (1+kappa)

    r = Lcrack-xs

    KI_est = np.zeros(npoints)
    KII_est = np.zeros(npoints)

    for i, r_x in enumerate(r):
        KI_est[i] = factor * np.sqrt(1/r_x) * dU[i]
        KII_est[i] = factor * np.sqrt(1/r_x) * dV[i]

    if save_vtu:
        file = File('caso2.pvd')
        file << u_sol

    return r, KI_est, KII_est