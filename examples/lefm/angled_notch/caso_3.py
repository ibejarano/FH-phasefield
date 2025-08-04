from dolfin import *
import numpy as np
from math import sin, cos
from variational_forms.linear_static import elastic_energy_funcional

beta = 45
alpha = 0
Lcrack = 0.1

mesh = Mesh("caso_3_beta.xml")
boundaries = MeshFunction("size_t", mesh, "caso_3_beta_facet_region.xml")

V = VectorFunctionSpace(mesh, "P", 1)

tol = 1E-5


bc_bottom = DirichletBC(V.sub(0), Constant(0), boundaries, 2)
bc_top = DirichletBC(V.sub(0), Constant(0), boundaries, 4)

bcs = [bc_bottom, bc_top]
E = 164.3e9
nu = 0.32

mu = E / (2 * (1 + nu))

lmbda = E*nu / ((1 + nu)*(1 - 2*nu))
# lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)

ds = ds(subdomain_data=boundaries)
u = TrialFunction(V)
v = TestFunction(V)

a = elastic_energy_funcional(u, v, lmbda, mu)

p1 = 183e6*0.02

upper_traction = Constant((0.0, p1))
L_form  = dot(upper_traction, v)*ds(4) - dot(upper_traction, v)*ds(2)
u_sol = Function(V, name="desplazamientos")
solve(a == L_form, u_sol, bcs)


# INICIO DEL POSTPROCESO
kappa = 3 - 4 *nu # Plane - strain 
# kappa = (3 - nu)/(1+nu) # Plane - stress
factor = np.sqrt(2 * np.pi) * mu / (1+kappa)
rel_KI = (p1 * np.sqrt(np.pi * Lcrack))


npoints = 500

beta_rad = np.deg2rad(beta)
print("Beta", beta, beta_rad)
xs = np.linspace(Lcrack*0.5*np.cos(beta_rad), Lcrack*np.cos(beta_rad), npoints)
ys = np.linspace(Lcrack*0.5*np.sin(beta_rad), Lcrack*np.sin(beta_rad), npoints)

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

    r_crack[i] = Lcrack -  np.sqrt(x**2 + y**2)
    KI_calc[i] = factor * abs(Un - Vn)/np.sqrt(r_crack[i])
    KII_calc[i] = factor * abs(Ut - Vt)/np.sqrt(r_crack[i])

dU = np.abs(uplus_res[:, 1] - uminus_res[:, 1]) # No TOCAR
dV = np.abs(uplus_res[:, 0] - uminus_res[:, 0])


KI_teo = p1 * np.sqrt(np.pi * Lcrack) * (np.cos(beta_rad)**2 + alpha * np.sin(beta_rad)**2)
KII_teo = p1 * np.sqrt(np.pi * Lcrack) * (1-alpha) * np.sin(beta_rad) * np.cos(beta_rad)
print("Teo KI" , KI_teo / rel_KI)
print("-------")
print("Teo KII" , KII_teo / rel_KI)

file = File('caso3.pvd')
file << u_sol
# np.savetxt(f"caso_1.csv", np.array([r, KI_est, KII_est]), header="r,KI,KII", delimiter=",", comments='')