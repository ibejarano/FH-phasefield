from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import subprocess

Lx= 2.0
Ly= 2.0
Lcrack=0.3

mesh = Mesh("mesh.xml")
boundaries = MeshFunction("size_t", mesh, "mesh_facet_region.xml")

V = VectorFunctionSpace(mesh, "P", 1)

tol = 1E-12

bc_bottom = DirichletBC(V.sub(1), Constant(0), boundaries, 2)
bc_left = DirichletBC(V.sub(0), Constant(0), boundaries, 1)
bc_right = DirichletBC(V.sub(0), Constant(0), boundaries, 3)
bc_top = DirichletBC(V.sub(1), Constant(0), boundaries, 4)

bcs = [bc_bottom, bc_left, bc_right]

E = 2e8
nu = 0.3

mu = E / (2 * (1 + nu))

lmbda = E*nu / ((1 + nu)*(1 - 2*nu))

def epsilon(u):
    return sym(grad(u))

def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(2) + 2*mu*epsilon(u)

p1 = 1e4

ds = ds(subdomain_data=boundaries)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(sigma(u), epsilon(v))*dx

internal_pressure = Constant((0.0, p1))
L_form  = dot(internal_pressure, v)*ds(10) - dot(internal_pressure, v)*ds(11)
u_sol = Function(V)
solve(a == L_form, u_sol, bcs)


npoints = 100
xs = np.linspace(0, Lcrack, npoints)
uplus_res = np.zeros((npoints, 2))
uminus_res = np.zeros((npoints, 2))
d_offset = 1e-6

for i, x in enumerate(xs):
    uplus_res[i] = u_sol(x, d_offset)
    uminus_res[i] = u_sol(x, -d_offset)

# plt.plot(xs, uplus_res[:, 1]  - uminus_res[:, 1], "-o", label=f"Opening")
# plt.plot(xs, uminus_res[:, 0],"-o", label=f"Inferior")



#plt.legend()
#plt.show()


dU = np.abs(uplus_res[:, 1] - uminus_res[:, 1]) # No TOCAR
dV = np.abs(uplus_res[:, 0] - uminus_res[:, 0])

#plt.plot(xs, uplus_res[:, 0], "r--", label="top face")
#plt.plot(xs, uminus_res[:, 0], "b--", label="bottom face")
# plt.semilogx(xs, dV, "k--")

G = E / (2.0 * (1.0 + nu)) # Shear modulus
kappa = 3 - 4 *nu # Plane - strain 
factor = np.sqrt(2 * np.pi) * G / (1+kappa)

r = (Lcrack-xs)

KI_est = np.zeros(npoints)
KII_est = np.zeros(npoints)

for i, r_x in enumerate(r):
    if r_x <= 0:
        continue
    KI_est[i] = factor * np.sqrt(1/r_x) * dU[i]
    KII_est[i] = factor * np.sqrt(1/r_x) * dV[i]

file = File('poisson.pvd')
file << u_sol
np.savetxt(f"output_y_{int(Ly*100)}.csv", np.array([r, KI_est, KII_est]), header="r,KI,KII", delimiter=",", comments='')