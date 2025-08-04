from dolfin import *
import numpy as np
import matplotlib.pyplot as plt


# Este primer caso esta enfocado en reproducir la placa con crack central
# La placa esta sometida a traccion 
# Calculamos el KI y KII

Lcrack = 0.2

mesh = Mesh("caso_2.xml")
boundaries = MeshFunction("size_t", mesh, "caso_2_facet_region.xml")

V = VectorFunctionSpace(mesh, "P", 1)

tol = 1E-2

def midpoint(x, on_boundary):
    return near(x[0], 0.0) and (x[0] > Lcrack) 

bc_bottom = DirichletBC(V.sub(0), Constant(0), boundaries, 2)
bc_top = DirichletBC(V.sub(0), Constant(0), boundaries, 4)
bc2 = DirichletBC(V.sub(0), Constant(0), boundaries, 10)
bc3 = DirichletBC(V.sub(0), Constant(0), boundaries, 11)
bcs = [bc_bottom, bc_top]
E = 164.3e9
nu = 0.32

mu = E / (2 * (1 + nu))

lmbda = E*nu / ((1 + nu)*(1 - 2*nu))
#lmbda = 2 * mu * lmbda / (lmbda + 2 * mu)

def epsilon(u):
    return sym(grad(u))

def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(2) + 2*mu*epsilon(u)


ds = ds(subdomain_data=boundaries)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(sigma(u), epsilon(v))*dx

p1 = 183e6*0.02

upper_traction = Constant((0.0, p1))
L_form  = dot(upper_traction, v)*ds(4) - dot(upper_traction, v)*ds(2)
u_sol = Function(V, name="desplazamientos")
solve(a == L_form, u_sol, bcs)


npoints = 200

xs = np.linspace(Lcrack*0.6, Lcrack*0.999, npoints)
uplus_res = np.zeros((npoints, 2))
uminus_res = np.zeros((npoints, 2))
d_offset = 1e-4

for i, x in enumerate(xs):
    uplus_res[i] = u_sol(x, d_offset)
    uminus_res[i] = u_sol(x, -d_offset)

dU = np.abs(uplus_res[:, 1] - uminus_res[:, 1]) # No TOCAR
dV = np.abs(uplus_res[:, 0] - uminus_res[:, 0])


#plt.plot(xs, uplus_res[:, 0], "r--", label="top face")
#plt.plot(xs, uminus_res[:, 0], "b--", label="bottom face")
# plt.semilogx(xs, dV, "k--")

kappa = 3 - 4 *nu # Plane - strain 
kappa = (3 - nu)/(1+nu) # Plane - stress
factor = np.sqrt(2 * np.pi) * mu / (1+kappa)

r = Lcrack-xs

KI_est = np.zeros(npoints)
KII_est = np.zeros(npoints)

print("xs, r , DU")
for i, r_x in enumerate(r):
    #print(xs[i], r_x, dU[i])
    KI_est[i] = factor * np.sqrt(1/r_x) * dU[i]
    KII_est[i] = factor * np.sqrt(1/r_x) * dV[i]


den_KI = (p1 * np.sqrt(np.pi * Lcrack))
norm_KI = KI_est / den_KI

KI_eq25 = 1.12 * p1 * np.sqrt(np.pi * Lcrack)
print("Teo eq 25" , KI_eq25 / den_KI)
print("Calculado", np.max(norm_KI))
plt.semilogx(r, norm_KI)
# plt.plot(r, KII_est)

plt.show()
file = File('caso2.pvd')
file << u_sol
# np.savetxt(f"caso_1.csv", np.array([r, KI_est, KII_est]), header="r,KI,KII", delimiter=",", comments='')