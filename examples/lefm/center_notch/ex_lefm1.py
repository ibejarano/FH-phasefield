from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

Lx=2.0
Ly=2.0
Lcrack=0.3

mesh = Mesh("lefm.xml")

V = VectorFunctionSpace(mesh, "P", 2)

tol = 1E-14

def bottom_no_crack(x, on_boundary):
    return on_boundary and near(x[1], 0.0, tol) and np.greater(np.abs(x[0]), Lcrack)

def left(x, on_boundary):
    return on_boundary and near(x[0], 0, tol)

def right(x, on_boundary):
    return on_boundary and near(x[0], Lx, tol)

def top(x, on_boundary):
    return on_boundary and near(x[1], Ly, tol)

bc_bottom = DirichletBC(V.sub(1), Constant(0), bottom_no_crack)
bc_left = DirichletBC(V.sub(0), Constant(0), left)
bc_right = DirichletBC(V.sub(0), Constant(0), right)
bc_top = DirichletBC(V.sub(1), Constant(0), top)
bcs = [bc_bottom, bc_right, bc_top]

E = 2e8
nu = 0.3

mu = E / (2 * (1 + nu))

lmbda = E*nu / ((1 + nu)*(1 - 2*nu))

def epsilon(u):
    return sym(grad(u))

def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(2) + 2*mu*epsilon(u)

p0 = 1e5
p0 = 0.0
p1 = 1e2


n = FacetNormal(mesh)

class TopFace(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], Ly, tol)

class CrackFace(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and np.less(x[0], Lcrack)

crack = CrackFace()
top_bound = TopFace()

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)
crack.mark(boundaries, 1)
top_bound.mark(boundaries, 2)

ds = ds(subdomain_data=boundaries)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(sigma(u), epsilon(v))*dx

internal_pressure = Constant((0.0, p1))
L_form  =  dot(internal_pressure, v)*ds(1)
u_sol = Function(V)
solve(a == L_form, u_sol, bcs)

S = sigma(u_sol)

W = TensorFunctionSpace(mesh, 'DG', 0)

plot(u_sol*10, mode="displacement")
S_proj = project(S, W)

plt.figure()

## Evaluando el opening
npoints = 100
xs = np.linspace(0,Lcrack, npoints)
uys = np.zeros(npoints)
uxs = np.zeros(npoints)
for i, points in enumerate(xs):
    uys[i] = u_sol(points, 0.0)[1]
    uxs[i] = u_sol(points, 0.0)[0]

plt.plot(xs, uys, "-o")
#plt.show()

plt.figure()
plt.plot(xs, uxs, "-o")
plt.title("Desplazamientos X")

plt.figure()
r = (Lcrack-xs)

G = E / (2.0 * (1.0 + nu)) # Shear modulus
kappa = 3 - 4 *nu # Plane - strain 

factor = np.sqrt(2*np.pi) * G / (kappa + 1)

KI_est = np.zeros(npoints)
KII_est = np.zeros(npoints)

for i, r_x in enumerate(r):
    if r_x <= 0:
        continue
    KI_est[i] = uys[i] * np.sqrt(1/r_x)*factor
    KII_est[i] = abs(uxs[i]) * np.sqrt(1/r_x)*factor

plt.figure()

plt.semilogx(r, KI_est,".", label="KI")
plt.semilogx(r, KII_est,"--", label="KII")

plt.legend()

KI_an = abs(p1) * sqrt(np.pi * Lcrack)
print(KI_an)

plt.xlabel("r")
plt.ylabel(r"${u_y} \,\frac{2\mu}{k+1} \,\sqrt{2\pi/r}$")

plt.show()
