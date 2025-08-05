from dolfin import *
import numpy as np

# Parámetros del modelo
l = 0.01  # Longitud característica
Gc = 1.0  # Energía crítica de fractura
E = 1e4   # Módulo de Young
nu = 0.3  # Relación de Poisson
mu = E / (2 * (1 + nu))  # Módulo de cizalla
lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))  # Constante de Lamé

# Crear malla (por ejemplo, un cuadrado 2D)
L = 1.0  # Tamaño del dominio
mesh = RectangleMesh(Point(0, 0), Point(L, L), 100, 100)

# Espacios de funciones
el_u = FiniteElement("CG", mesh.ufl_cell(), 2)
el_phi = FiniteElement("CG", mesh.ufl_cell(), 1)

V = FunctionSpace(mesh, MixedElement([el_u, el_phi]))

w = Function(V)

(u, phi_n) = split(w)

print(type(phi_n))
phi = Function(phi_n)

(v, psi) = TestFunctions(V)

# Condición inicial para phi (fractura inicial opcional)
phi_init = Expression("exp(-pow(x[0] - 0.5, 2) / (2 * pow(0.01, 2)))", degree=2)
phi_n = interpolate(phi_init, w.sub(1))



# Tensor de deformación
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

# Energía elástica
def psi_0(eps):
    return 0.5 * lambda_ * tr(eps)**2 + mu * inner(eps, eps)

# Energía elástica degradada
kappa = 1e-6
psi = (1 - phi)**2 * psi_0(epsilon(u)) + kappa * psi_0(epsilon(u))

# Energía de fractura
gamma = Gc / 2 * (phi**2 / l + l * inner(grad(phi), grad(phi)))

# Energía total
F = psi * dx + gamma * dx

# Derivada variacional
dF = derivative(F, w, TestFunction(V))


p = Constant(1.0)  # Presión del fluido (puedes definirla como un campo)
F += -p * phi * div(u) * dx  # Término de acoplamiento hidráulico

# Condiciones de contorno (por ejemplo, desplazamiento fijo en el borde inferior)
bc_u = DirichletBC(V.sub(0), Constant((0, 0)), "near(x[1], 0.0)")
bc_phi = DirichletBC(V.sub(1), Constant(0.0), "near(x[0], 0.0) || near(x[0], 1.0)")
bcs = [bc_u, bc_phi]

# Resolver el problema
problem = NonlinearVariationalProblem(dF, w, bcs)
solver = NonlinearVariationalSolver(problem)
solver.solve()


# Calcular la energía de fractura
gamma_form = Gc / 2 * (phi**2 / l + l * inner(grad(phi), grad(phi))) * dx
Gamma = assemble(gamma_form)

# Longitud de la fractura
L_f = Gamma / Gc
print(f"Longitud de la fractura: {L_f:.4f}")


# Extraer región fracturada
phi_values = w.sub(1).vector().get_local()
fractured_region = phi_values > 0.9
num_fractured_points = np.sum(fractured_region)

# Aproximar la longitud basado en el área de la región fracturada
h = mesh.hmin()  # Tamaño característico de la malla
L_f_approx = num_fractured_points * h
print(f"Longitud aproximada de la fractura (geométrica): {L_f_approx:.4f}")


V_p = FunctionSpace(mesh, "P", 1)
p = Function(V_p)
q = TestFunction(V_p)

# Permeabilidad dependiente de phi
k = Constant(1e-12) * phi**2
F_p = phi * p * q * dx + dt * inner(k / mu_f * grad(p), grad(q)) * dx