from fenics import *
from ufl import grad, inner, dot, tr, Identity
import numpy as np
import matplotlib.pyplot as plt
from material_model import ModuloPorCapas

# 1. Parámetros del material, capas y rotación
E_capa1 = 210e4
E_capa2 = 210e7
nu_constante = 0.3
espesor_capa1 = 0.02
espesor_capa2 = 0.04
e_total = espesor_capa1 + espesor_capa2

# ¡NUEVO PARÁMETRO! Ángulo de rotación en grados
# angulo_grados = 90.0

deformaciones = []
angulos = np.linspace(0, 90, 20)

# 2. Malla y espacio de funciones
mesh = RectangleMesh(Point(0,0), Point(1, 1), 100, 100)
V = VectorFunctionSpace(mesh, 'P', 2)

for angulo_grados in angulos:

    E = ModuloPorCapas(E1=E_capa1, E2=E_capa2, e1=espesor_capa1, e2=espesor_capa2, degree=0, angle=angulo_grados)

    nu = Constant(nu_constante)
    mu = E / (2*(1 + nu))   
    lmbda = E*nu / ((1 + nu)*(1 - 2*nu))

    def epsilon(u):
        return 0.5*(grad(u) + grad(u).T)

    def sigma(u):
        return lmbda*tr(epsilon(u))*Identity(2) + 2*mu*epsilon(u)

    # Definir un problema de ejemplo (viga en voladizo)
    u = TrialFunction(V)
    v = TestFunction(V)

    # Empotrar el lado izquierdo
    bc = DirichletBC(V.sub(0), Constant(0.0), lambda x, on_boundary: on_boundary and near(x[0], 0))

    # Carga puntual en la esquina superior derecha
    a = inner(sigma(u), epsilon(v))*dx

    # 1. Definir la presión
    p = Constant(-5e3)  # Valor de la presión (ajusta el signo según la dirección deseada)

    # 2. Definir la normal y la frontera derecha
    n = FacetNormal(mesh)

    # 3. Definir la condición para la frontera derecha (x = 1)
    def right_boundary(x, on_boundary):
        return on_boundary and near(x[0], 1.0)

    # 4. Medida de integración sobre la frontera derecha

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 1.0)
    Right().mark(boundaries, 1)
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    L = dot(p * n, v) * ds(1)


    u_sol = Function(V, name="Desplazamiento")
    solve(a == L, u_sol, bc) # Comentado para correr sin errores de carga puntual


    # plot(u_sol, mode='displacement')
    # plt.show()

    # 1. Calcular el tensor de tensiones en la solución
    sigma_val = sigma(u_sol)

    # 2. Definir la normal y la frontera izquierda
    n = FacetNormal(mesh)

    # 3. Marcar la frontera izquierda (x = 0)
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 0.0)
    left_boundary = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    Left().mark(left_boundary, 1)
    ds_left = Measure('ds', domain=mesh, subdomain_data=left_boundary)

    # 4. Calcular la fuerza de reacción integrando la tracción sobre la frontera izquierda
    # Fuerza total (vectorial)
    reaction_force = assemble(dot(sigma_val * n, as_vector((1.0, 0.0))) * ds_left(1))
    print(f"Angulo: {angulo_grados}°, Fuerza de reacción en x=0: {reaction_force:.4e} N")
    despl_x = u_sol(1,0)[0]
    e_deform_compr = despl_x/1
    print("Desplazamientos: " , despl_x)
    print("Deformación Compresión: " , e_deform_compr)
    deformaciones.append(e_deform_compr)



plt.plot(angulos, deformaciones)
plt.show()