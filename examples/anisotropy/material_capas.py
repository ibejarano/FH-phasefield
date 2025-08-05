# Utilizar el modulo que define el modulo por capas
from material_model import ModuloPorCapas
from fenics import *
import matplotlib.pyplot as plt

# 1. Parámetros del material y las capas
E_capa1 = 210e7
E_capa2 = 210e4
nu_constante = 0.3
espesor_capa1 = 0.01
espesor_capa2 = 0.04

E = ModuloPorCapas(E1=E_capa1, E2=E_capa2, e1=espesor_capa1, e2=espesor_capa2, degree=0, angle=0)

# 2. Definir la malla y el espacio de funciones
# Un rectángulo de 1x2 para ver varias capas
mesh = RectangleMesh(Point(0,0), Point(1, 1), 200, 200)
V = VectorFunctionSpace(mesh, 'P', 1)

V_scalar = FunctionSpace(mesh, 'DG', 0) # Discontinuous Galerkin de grado 0

E_plot = project(E, V_scalar)
c = plot(E_plot, title="Módulo de Young por Capas")
plt.colorbar(c)
plt.show()