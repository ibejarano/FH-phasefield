# from dolfin import Function, VectorFunctionSpace, TrialFunction, TestFunction

import numpy as np
import ufl
from dolfinx import fem
from basix.ufl import element
from dolfinx.fem.petsc import LinearProblem

class DisplacementField:
    def __init__(self, mesh, V=None):
        if V is None:
            # Si no se proporciona un espacio de funciones, se crea uno por defecto
            u_cg2 = element("Lagrange", mesh.topology.cell_name(), 2, shape=(mesh.geometry.dim, ))
            self.V = fem.functionspace(mesh, u_cg2)
        self.new = fem.Function(self.V, name="displacement")
        self.old = fem.Function(self.V)
        self.temp = fem.Function(self.V)

    def update(self):
        self.old.x.array[:] = self.new.x.array[:]

    def get_trialfunction(self):
        return ufl.TrialFunction(self.V)

    def get_testfunction(self):
        return ufl.TestFunction(self.V)

    def get(self):
        return self.new
    
    def setup_solver(self, E_du, bc_u):
        """
        Configura el Solver para el campo de desplazamiento.
        """
        a = fem.form(ufl.lhs(E_du))
        L = fem.form(ufl.rhs(E_du))
        problem = LinearProblem(a, L, bc_u, self.new,
                                petsc_options={
                                        "ksp_type": "gmres",
                                        "pc_type": "hypre",
                                        "pc_hypre_type": "boomeramg",
                                        "ksp_rtol": 1e-8,
                                        "ksp_max_it": 1000
                                    })
        self.problem = problem

    def solve(self):
        """
        Resuelve el problema de desplazamiento.
        """
        self.problem.solve()
        self.update()