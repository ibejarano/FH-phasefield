#from dolfin import FunctionSpace, Function, project, conditional, gt
import ufl
import numpy as np
from dolfinx import fem
from basix.ufl import element

class HistoryField:
    def __init__(self, mesh, psi_func, E_expr, nu, data):
        """
        V: espacio de funciones
        psi_func: función de energía (por ejemplo, psi_linear)
        E_expr: módulo de Young espacial (Expression)
        nu: coeficiente de Poisson
        data: diccionario de configuración
        """
        element_type = element("DG", mesh.topology.cell_name(), 0)
        self.V = fem.functionspace(mesh, element_type)
        self.field = fem.Function(self.V)
        self.psi_func = psi_func
        self.E_expr = E_expr
        self.nu = nu
        self.data = data

    def update(self, u):
        """
        Actualiza el campo de historia con el desplazamiento actual u.
        """
        psi_val = self.psi_func(u, self.E_expr, self.nu)
        # Actualiza el campo de historia: H = max(H, psi)
        new_H_expr = ufl.conditional(ufl.gt(psi_val, self.field), psi_val, self.field)
        self.field.interpolate(fem.Expression(new_H_expr, self.V.element.interpolation_points()))

    def get(self):
        """
        Devuelve el campo de historia actual (Function).
        """
        return self.field