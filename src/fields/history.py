from dolfin import FunctionSpace, Function, project, conditional, gt
from variational_forms.phase_field import psi

class HistoryField:
    def __init__(self, mesh, _lambda, _mu):
        """
        V: espacio de funciones
        psi_func: función de energía (por ejemplo, psi_linear)
        E_expr: módulo de Young espacial (Expression)
        nu: coeficiente de Poisson
        data: diccionario de configuración
        """
        self.V = FunctionSpace(mesh, "DG", 0)
        self.field = Function(self.V)
        self._lambda = _lambda
        self._mu = _mu

    def update(self, u: Function):
        """
        Actualiza el campo de historia con el desplazamiento actual u.
        """
        psi_val = psi(u, self._lambda, self._mu)
        # Actualiza el campo de historia: H = max(H, psi)
        new_H = conditional(gt(psi_val, self.field), psi_val, self.field)
        self.field.assign(project(new_H, self.field.function_space()))

    def get(self):
        """
        Devuelve el campo de historia actual (Function).
        """
        return self.field