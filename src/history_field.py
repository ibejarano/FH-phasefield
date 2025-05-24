from dolfin import Function, project, conditional, gt

class HistoryField:
    def __init__(self, V, psi_func, E_expr, nu, data):
        """
        V: espacio de funciones
        psi_func: función de energía (por ejemplo, psi_linear)
        E_expr: módulo de Young espacial (Expression)
        nu: coeficiente de Poisson
        data: diccionario de configuración
        """
        self.field = Function(V)
        self.psi_func = psi_func
        self.E_expr = E_expr
        self.nu = nu
        self.data = data

    def update(self, u):
        """
        Actualiza el campo de historia con el desplazamiento actual u.
        """
        psi_val = self.psi_func(u, self.E_expr, self.nu)
        # Proyecta psi al espacio de funciones
        psi_proj = project(psi_val, self.field.function_space())
        # Actualiza el campo de historia: H = max(H, psi)
        new_H = conditional(gt(psi_proj, self.field), psi_proj, self.field)
        self.field.assign(project(new_H, self.field.function_space()))

    def get(self):
        """
        Devuelve el campo de historia actual (Function).
        """
        return self.field