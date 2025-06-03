from dolfinx import fem
from basix.ufl import element

class StressField:
    def __init__(self, mesh, data, E_expr):
        tensor_elem = element("DG", mesh.topology.cell_name(), 0, shape=(mesh.geometry.dim, mesh.geometry.dim))
        self.V = fem.functionspace(mesh, tensor_elem)
        self.sigmat = fem.Function(self.V, name="stress")
        self.E_expr = E_expr
        self.data = data

    def compute(self, sigma_func, pnew, unew):
        """Compute the stress field based on displacement, phase, and history."""
        # Assuming `select_sigma` is a function that computes the stress tensor
        nu = self.data.get("nu", 0.3)
        stress_expr = (1-pnew)**2 * sigma_func(unew, self.E_expr, nu)
        self.sigmat.interpolate(stress_expr)
        return self.sigmat