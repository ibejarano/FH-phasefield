#from dolfin import Function, FunctionSpace, TrialFunction, TestFunction, errornorm
import ufl
from dolfinx import fem
from basix.ufl import element
from dolfinx.fem.petsc import LinearProblem

class PhaseField:
    def __init__(self, mesh, V=None):
        self.mesh = mesh
        if V is None:
            # Si no se proporciona un espacio de funciones, se crea uno por defecto
            phi_cg1 = element("CG", mesh.topology.cell_name(), 1)
            self.V = fem.functionspace(mesh, phi_cg1)
        else:
            self.V = V

        self.new = fem.Function(self.V, name="phi")
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
    
    def get_old(self):
        return self.old
    
    def get_error(self):
        diff = self.new.x.array - self.old.x.array
        return float((diff @ diff)**0.5)
    
    def setup_solver(self, E_phi, bc_phi):
        """
        Configura el Solver para el campo de fase.
        """
        a = fem.form(ufl.lhs(E_phi))
        L = fem.form(ufl.rhs(E_phi))
        print("BCs phi:", bc_phi)
        self.problem = LinearProblem(a, L, bc_phi, self.new,
                                 petsc_options={"ksp_type": "gmres", "pc_type": "ilu"}
                                 )

    def solve(self):
        """
        Resuelve el problema de campo de fase.
        """
        self.problem.solve()
        error = self.get_error()
        self.update()
        return error