import ufl
import logging
from dolfinx import fem
from basix.ufl import element
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI

logger = logging.getLogger(__name__)

class PhaseField:
    def __init__(self, mesh, V=None):
        self.mesh = mesh
        if V is None:
            # Si no se proporciona un espacio de funciones, se crea uno por defecto
            phi_cg1 = element("Lagrange", mesh.topology.cell_name(), 1)
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
        local_norm2 = float(diff @ diff)
        # Suma global de la norma cuadrada
        global_norm2 = self.mesh.comm.allreduce(local_norm2, op=MPI.SUM)
        return global_norm2**0.5
    
    def setup_solver(self, E_phi, bc_phi):
        """
        Configura el Solver para el campo de fase.
        """
        try:
            a = fem.form(ufl.lhs(E_phi))
            L = fem.form(ufl.rhs(E_phi))
            problem = LinearProblem(a, L, bc_phi, self.new,
                                    petsc_options={
                                        "ksp_type": "gmres",
                                        "pc_type": "hypre",
                                        "pc_hypre_type": "boomeramg",
                                        "ksp_rtol": 1e-8,
                                        "ksp_max_it": 1000
                                    })

            self.problem = problem
            logger.info("Phase field solver setup completed successfully")

        except Exception as e:
            logger.error(f"Error in phase field solver setup: {e}")
            raise

    def solve(self):
        """
        Resuelve el problema de campo de fase.
        """
        self.problem.solve()
        logger.debug("Phi solved")
        error = self.get_error()
        logger.debug(f"Error get {error:.2e}")
        self.update()
        logger.debug("Phi Updated")
        return error