import logging
from dolfin import TensorFunctionSpace, Function

from .config import Config
from .mesh_setup import setup_gmsh
from .material_model import epsilon, select_sigma, select_psi, get_E_expression
from .variational_forms import define_variational_forms
from .boundary_conditions import setup_shallow_bc, setup_deep_bc, create_markers
from .fields.history import HistoryField
from .fields.phase import PhaseField
from .fields.displacement import DisplacementField

logger = logging.getLogger(__name__)

class Problem:
    """
    Sets up the physical problem, including mesh, function spaces,
    material models, fields, and variational forms.
    """
    def __init__(self, config: Config):
        self.config = config
        logger.info("Setting up the simulation problem...")

        ## MESHING ##        
        self.mesh = setup_gmsh(self.config.case_dir, self.config.params)

        self.Vsig = TensorFunctionSpace(self.mesh, "DG", 0)

        # Material models
        self.E_expr = get_E_expression(self.config.params["material_parameters"])
        self.psi = select_psi(self.config.get("psi_model", "linear"))
        self.sigma = select_sigma(self.config.get("psi_model", "linear"))
        logger.info(f"Using energy model: {self.config.get('psi_model', 'linear')}")

        # Fields definitions
        self.history = HistoryField(self.mesh, self.psi, self.E_expr, self.config.params["material_parameters"]["nu"], self.config.params)
        self.displacement = DisplacementField(self.mesh)
        self.phase = PhaseField(self.mesh)

        # Boundary Conditions
        # TODO: Make BC selection more robust based on config
        symmetric = self.config.get("symmetric", False)
        if self.config.params["case_type"] == "shallow":
            self.bc_u, self.bc_phi = setup_shallow_bc(self.phase, self.displacement, self.config.params, symmetric)
        else:
            self.bc_u, self.bc_phi = setup_deep_bc(self.phase, self.displacement, self.config.params, symmetric)

        self.markers = create_markers(self.mesh)

        # Functions
        self.sigt = Function(self.Vsig, name="stress")

        # Variational Forms
        E_du, E_phi, self.pressure = define_variational_forms(
            epsilon, self.sigma, self.history.get(), self.phase, self.displacement,
            self.config.params, self.E_expr, self.markers
        )

        # Setup solvers within the field objects
        self.displacement.setup_solver(E_du, self.bc_u)
        self.phase.setup_solver(E_phi, self.bc_phi)

        logger.info("Problem setup complete.") 