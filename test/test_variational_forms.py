import pytest
import ufl
from dolfinx.fem import Constant
from dolfinx.mesh import create_unit_square
from dolfinx.fem import functionspace
from dolfinx.fem.function import Function
from mpi4py import MPI

from variational_forms import define_elastic_energy, define_fracture_energy, define_variational_forms

# Mock class to simulate the basic functionality needed for testing
class MockFunction:
    def __init__(self, mesh):
        V = functionspace(mesh, ("Lagrange", 1))
        self.V = V

    def get_trialfunction(self):
        return ufl.TrialFunction(self.V)

    def get_testfunction(self):
        return ufl.TestFunction(self.V)

    def get_old(self):
        return Function(self.V)

class MockVectorFunction:
    def __init__(self, mesh):
        V = functionspace(mesh, ("Lagrange", 1, (2,)))
        self.V = V

    def get_trialfunction(self):
        return ufl.TrialFunction(self.V)

    def get_testfunction(self):
        return ufl.TestFunction(self.V)

    def get_old(self):
        return Function(self.V)

# Fixtures for common test setup
@pytest.fixture
def basic_test_setup():
    mesh = create_unit_square(MPI.COMM_WORLD, 10, 10)

    displacement = MockVectorFunction(mesh)
    phase = MockFunction(mesh)

    def epsilon(v):
        return ufl.sym(ufl.grad(v))

    def sigma_func(u, E, nu):
        return E/(2*(1+nu)) * (ufl.grad(u) + ufl.grad(u).T)

    data = {
        "Gc": 1.0,
        "aspect_hl": 0.1,
        "h": 1.0,
        "p_init": 100.0,
        "px": 0.5
    }

    boundary_markers = mesh
    elastic_expr = Constant(mesh, 1.0)
    nu = 0.3
    history = 1.0

    return {
        "mesh": mesh,
        "displacement": displacement,
        "phase": phase,
        "epsilon": epsilon,
        "sigma_func": sigma_func,
        "data": data,
        "boundary_markers": boundary_markers,
        "elastic_expr": elastic_expr,
        "nu": nu,
        "history": history
    }

# Tests for define_variational_forms
def test_define_variational_forms_basic(basic_test_setup):
    setup = basic_test_setup

    el_energy_u, fr_energy_phi, pressure = define_variational_forms(
        setup["epsilon"],
        setup["sigma_func"],
        setup["history"],
        setup["phase"],
        setup["displacement"],
        setup["data"],
        setup["boundary_markers"],
        setup["elastic_expr"],
        setup["nu"]
    )

    assert isinstance(el_energy_u, ufl.form.Form)
    assert isinstance(fr_energy_phi, ufl.form.Form)
    assert isinstance(pressure, Constant)
    assert float(pressure) == 100.0

def test_define_variational_forms_default_p_init(basic_test_setup):
    setup = basic_test_setup
    setup["data"].pop("p_init")  # Remove p_init to test default value

    _, _, pressure = define_variational_forms(
        setup["epsilon"],
        setup["sigma_func"],
        setup["history"],
        setup["phase"],
        setup["displacement"],
        setup["data"],
        setup["boundary_markers"],
        setup["elastic_expr"],
        setup["nu"]
    )

    assert float(pressure) == 100.0  # Default value

# Tests for define_elastic_energy
def test_define_elastic_energy_with_phase(basic_test_setup):
    setup = basic_test_setup
    pressure = Constant(setup["mesh"], 100.0)

    el_energy = define_elastic_energy(
        setup["epsilon"],
        setup["sigma_func"],
        setup["displacement"],
        setup["data"],
        setup["boundary_markers"],
        setup["elastic_expr"],
        setup["nu"],
        setup["phase"],
        pressure
    )

    assert isinstance(el_energy, ufl.form.Form)

def test_define_elastic_energy_without_phase(basic_test_setup):
    setup = basic_test_setup

    el_energy = define_elastic_energy(
        setup["epsilon"],
        setup["sigma_func"],
        setup["displacement"],
        setup["data"],
        setup["boundary_markers"],
        setup["elastic_expr"],
        setup["nu"]
    )

    assert isinstance(el_energy, ufl.form.Form)

# Tests for define_fracture_energy
def test_define_fracture_energy_basic(basic_test_setup):
    setup = basic_test_setup

    fr_energy = define_fracture_energy(
        setup["history"],
        setup["phase"],
        setup["data"]
    )

    assert isinstance(fr_energy, ufl.form.Form)

def test_define_fracture_energy_different_params(basic_test_setup):
    setup = basic_test_setup
    setup["data"]["Gc"] = 2.0
    setup["data"]["aspect_hl"] = 0.2
    setup["data"]["h"] = 0.5

    fr_energy = define_fracture_energy(
        setup["history"],
        setup["phase"],
        setup["data"]
    )

    assert isinstance(fr_energy, ufl.form.Form)

# Error handling tests
def test_define_variational_forms_missing_required_data(basic_test_setup):
    setup = basic_test_setup
    del setup["data"]["Gc"]

    with pytest.raises(KeyError):
        define_variational_forms(
            setup["epsilon"],
            setup["sigma_func"],
            setup["history"],
            setup["phase"],
            setup["displacement"],
            setup["data"],
            setup["boundary_markers"],
            setup["elastic_expr"],
            setup["nu"]
        )

def test_define_elastic_energy_invalid_boundary_marker(basic_test_setup):
    setup = basic_test_setup
    setup["boundary_markers"] = None

    with pytest.raises(Exception):  # Adjust the exception type based on actual behavior
        define_elastic_energy(
            setup["epsilon"],
            setup["sigma_func"],
            setup["displacement"],
            setup["data"],
            setup["boundary_markers"],
            setup["elastic_expr"],
            setup["nu"],
            setup["phase"],
            Constant(setup["mesh"], 100.0)
        )