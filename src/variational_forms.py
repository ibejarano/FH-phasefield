import ufl
from dolfinx.fem import Constant

def define_variational_forms(epsilon, sigma_func, history, phase, displacement, data, boundary_markers, elastic_expr, nu):
    p_init = data.get("p_init", 100.0)
    pressure = Constant(displacement.V.mesh, p_init)

    el_energy_u = define_elastic_energy(epsilon, sigma_func, displacement, data, boundary_markers, elastic_expr, nu, phase, pressure)
    fr_energy_phi = define_fracture_energy(history, phase, data)
    
    return el_energy_u, fr_energy_phi, pressure


def define_elastic_energy(epsilon, sigma_func, displacement, data, boundary_markers, elastic_expr, nu, phase=None, pressure=None):
    u = displacement.get_trialfunction()
    v = displacement.get_testfunction()
    sigma = sigma_func(u, elastic_expr, nu)

    if phase and pressure:
        if not boundary_markers:
            raise ValueError("Boundary markers not defined")
        pxx = data.get("px", 0.0)
        px_vec = Constant(displacement.V.mesh, (pxx, 0.0))
        ds = ufl.Measure("ds", domain=displacement.V.mesh, subdomain_data=boundary_markers)
        p_old = phase.get_old()
        el_energy_du = (1 - p_old)**2 * ufl.inner(epsilon(v), sigma) * ufl.dx \
                + pressure * ufl.inner(v, ufl.grad(p_old)) * ufl.dx \
                + ufl.dot(px_vec, v) * ds(10) - ufl.dot(px_vec, v) * ds(30)
    else:
        el_energy_du = ufl.inner(epsilon(v), sigma) * ufl.dx
    return el_energy_du


def define_fracture_energy(history, phase, data):
    Gc = data["Gc"]
    l = data["aspect_hl"] * data["h"]

    p = phase.get_trialfunction()
    q = phase.get_testfunction()

    fracture_energy = (Gc * l * ufl.inner(ufl.grad(p), ufl.grad(q))
             + ((Gc / l) + 2.0 * history) * ufl.inner(p, q)
             - 2.0 * history * q) * ufl.dx
    return fracture_energy