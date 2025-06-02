import ufl
from dolfinx.fem import Constant

def define_variational_forms(epsilon, sigma_func, H, phase, displacement, data, boundary_markers, E_expr, nu):
    Gc = data["Gc"]
    l = data["aspect_hl"] * data["h"]
    p_init = data.get("p_init", 100.0)
    pxx = data["px"]

    pressure = Constant(displacement.function_space.mesh, p_init)
    ds = ufl.Measure("ds", domain=displacement.function_space.mesh, subdomain_data=boundary_markers)
    px_vec = Constant(displacement.function_space.mesh, (pxx, 0.0))

    p = phase.get_trialfunction()
    q = phase.get_testfunction()
    u = displacement.get_trialfunction()
    v = displacement.get_testfunction()
    pold = phase.get_old()

    sigma = sigma_func(u, E_expr, nu)

    E_du = (1 - pold)**2 * ufl.inner(epsilon(v), sigma) * ufl.dx \
           + pressure * ufl.inner(v, ufl.grad(pold)) * ufl.dx \
           + ufl.dot(px_vec, v) * ds(10) - ufl.dot(px_vec, v) * ds(30)

    E_phi = (Gc * l * ufl.inner(ufl.grad(p), ufl.grad(q))
             + ((Gc / l) + 2.0 * H) * ufl.inner(p, q)
             - 2.0 * H * q) * ufl.dx

    return E_du, E_phi, pressure