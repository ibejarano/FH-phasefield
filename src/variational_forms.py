from dolfin import inner, grad, dx, dot, Measure, Constant

def define_variational_forms(epsilon, sigma_func, H, phase, displacement, data, E_expr, markers):
    nu = data["material_parameters"]["nu"]
    Gc = data["material_parameters"]["Gc"]
    l = data["aspect_hl"] * data["meshing_parameters"]["h"]
    p_init = data.get("p_init", 100.0)
    pxx = data["px"]

    pressure = Constant(p_init)
    ds = Measure("ds", subdomain_data=markers)
    px_vec = Constant((pxx, 0.0))


    p = phase.get_trialfunction()
    q = phase.get_testfunction()
    u = displacement.get_trialfunction()
    v = displacement.get_testfunction()
    pold = phase.get_old()

    sigma = sigma_func(u, E_expr, nu)

    E_du = (1 - pold)**2 * inner(epsilon(v), sigma) * dx \
           + pressure * inner(v, grad(pold)) * dx \
           + dot(px_vec, v) * ds(20) - dot(px_vec, v) * ds(10)

    E_phi = (Gc * l * inner(grad(p), grad(q)) \
             + ((Gc / l) + 2.0 * H) * inner(p, q) \
             - 2.0 * H * q) * dx

    return E_du, E_phi, pressure


def fracture_phasefield_problem(epsilon, sigma_func, H, phase, displacement, data, E_expr, markers):
    nu = data["material_parameters"]["nu"]
    Gc = data["material_parameters"]["Gc"]
    l = data["aspect_hl"] * data["meshing_parameters"]["h"]
    p_init = data.get("p_init", 100.0)
    pxx = data["px"]

    pressure = Constant(p_init)
    ds = Measure("ds", subdomain_data=markers)
    px_vec = Constant((pxx, 0.0))


    p = phase.get_trialfunction()
    q = phase.get_testfunction()
    u = displacement.get_trialfunction()
    v = displacement.get_testfunction()
    pold = phase.get_old()

    sigma = sigma_func(u, E_expr, nu)

    E_du = (1 - pold)**2 * inner(epsilon(v), sigma) * dx \
           + pressure * inner(v, grad(pold)) * dx \
           + dot(px_vec, v) * ds(20) - dot(px_vec, v) * ds(10)

    E_phi = (Gc * l * inner(grad(p), grad(q)) \
             + ((Gc / l) + 2.0 * H) * inner(p, q) \
             - 2.0 * H * q) * dx

    return E_du, E_phi, pressure