from dolfin import inner, grad, dx, dot, Measure, Constant

def define_variational_forms(epsilon, sigma_func, H, pold, u, v, p, q, data, boundary_markers, E_expr, nu):
    Gc = data["Gc"]
    l = data["aspect_hl"] * data["h"]
    p_init = data.get("p_init", 100.0)
    pxx = data["px"]

    pressure = Constant(p_init)
    ds = Measure("ds", subdomain_data=boundary_markers)
    px_vec = Constant((pxx, 0.0))

    sigma = sigma_func(u, E_expr, nu)

    E_du = (1 - pold)**2 * inner(epsilon(v), sigma) * dx \
           + pressure * inner(v, grad(pold)) * dx \
           + dot(px_vec, v) * ds(10) - dot(px_vec, v) * ds(30)

    E_phi = (Gc * l * inner(grad(p), grad(q)) \
             + ((Gc / l) + 2.0 * H) * inner(p, q) \
             - 2.0 * H * q) * dx

    return E_du, E_phi, pressure