from dolfin import inner, grad, dx, dot, Measure, Constant

def define_variational_forms(W, V, epsilon, sigma, H, psi, pold, u, v, p, q, unew, Hold, data, boundary_markers):
    Gc = data["Gc"]
    l = data["aspect_hl"] * data["h"]
    mu = data["mu"]
    lmbda = data["lmbda"]
    Q0 = data["Qo"]
    p_init = data.get("p_init", 100.0)
    pxx = data["px"]

    pressure = Constant(p_init)
    ds = Measure("ds", subdomain_data=boundary_markers)
    px_vec = Constant((pxx, 0.0))

    E_du = (1 - pold)**2 * inner(epsilon(v), sigma(u, mu, lmbda)) * dx \
           + pressure * inner(v, grad(pold)) * dx \
           + dot(px_vec, v) * ds(10) - dot(px_vec, v) * ds(30)

    E_phi = (Gc * l * inner(grad(p), grad(q)) \
             + ((Gc / l) + 2.0 * H(unew, Hold, data, psi)) * inner(p, q) \
             - 2.0 * H(unew, Hold, data, psi) * q) * dx

    return E_du, E_phi, pressure