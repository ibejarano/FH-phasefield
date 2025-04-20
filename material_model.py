from dolfin import inv, sym, inner, tr, Identity, det, ln, dev, grad, conditional, gt

def epsilon(u):
    return sym(grad(u))

def sigma_hyperelastic(u, mu, lmbda):
    I = Identity(len(u))
    F = I + grad(u)
    J = det(F)
    F_invT = inv(F).T
    return mu * (F * F.T - I) + lmbda * ln(J) * F_invT

def sigma_linear(u, mu, lmbda):
    return 2.0 * mu * epsilon(u) + lmbda * tr(epsilon(u)) * Identity(len(u))

def psi_linear(u, data):
    mu = data["mu"]
    lmbda = data["lmbda"]
    return 0.5*(lmbda+mu)*(0.5*(tr(epsilon(u)) + abs(tr(epsilon(u)))))**2 + mu*inner(dev(epsilon(u)), dev(epsilon(u)))

def psi_hyperelastic(u, data):
    mu = data["mu"]
    lmbda = data["lmbda"]
    I = Identity(len(u))
    F = I + grad(u)
    C = F.T * F
    Ic = tr(C)
    J = det(F)
    return (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * ln(J)**2

def select_psi(model="linear"):
    if model == "linear":
        return psi_linear
    elif model == "hyperelastic":
        return psi_hyperelastic
    else:
        raise ValueError(f"Modelo psi desconocido: {model}")

def select_sigma(model="linear"):
    if model == "linear":
        return sigma_linear
    elif model == "hyperelastic":
        return sigma_hyperelastic
    else:
        raise ValueError(f"Modelo sigma desconocido: {model}")

def H(u, Hold, data, psi):
    delta_H = data.get("delta_H", 0.0)
    psi_u = psi(u, data)
    return conditional(gt(psi_u, Hold), psi_u, Hold + delta_H)