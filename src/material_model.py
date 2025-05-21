from dolfin import inv, sym, inner, tr, Identity, det, ln, dev, grad, conditional, gt, assemble, dx, TensorFunctionSpace, Expression

def epsilon(u):
    return sym(grad(u))

def sigma_hyperelastic(u, mu, lmbda):
    I = Identity(len(u))
    F = I + grad(u)
    J = det(F)
    F_invT = inv(F).T
    return mu * (F * F.T - I) + lmbda * ln(J) * F_invT

def sigma_linear(u, E, nu):
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return 2.0 * mu * epsilon(u) + lmbda * tr(epsilon(u)) * Identity(len(u))

def psi_linear(u, E, nu):
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
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

def H(Hold, data, psi):
    delta_H = data.get("delta_H", 0.0)
    return conditional(gt(psi, Hold), psi, Hold + delta_H)

def compute_fracture_volume(phi, u):
    vol_frac = assemble( inner(grad(phi), -u) * dx )
    return vol_frac

def get_E_expression(data):
    if "E_regions" in data:
        regions = sorted(data["E_regions"], key=lambda r: r["x_max"])
        expr = ""
        prev = None
        for i, region in enumerate(regions):
            if i == 0:
                expr += f"(x[0]<={region['x_max']})*{region['E']} + "
            else:
                prev = regions[i-1]['x_max']
                expr += f"((x[0]>{prev}) && (x[0]<={region['x_max']}))*{region['E']} + "
        expr += f"((x[0]>{regions[-1]['x_max']}))*{regions[-1]['E']}"

        return Expression(expr, degree=0)
    else:
        return Expression(str(data["E"]), degree=0)