import ufl
import numpy as np
from dolfinx.fem import assemble_scalar, Constant
from dolfinx.fem import Function
from ufl import grad, sym, inner, tr, Identity, det, ln, dev, conditional, gt, dx


def epsilon(u):
    return sym(grad(u))

def sigma_hyperelastic(u, mu, lmbda):
    I = Identity(len(u))
    F = I + grad(u)
    J = det(F)
    F_invT = ufl.inv(F).T
    return mu * (F * F.T - I) + lmbda * ln(J) * F_invT

def sigma_linear(u, E, nu):
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return 2.0 * mu * epsilon(u) + lmbda * tr(epsilon(u)) * Identity(len(u))

def psi_linear(u, E, nu):
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return 0.5*(lmbda+mu)*(0.5*(tr(epsilon(u)) + abs(tr(epsilon(u)))))**2 + mu*inner(dev(epsilon(u)), dev(epsilon(u)))

def psi_linear_epsilon(eps, E, nu):
    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    return 0.5*(lmbda+mu)*(0.5*(tr(eps) + abs(tr(eps))))**2 + mu*inner(dev(eps), dev(eps))

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
    v = ufl.TestFunction(phi.function_space)
    expr = ufl.inner(grad(phi), -u) * v * dx
    return fem.assemble_scalar(fem.form(expr))

def get_E_expression(data):
    # Devuelve una expresión UFL para E(x) según regiones
    if "E_regions" in data:
        regions = sorted(data["E_regions"], key=lambda r: r["x_max"])
        expr = None
        prev = None
        x = ufl.SpatialCoordinate(data["mesh"])
        for i, region in enumerate(regions):
            cond = (x[0] <= region['x_max']) if i == 0 else ufl.And(x[0] > prev, x[0] <= region['x_max'])
            region_expr = ufl.conditional(cond, region['E'], 0.0)
            expr = region_expr if expr is None else expr + region_expr
            prev = region['x_max']
        # Para x > último x_max
        expr = expr + ufl.conditional(x[0] > regions[-1]['x_max'], regions[-1]['E'], 0.0)
        return expr
    else:
        return Constant(data["mesh"], float(data["E"]))
    

if __name__ == "__main__":
    # Ejemplo de uso
    from dolfinx import mesh, fem
    from mpi4py import MPI
    
    from mesh_setup import setup_gmsh
    from basix.ufl import element

    # Crear una malla de ejemplo
    domain, boundary_markers = setup_gmsh(
        case_dir="./results/fenicsx_tests",
        data = {
        "mesh_data": {
            "file_dir": "meshes",
            "file_name": "demo_fenicsx"
        }
    }
    )

    u_cg2 = element("CG", domain.topology.cell_name(), 2, shape=(domain.geometry.dim, ))
    phi_cg1 = element("CG", domain.topology.cell_name(), 1)

    V = fem.functionspace(domain, phi_cg1)
    W = fem.functionspace(domain, u_cg2)

    # Crear un campo de desplazamiento ficticio
    u = fem.Function(W)
    u.interpolate(lambda x: np.vstack((0.1 * x[0], 0.2 * x[1])))

    # Crear un campo de fase ficticio
    phi = fem.Function(V)
    phi.interpolate(lambda x: 0.5 + 0.5 * x[0])
    vol = compute_fracture_volume(phi, u)

    # Datos de ejemplo
    data = {
        "boundary_conditions": {
        "displacement": {
            "20": {"value": [0.0, 0.0]}
        },
        "initial_crack" : {
            "center": [0.0, 0.0],
            "l0": 0.0075,
            "w0": 0.5e-3
        }
    },
        "linit": 0.1,
        "h": 0.05
    }

        # Probar compute_fracture_volume
    vol = compute_fracture_volume(phi, u)
    if domain.comm.rank == 0:
        print(f"Volumen de fractura (ficticio): {vol:.6e}")


        msh = mesh.create_unit_square(MPI.COMM_WORLD, 8, 8)
    V = fem.functionspace(domain, phi_cg1)
    W = fem.functionspace(domain, u_cg2)

    # Crear un campo de desplazamiento ficticio
    u = fem.Function(W)
    u.interpolate(lambda x: np.vstack((0.1 * x[0], 0.2 * x[1])))

    # Crear un campo de fase ficticio
    phi = fem.Function(V)
    phi.interpolate(lambda x: 0.5 + 0.5 * x[0])

    # Probar compute_fracture_volume
    vol = compute_fracture_volume(phi, u)
    if msh.comm.rank == 0:
        print(f"Volumen de fractura (ficticio): {vol:.6e}")

