from dolfin import inner, grad, dx, dot, Measure, Constant
from solvers import pressure_solver
from variational_forms.common import compute_fracture_volume, sigma, epsilon

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

def phase_field_problem(
        phase,
        displacement,
        history, 
        lmbda: float,
        mu: float,
        Gc: float,
        p_init: float,
        l_c: float
):
    pressure = Constant(p_init)
    p = phase.get_trialfunction()
    q = phase.get_testfunction()
    u = displacement.get_trialfunction()
    v = displacement.get_testfunction()
    pold = phase.get_old()

    sigma_u = sigma(u, lmbda, mu)

    H = history.get()

    E_du = (1 - pold)**2 * inner(epsilon(v), sigma_u) * dx \
           + pressure * inner(v, grad(pold)) * dx

    E_phi = (Gc * l_c * inner(grad(p), grad(q)) \
             + ((Gc / l_c) + 2.0 * H) * inner(p, q) \
             - 2.0 * H * q) * dx

    return E_du, E_phi, pressure

def compute_pressure(displacement, 
                     phase, 
                     history, 
                     pressure: float, 
                     vol_target: float, 
                     method="root_scalar"):
    
    ite_p, pn = pressure_solver(
        Vtarget=vol_target,
        phase=phase,
        displacement=displacement,
        history=history,
        pressure=pressure,
        vol_tol=1e-6,
        method=method
    )

    return ite_p, pn

def solve_step_staggered(displacement, phase, history, pressure, dV:float, phi_tol=1e-3):
    err_phi = 1.0
    outer_ite = 0
    V0 = compute_fracture_volume(phase.get_old(), displacement.get())
    vol_target = V0 + dV

    while err_phi > phi_tol:
        outer_ite += 1
        ite_p, pn = compute_pressure(displacement, phase, history, pressure, vol_target)
        if ite_p < 0:
            raise RuntimeError("Pressure adjustment failed to converge.")
        err_phi = phase.solve()
        if outer_ite > 15:
            raise RuntimeError(f"Outer staggered loop failed to converge (err_phi={err_phi:.2e})")

    displacement.update()
    phase.update()
    history.update(displacement.get())

    vol = compute_fracture_volume(phase.get(), displacement.get())

    return pn, vol