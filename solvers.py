from dolfin import lhs, rhs, LinearVariationalProblem, LinearVariationalSolver

def setup_solvers(E_du, E_phi, unew, pnew, bc_u, bc_phi):
    p_disp = LinearVariationalProblem(lhs(E_du), rhs(E_du), unew, bc_u)
    p_phi = LinearVariationalProblem(lhs(E_phi), rhs(E_phi), pnew, bc_phi)

    solver_disp = LinearVariationalSolver(p_disp)
    solver_phi = LinearVariationalSolver(p_phi)

    solver_phi.parameters["linear_solver"] = "gmres"
    solver_phi.parameters["preconditioner"] = "ilu"

    return solver_disp, solver_phi