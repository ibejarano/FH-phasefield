import numpy as np
import ufl
from dolfinx import fem, mesh, default_scalar_type
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem

def get_rotated_stiffness_matrix(E1, E2, nu12, G12, theta_deg):
    """
    Calcula la matriz de rigidez para tensión plana y la rota a un ángulo dado.
    """
    # Convertir ángulo a radianes
    theta = np.deg2rad(theta_deg)
    c = np.cos(theta)
    s = np.sin(theta)

    # Calcular nu21
    nu21 = nu12 * E2 / E1

    # Matriz de rigidez en el sistema de coordenadas del material (1-2)
    Q = np.zeros((3, 3))
    Q[0, 0] = E1 / (1 - nu12 * nu21)
    Q[0, 1] = nu12 * E2 / (1 - nu12 * nu21)
    Q[1, 0] = Q[0, 1]
    Q[1, 1] = E2 / (1 - nu12 * nu21)
    Q[2, 2] = G12

    # Matriz de transformación (de Reuter)
    T = np.array([
        [c**2, s**2, 2*s*c],
        [s**2, c**2, -2*s*c],
        [-s*c, s*c, c**2 - s**2]
    ])
    
    # Matriz de rigidez rotada en el sistema global (x-y)
    Q_bar = T.T @ Q @ T

    
    # Convertir a tensor de 4to orden para FEniCSx
    C = np.zeros((2, 2, 2, 2))
    C[0, 0, 0, 0] = Q_bar[0, 0]
    C[0, 0, 1, 1] = C[1, 1, 0, 0] = Q_bar[0, 1]
    C[1, 1, 1, 1] = Q_bar[1, 1]
    C[0, 1, 0, 1] = C[0, 1, 1, 0] = C[1, 0, 0, 1] = C[1, 0, 1, 0] = Q_bar[2, 2]
    
    return ufl.as_matrix(Q_bar)

def epsilon(u):
    """Calcula el tensor de deformaciones."""
    return ufl.sym(ufl.grad(u))

def sigma(u, C_mat):
    """Calcula el tensor de esfuerzos usando la matriz de rigidez de Voigt."""
    eps = epsilon(u)
    # Vector de deformaciones en notación de Voigt: [eps_xx, eps_yy, gamma_xy]
    # gamma_xy = 2 * eps_xy
    eps_vec = ufl.as_vector([eps[0, 0], eps[1, 1], 2 * eps[0, 1]])
    
    # Calcular el vector de esfuerzos: sigma_vec = C_mat * eps_vec
    sigma_vec = C_mat * eps_vec
    
    # Reconstruir el tensor de esfuerzos de 2x2 a partir del vector de Voigt
    return ufl.as_tensor([[sigma_vec[0], sigma_vec[2]],
                          [sigma_vec[2], sigma_vec[1]]])

def run_numerical_test(theta_deg, material_props):
    """
    Ejecuta el ensayo de tracción numérico para un ángulo dado.
    """
    E1, E2, nu12, G12 = material_props
    
    # 1. Malla y Espacio de Funciones
    domain = mesh.create_unit_square(MPI.COMM_WORLD, 32, 32)
    V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))
    
    # 2. Condiciones de Contorno
    L = 1.0  # Longitud del dominio
    delta = 0.0001  # Desplazamiento aplicado
    
    # Lado izquierdo (x=0) fijo en X
    def left_boundary(x):
        return np.isclose(x[0], 0)
    
    fdim = domain.topology.dim - 1
    facets_left = mesh.locate_entities_boundary(domain, fdim, left_boundary)    
    u_left = np.array([0, 0], dtype=default_scalar_type)
    bc_left = fem.dirichletbc(u_left, fem.locate_dofs_topological(V, fdim, facets_left), V)
    
    # Lado derecho (x=L) con desplazamiento aplicado en X
    def right_boundary(x):
        return np.isclose(x[0], L)

    u_right = np.array([delta, 0], dtype=np.float64)
    # Necesitamos aplicar el BC solo en la componente X
    facets_right = mesh.locate_entities_boundary(domain, fdim, right_boundary)    
    u_right = np.array([delta, 0], dtype=default_scalar_type)
    bc_right = fem.dirichletbc(u_right, fem.locate_dofs_topological(V, fdim, facets_right), V)


    bcs = [bc_left, bc_right]
    
    # 3. Definición del problema variacional
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    
    # Obtener matriz de rigidez rotada
    C = get_rotated_stiffness_matrix(E1, E2, nu12, G12, theta_deg)
    
    # Forma débil
    a = ufl.inner(sigma(u, C), epsilon(v)) * ufl.dx
    f = fem.Constant(domain, (0.0, 0.0))
    L_form = ufl.inner(f, v) * ufl.dx
    
    problem = fem.petsc.LinearProblem(a, L_form, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    uh = problem.solve()
    
    # 4. Post-procesamiento
    # Calcular esfuerzo promedio en X
    avg_stress_xx_form = fem.form(sigma(uh, C)[0, 0] * ufl.dx)
    total_stress_xx = fem.assemble_scalar(avg_stress_xx_form)
    avg_stress_xx = total_stress_xx / (L * L) # Dividir por el área
    
    # Deformación promedio en X es conocida
    avg_strain_xx = delta / L
    
    # Módulo de Young
    E_eff = avg_stress_xx / avg_strain_xx
    
    return E_eff

# --- Programa Principal ---

if __name__ == "__main__":
    # Propiedades de un compuesto típico de fibra de vidrio (unidireccional) en GPa
    material_props = {
        "E1": 38.6,   # Módulo a 0 grados
        "E2": 8.27,   # Módulo a 90 grados
        "nu12": 0.26, # Coeficiente de Poisson
        "G12": 0.01   # Módulo de cortante
    }
    
    # Lista de ángulos a probar
    angles = [0, 15,20, 30, 45,55, 60, 75, 90]
    print("Calculando Módulo de Young efectivo para diferentes ángulos:")
    import matplotlib.pyplot as plt
    for gs in [0.01, 0.1, 1, 10, 100, 1000]:
        E_effs = list() 
        material_props["G12"] = gs
        
        for angle in angles:
            E_effective = run_numerical_test(angle, list(material_props.values()))
            print(f"  Ángulo = {angle:2d}°  ->  E_efectivo = {E_effective:.2f} GPa")
            E_effs.append(E_effective)


        plt.plot(angles, E_effs, label=f"G12 = {gs} GPa")
    plt.legend()
    plt.show()
        
    # Verificación con la fórmula analítica (el caso donde el máximo podría no estar en 0)
    # Este material no tiene un G12 tan alto como para causar un máximo off-axis
    # pero el código lo detectaría si fuera el caso.