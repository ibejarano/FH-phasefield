from dolfin import DirichletBC, SubDomain, Constant, near, MeshFunction

def create_crack_domain(center, l0, w0):
    class CrackDomain(SubDomain):
        def inside(self, x, on_boundary):
            # Usa center, l0, w0 pasados como argumentos
            return abs(x[0] - center[0]) <= l0 and abs(x[1] - center[1]) <= w0
    return CrackDomain()

def setup_shallow_bc(phase, displacement, data, symmetric=False):
    print("Setting up shallow bc", symmetric)
    V = phase.V
    W = displacement.V
    mesh = W.mesh()
    
    bcs_u = []
    bcs_phi = []

    def bottom_side(x, on_boundary):
        return near(x[1], mesh.coordinates()[:, 1].min())

    bc_bottom = DirichletBC(W, Constant((0.0, 0.0)), bottom_side)
    bcs_u.append(bc_bottom)

    if symmetric:
        bc_symmetry = setup_symmetric_bc(displacement)
        bcs_u.append(bc_symmetry)

    crack_config = data.get("initial_crack", {})
    if crack_config:
        center = crack_config.get("center", [0.0, 0.0])
        l0 = crack_config.get("l0", data.get("linit", 0.0)) # Usar linit si no está especificado aquí
        w0 = crack_config.get("w0", data.get("h", 0.0))     # Usar h si no está especificado aquí
        crack_subdomain = create_crack_domain(center, l0, w0)
        bc_phi_crack = DirichletBC(V, Constant(1.0), crack_subdomain)
        bcs_phi.append(bc_phi_crack)

    return bcs_u, bcs_phi

def setup_deep_bc(phase, displacement, data, symmetric=False):
    """
    Aqui voy a colocar como deep pero en realidad estoy restringuiendo el movimiento en la superficie libre
    """
    print("Setting up deep bc", symmetric)
    V = phase.V
    W = displacement.V
    mesh = W.mesh()

    bcs_u = []
    bcs_phi = []

    def upper_side(x, on_boundary):
        return near(x[1], mesh.coordinates()[:, 1].max())

    def bottom_side(x, on_boundary):
        return near(x[1], mesh.coordinates()[:, 1].min())

    bc_upper = DirichletBC(W, Constant((0.0, 0.0)), upper_side)
    bc_bottom = DirichletBC(W, Constant((0.0, 0.0)), bottom_side)
    bcs_u.append(bc_bottom)
    bcs_u.append(bc_upper)

    if symmetric:
        bc_symmetry = setup_symmetric_bc(displacement)
        bcs_u.append(bc_symmetry)

    crack_config = data.get("initial_crack", {})
    if crack_config:
        center = crack_config.get("center", [0.0, 0.0])
        l0 = crack_config.get("l0", data.get("linit", 0.0)) # Usar linit si no está especificado aquí
        w0 = crack_config.get("w0", data.get("h", 0.0))     # Usar h si no está especificado aquí
        crack_subdomain = create_crack_domain(center, l0, w0)
        bc_phi_crack = DirichletBC(V, Constant(1.0), crack_subdomain)
        bcs_phi.append(bc_phi_crack)

    return bcs_u, bcs_phi

def create_markers(mesh):
    markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    markers.set_all(0)

    class LeftBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], mesh.coordinates()[:, 0].min())
        
    class RightBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], mesh.coordinates()[:, 0].max())

    left = LeftBoundary()
    right = RightBoundary()

    left.mark(markers, 10)
    right.mark(markers, 20)
    return markers

def setup_symmetric_bc(displacement):
    """
    Configura condiciones de frontera para simetría en x=0 (u_x=0).
    """
    W = displacement.V

    def symmetry_plane(x, on_boundary):
        return on_boundary and near(x[0], 0.0)

    # u_x = 0 en x=0 (plano de simetría)
    bc_symmetry = DirichletBC(W.sub(0), Constant(0.0), symmetry_plane)
    return bc_symmetry