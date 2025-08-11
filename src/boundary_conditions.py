from dolfin import DirichletBC, SubDomain, Constant, near, MeshFunction, Mesh

def create_crack_domain(center, l0, w0):
    class CrackDomain(SubDomain):
        def inside(self, x, on_boundary):
            # Usa center, l0, w0 pasados como argumentos
            return abs(x[0] - center[0]) <= l0 and abs(x[1] - center[1]) <= w0
    return CrackDomain()


def setup_bc(phase, 
                  displacement, 
                  linit: float,
                  h_elem: float,
                  crack_center= [0, 0],
                  upper_face_free = False,
                  symmetric = False):

    V = phase.V
    W = displacement.V
    mesh = W.mesh()

    bcs_u = []
    bcs_phi = []

    def upper_side(x, on_boundary):
        return near(x[1], mesh.coordinates()[:, 1].max())

    def bottom_side(x, on_boundary):
        return near(x[1], mesh.coordinates()[:, 1].min())

    bc_bottom = DirichletBC(W, Constant((0.0, 0.0)), bottom_side)
    bcs_u.append(bc_bottom)

    if not upper_face_free:
        bc_upper = DirichletBC(W, Constant((0.0, 0.0)), upper_side)
        bcs_u.append(bc_upper)

    if symmetric:
        bc_symmetry = setup_symmetric_bc(displacement)
        bcs_u.append(bc_symmetry)

    crack_subdomain = create_crack_domain(crack_center, linit, h_elem)
    bc_phi_crack = DirichletBC(V, Constant(1.0), crack_subdomain)
    # bcs_map = bc_phi_crack.get_boundary_values()

    # bcs_phi.append(bc_phi_crack)
    # dof_indices = sorted(bcs_map.keys())
    # print("Dofs afectados:", dof_indices)

    # coords = V.tabulate_dof_coordinates().reshape((-1, mesh.geometry().dim()))
    # boundary_node_coords = coords[dof_indices]

    # print("Coordenadas de los nodos con DirichletBC:")
    # for idx, xy in zip(dof_indices, boundary_node_coords):
    #     print(f"  dof {idx} → {xy}")


    bcs_phi.append(bc_phi_crack)
    return bcs_u, bcs_phi

def create_markers(mesh: Mesh, left_mark: int =10, right_mark: int =20):
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

    left.mark(markers, left_mark)
    right.mark(markers, right_mark)
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