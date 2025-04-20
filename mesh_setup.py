# === Módulo de configuración de malla y espacios funcionales ===
def setup_mesh_and_spaces(mesh_path):
    from dolfin import Mesh, MeshFunction, FunctionSpace, VectorFunctionSpace

    mesh = Mesh(mesh_path + ".xml")
    subdomains = MeshFunction('size_t', mesh, mesh_path + "_physical_region.xml")
    boundary_markers = MeshFunction('size_t', mesh, mesh_path + "_facet_region.xml")

    V = FunctionSpace(mesh, 'CG', 1)
    W = VectorFunctionSpace(mesh, 'CG', 1)
    WW = FunctionSpace(mesh, 'DG', 0)

    return mesh, subdomains, boundary_markers, V, W, WW
