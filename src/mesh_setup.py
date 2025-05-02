# === Módulo de configuración de malla y espacios funcionales ===
from dolfin import Mesh, MeshFunction, FunctionSpace, VectorFunctionSpace, Point, RectangleMesh, HDF5File, XDMFFile
from os import path
import meshio
import gmsh

def setup_gmsh(caseDir, data):
    mesh_file = data["mesh_data"]["xml_file"]
    mesh_path = path.join(caseDir, mesh_file)
    mesh = Mesh(mesh_path + ".xml")
    boundary_markers = MeshFunction('size_t', mesh, mesh_path + "_facet_region.xml")

    return mesh, boundary_markers

# === Módulo de configuración de malla y espacios funcionales ===
def setup_rect_mesh(data):
    mesh_data = data["mesh_data"]
    p0 = mesh_data["p0"]
    p1 = mesh_data["p1"]
    nelem = mesh_data["nelem"]
    mesh = RectangleMesh(Point(p0[0], p0[1]), Point(p1[0], p1[1]), nelem[0], nelem[1])

    return mesh

def set_function_spaces(mesh):
    V = FunctionSpace(mesh, 'CG', 1)
    W = VectorFunctionSpace(mesh, 'CG', 1)
    WW = FunctionSpace(mesh, 'DG', 0)
    return V, W, WW


def generate_gmsh(case_dir, mesh_name):

    # Inicializar Gmsh
    gmsh.initialize()

    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # Versión 2.2 es compatible con ASCII v2
    gmsh.option.setNumber("Mesh.Binary", 0)  # Asegura formato ASCII
    
    gmsh.open(f"./meshes/{mesh_name}.geo")
    gmsh.model.mesh.generate(dim=2)  # Cambia dim a 3 si es un mallado 3D
    gmsh.write(f"{case_dir}/{mesh_name}.msh")
    gmsh.finalize()