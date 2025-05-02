# === Módulo de configuración de malla y espacios funcionales ===
from dolfin import Mesh, MeshFunction, FunctionSpace, VectorFunctionSpace, Point, RectangleMesh
import numpy as np # Import numpy
from os import path

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
