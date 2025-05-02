# === Módulo de configuración de malla y espacios funcionales ===
from dolfin import Mesh, MeshFunction, FunctionSpace, VectorFunctionSpace, Point, RectangleMesh, DOLFIN_EPS, near, SubDomain, File
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
    boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)

    class BoundaryLeft(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], p0[0], DOLFIN_EPS)

    class BoundaryRight(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], p1[0], DOLFIN_EPS)

    class BoundaryBottom(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], p0[1], DOLFIN_EPS)

    class BoundaryTop(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], p1[1], DOLFIN_EPS)

    left = BoundaryLeft()
    right = BoundaryRight()
    bottom = BoundaryBottom()
    top = BoundaryTop()

    left.mark(boundary_markers, 30)
    right.mark(boundary_markers, 10)
    bottom.mark(boundary_markers, 20)
    top.mark(boundary_markers, 40)

    return mesh, boundary_markers

def set_function_spaces(mesh):
    V = FunctionSpace(mesh, 'CG', 1)
    W = VectorFunctionSpace(mesh, 'CG', 1)
    WW = FunctionSpace(mesh, 'DG', 0)
    return V, W, WW
