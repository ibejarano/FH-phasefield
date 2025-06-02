from dolfinx.mesh import create_mesh, CellType
from dolfinx.io import gmshio
from mpi4py import MPI
import gmsh
import os

def setup_gmsh(case_dir, data):
    mesh_dir = data["mesh_data"]["file_dir"]
    mesh_name = data["mesh_data"]["file_name"]
    geo_path = os.path.join(mesh_dir, mesh_name)
    geo_file = geo_path + ".geo"

    msh_path = os.path.join(case_dir, mesh_name)
    msh_file = msh_path + ".msh"

    # 1. Generar el .msh desde el .geo usando gmsh Python API
    gmsh.initialize()
    gmsh.open(geo_file)
    gmsh.model.mesh.generate(2)  # 2D, usa 3 para 3D
    gmsh.write(msh_file)
    gmsh.finalize()

    # 2. Leer el .msh con dolfinx
    mesh, cell_tags, facet_tags = gmshio.read_from_msh(msh_file, MPI.COMM_WORLD, gdim=2)
    boundary_markers = facet_tags

    return mesh, boundary_markers


def set_function_spaces(mesh):
    pass
    #V = FunctionSpace(mesh, 'CG', 1)
    #W = VectorFunctionSpace(mesh, 'CG', 1)
    #return V, W


if __name__ == "__main__":
    # Ejemplo de uso
    case_dir = "./results/fenicsx_tests"
    data = {
        "mesh_data": {
            "file_dir": "meshes",
            "file_name": "demo_fenicsx"
        }
    }
    
    mesh, boundary_markers = setup_gmsh(case_dir, data)
    #V, W = set_function_spaces(mesh)
    
    print("Malla y espacios funcionales configurados correctamente.")