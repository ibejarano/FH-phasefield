# from dolfin import DirichletBC, SubDomain, Constant
import numpy as np
from dolfinx.mesh import locate_entities
from basix.ufl import element
from dolfinx import default_scalar_type
from dolfinx import fem

def create_crack_domain(center, l0, w0):
    def crack_indicator(x):
        return np.logical_and(
            np.abs(x[0] - center[0]) <= l0,
            np.abs(x[1] - center[1]) <= w0
        )
    return crack_indicator

def setup_boundary_conditions(phase, displacement, mesh, boundary_markers, data):
    V = phase.V
    W = displacement.V
    bcs_u = []
    bcs_phi = []
    bc_config = data.get("boundary_conditions", {})

    # Desplazamiento
    for marker_id, bc_data in bc_config.get("displacement", {}).items():
        value = bc_data.get("value", [0.0, 0.0])
        marker_id_int = int(marker_id)
        # Encuentra las entidades (caras) con ese marcador
        entity_indices = boundary_markers.find(marker_id_int)
        for i, v in enumerate(value):
            if v != "None" and v is not None:
                dofs = fem.locate_dofs_topological(W.sub(i), mesh.topology.dim - 1, entity_indices)
                bc = fem.dirichletbc(default_scalar_type(0), dofs, W.sub(i))
                bcs_u.append(bc)

    # Fractura inicial
    crack_config = bc_config.get("initial_crack", {})

    def left_boundary(x):
        return np.isclose(x[0], 0.0)

    bc_phi_test = fem.dirichletbc(
        fem.Constant(mesh, 1.0),
        fem.locate_dofs_geometrical(V, left_boundary),
        V
    )

    if crack_config:
        center = crack_config.get("center", [0.0, 0.0])
        l0 = crack_config.get("l0", data.get("linit", 0.0))
        w0 = crack_config.get("w0", data.get("h", 0.0))
        crack_indicator = create_crack_domain(center, l0, w0)
        dofs = locate_entities(
            mesh, mesh.topology.dim,
            crack_indicator
        )
        bc_phi_crack = fem.dirichletbc(default_scalar_type(1.0), dofs, V)
        bcs_phi.append(bc_phi_crack)

    bcs_phi.append(bc_phi_test)
    return bcs_u, bcs_phi

if __name__ == "__main__":
    # Ejemplo de uso
    from dolfinx import mesh, fem
    from mpi4py import MPI
    
    from mesh_setup import setup_gmsh

    # Crear una malla de ejemplo
    mesh, boundary_markers = setup_gmsh(
        case_dir="./results/fenicsx_tests",
        data = {
        "mesh_data": {
            "file_dir": "meshes",
            "file_name": "demo_fenicsx"
        }
    }
    )

    u_cg2 = element("CG", mesh.topology.cell_name(), 2, shape=(mesh.geometry.dim, ))
    phi_cg1 = element("CG", mesh.topology.cell_name(), 1)

    V = fem.functionspace(mesh, phi_cg1)
    W = fem.functionspace(mesh, u_cg2)

    # Datos de ejemplo
    data = {
        "boundary_conditions": {
        "displacement": {
            "20": {"value": [0.0, 0.0]}
        },
        "initial_crack" : {
            "center": [0.0, 0.0],
            "l0": 0.0075,
            "w0": 0.5e-3
        }
    },
        "linit": 0.1,
        "h": 0.05
    }

    bcs_u, bcs_phi = setup_boundary_conditions(V, W, mesh, boundary_markers, data)
    print("Condiciones de frontera configuradas correctamente.")