# from dolfin import DirichletBC, SubDomain, Constant
import numpy as np
from dolfinx.mesh import locate_entities
from basix.ufl import element
from dolfinx import default_scalar_type
from dolfinx import fem
from dolfinx.fem import Constant, locate_dofs_geometrical
from petsc4py import PETSc

def create_crack_domain(center, l0, w0):
    def crack_indicator(x):
        return np.logical_and(
            np.abs(x[0] - center[0]) <= l0,
            np.abs(x[1] - center[1]) <= w0
        )
    return crack_indicator

def setup_displacement_bc(displacement, mesh, boundary_markers, data):
    W = displacement.V
    bcs_u = []
    bc_config = data.get("boundary_conditions", {}).get("displacement", {})

    # Desplazamiento
    for marker_id, bc_data in bc_config.items():
        value = bc_data.get("value", [0.0, 0.0])
        marker_id_int = int(marker_id)
        # Encuentra las entidades (caras) con ese marcador
        entity_indices = boundary_markers.find(marker_id_int)
        for i, v in enumerate(value):
            if v != "None" and v is not None:
                dofs = fem.locate_dofs_topological(W.sub(i), mesh.topology.dim - 1, entity_indices)
                bc = fem.dirichletbc(default_scalar_type(v), dofs, W.sub(i))
                bcs_u.append(bc)

    return bcs_u

def setup_phi_bc(phi, domain, data):
    V = phi.V
    bcs_phi = []
    bc_config = data.get("boundary_conditions", {})
    crack_config = bc_config.get("initial_crack", {})

    if crack_config:
        center = crack_config.get("center", [0.0, 0.0])
        l0 = crack_config.get("l0", data.get("linit", 0.0))
        w0 = crack_config.get("w0", data.get("h", 0.0))
        crack_indicator = create_crack_domain(center, l0, w0)
        dofs = locate_dofs_geometrical(V, crack_indicator        )
        bc_phi_crack = fem.dirichletbc(Constant(domain, PETSc.ScalarType(1)), dofs, V)
        bcs_phi.append(bc_phi_crack)

    return bcs_phi

def setup_boundary_conditions(phase, displacement, domain, boundary_markers, data):
    bcs_u = setup_displacement_bc(displacement, domain, boundary_markers, data)
    bcs_phi = setup_phi_bc(phase, domain, data)
    return bcs_u, bcs_phi

if __name__ == "__main__":
    # Ejemplo de uso
    from dolfinx import mesh, fem
    from mpi4py import MPI
    
    from mesh_setup import setup_gmsh
    from fields.displacement import DisplacementField
    from fields.phase import PhaseField

    # Crear una malla de ejemplo
    domain, boundary_markers = setup_gmsh(
        case_dir="./results/fenicsx_testing",
        data = {
        "mesh_data": {
            "file_dir": "meshes",
            "file_name": "demo_fenicsx"
        }
    }
    )

    u_cg2 = element("Lagrange", domain.topology.cell_name(), 2, shape=(domain.geometry.dim, ))
    phi_cg1 = element("Lagrange", domain.topology.cell_name(), 1)

    displacement = DisplacementField(domain)
    phase = PhaseField(domain)

    # Datos de ejemplo
    data = {
        "boundary_conditions": {
        "displacement": {
            "20": {"value": [0.0, 0.0]}
        },
        "initial_crack" : {
            "center": [0.0, 0.0],
            "l0": 0.01,
            "w0": 1e-3
        }
    },
        "linit": 0.1,
        "h": 0.05
    }

    bcs_u, bcs_phi = setup_boundary_conditions(phase, displacement, domain, boundary_markers, data)
    print("Condiciones de frontera configuradas correctamente.")

    for bc in bcs_phi:
        dofs = bc.dof_indices()[0]
        value = bc.g.value
        phase.new.x.array[dofs] = value

    # Exporta phi a VTK para inspección visual
    from dolfinx.io import VTKFile
    with VTKFile(domain.comm, "./phi_bc_test.pvd", "w") as vtk:
        vtk.write_function(phase.new)

    print("Campo de fase exportado a phi_bc_test.pvd. Abre este archivo en Paraview para verificar las condiciones de frontera.")