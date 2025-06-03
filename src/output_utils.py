#from dolfin import File, TimeSeries, XDMFFile
from dolfinx.io import VTKFile
import logging

logger = logging.getLogger(__name__)

def create_output_files(case_dir, mesh):
    # Create VTK writer for the output
    vtk_file = VTKFile(mesh.comm, f"{case_dir}/results.pvd", "w")

    return vtk_file


def write_output(vtx_writer, t, phi=None, u=None, stress=None):
    """
    Writes the current state to output files.
    """
    try:
        # Write fields that are provided
        if phi is not None:
            vtx_writer.write_function(phi, t)

        if u is not None:
            vtx_writer.write_function(u, t)

        if stress is not None:
            vtx_writer.write_function(stress, t)

        vtx_writer.flush()  # Make sure data is written to disk

    except Exception as e:
        logger.error(f"Error writing output: {e}")


def store_snapshot(u, phi, t, caseDir):
    # Guarda los campos como archivos XDMF individuales (opcional)
    with XDMFFile(u.function_space.mesh.comm, f"{caseDir}/u_{t:.4f}.xdmf", "w") as xdmf_u:
        xdmf_u.write_mesh(u.function_space.mesh)
        xdmf_u.write_function(u, t)
    with XDMFFile(phi.function_space.mesh.comm, f"{caseDir}/phi_{t:.4f}.xdmf", "w") as xdmf_phi:
        xdmf_phi.write_mesh(phi.function_space.mesh)
        xdmf_phi.write_function(phi, t)