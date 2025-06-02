#from dolfin import File, TimeSeries, XDMFFile
from dolfin.io import XDMFFile, TimeSeries

def create_output_files(case_dir, comm):
    xdmf = XDMFFile(comm, f"{case_dir}/output.xdmf", "w")
    xdmf.write_mesh = True
    return xdmf

def write_output(xdmf, u, phi, stress, t):
    xdmf.write_function(u, t)
    xdmf.write_function(phi, t)
    xdmf.write_function(stress, t)

def store_snapshot(u, phi, t, caseDir):
    # Guarda los campos como archivos XDMF individuales (opcional)
    with XDMFFile(u.function_space.mesh.comm, f"{caseDir}/u_{t:.4f}.xdmf", "w") as xdmf_u:
        xdmf_u.write_mesh(u.function_space.mesh)
        xdmf_u.write_function(u, t)
    with XDMFFile(phi.function_space.mesh.comm, f"{caseDir}/phi_{t:.4f}.xdmf", "w") as xdmf_phi:
        xdmf_phi.write_mesh(phi.function_space.mesh)
        xdmf_phi.write_function(phi, t)