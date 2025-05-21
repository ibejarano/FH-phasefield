from dolfin import File, TimeSeries, XDMFFile

def create_output_files(caseDir):
    xdmf = XDMFFile(f"{caseDir}/output.xdmf")
    xdmf.parameters["flush_output"] = True
    xdmf.parameters["functions_share_mesh"] = True

    u_ts = TimeSeries(f"{caseDir}/u_series")
    phi_ts = TimeSeries(f"{caseDir}/phi_series")

    return xdmf, u_ts, phi_ts

def write_output(xdmf, u, phi, stress, t):
    xdmf.write(u, t)
    xdmf.write(phi, t)
    xdmf.write(stress, t)

def store_time_series(u_ts, phi_ts, u, phi, t):
    u_ts.store(u.vector(), t)
    phi_ts.store(phi.vector(), t)