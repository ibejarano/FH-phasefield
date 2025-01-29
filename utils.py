from dolfin import Mesh, MeshEditor, FunctionSpace, interpolate, assemble, dx, File, TimeSeries, XDMFFile, UserExpression
import numpy as np

def compute_opening_grad(uvec, phivec, grad_phivec, lelem, lfracmax, tol=0.01):
    #grad_phi.assign(project(grad(phit), WW))
    w_overx = list()
    TOL_PHI = tol

    xsn = np.arange(-lelem, -lfracmax/2, -lelem)
    xsp = np.arange(0, lfracmax/2, lelem)
    xs = np.concatenate([np.flip(xsn), xsp])


    for x in xs:
        ys = np.arange(0, 5e-2, lelem)
        wx = 0.0
        for y in ys:
            u_y = uvec(x, y)
            phi_y = phivec(x, y)
            grad_phi_y = grad_phivec(x, y)

            if phi_y < TOL_PHI:
                break
            wx += -u_y[1] * grad_phi_y[1] * lelem
            
        w_overx.append(wx)


    return xs, np.array(w_overx)

def compute_opening_cutoff(uvec, phivec, lelem, lfracmax, phi_val=0.5):
    xsn = np.arange(-lelem, -lfracmax/2, -lelem)
    xsp = np.arange(0, lfracmax/2, lelem)
    xs = np.concatenate([np.flip(xsn), xsp])
    w_overx = list()

    for x in xs:
        ys = np.flip(np.arange(0, 5e-2, lelem))
        wx_value = 0.0
        for y in ys:
            u_y = uvec(x, y)
            phi_y = phivec(x, y)

            if phi_y >= phi_val:
                wx_value = u_y[1]
                break

        w_overx.append(wx_value)
    return xs, np.array(w_overx)


def line_integral(u, A, B, n, cutoff = 0.0): 
    '''Integrate u over segment [A, B] partitioned into n elements'''
    A = np.array(A)
    B = np.array(B)
    assert u.value_rank() == 0
    assert len(A) == len(B) > 1 and np.linalg.norm(A-B) > 0
    assert n > 0

    # Mesh line for integration
    mesh_points = [A + t*(B-A) for t in np.linspace(0, 1, n+1)]
    tdim, gdim = 1, len(A)
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, "interval", tdim, gdim)
    editor.init_vertices(n+1)
    editor.init_cells(n)

    for vi, v in enumerate(mesh_points): editor.add_vertex(vi, v)

    for ci in range(n): editor.add_cell(ci, np.array([ci, ci+1], dtype='uintp'))

    editor.close()

    # Setup funcion space
    elm = u.function_space().ufl_element()
    family = elm.family()
    degree = elm.degree()
    V = FunctionSpace(mesh, family, degree)
    v = interpolate(u, V)

    newarr = np.array([ v if v > cutoff else 0 for v in v.vector()[:] ])
    v.vector()[:] = newarr[:]

    return assemble(v*dx)

def createSave(mesh, caseDir, fileType, varName=None):
    if fileType.lower() == "xml":
        File(f"{caseDir}/malla.xml.gz") << mesh
        xdmffile = XDMFFile(f"{caseDir}/output.xdmf")
        xdmffile.parameters["flush_output"] = True
        xdmffile.parameters["functions_share_mesh"] = True
        u_ts = TimeSeries(f"{caseDir}/u_series")
        p_ts = TimeSeries(f"{caseDir}/phi_series")
        return xdmffile, u_ts, p_ts
    elif fileType.lower() == "vtu":
        if varName == None:
            varName = "out_variable"
        out_f = File(f"{caseDir}/{varName}.pvd")
        return out_f, None
    else:
        raise Exception("file Type unknown: only XML or VTU")

def saveValues(saveFile, timesrs, vec, t):
    if timesrs:
        saveFile.write(vec, t)
        timesrs.store(vec.vector(), t)
    else:
        saveFile << vec
    return

def get_values_overline_x(x1, x2, y=0.0, vec=None , npoints=100):
    tol = 0.0001 # avoid hitting points outside the domain
    xs = np.linspace(x1 + tol, x2 - tol, npoints)
    points = [(x_, y) for x_ in xs] # 2D points
    values = np.array([vec(point) for point in points])
    return values, xs

def get_distance_from_cutoff(cutoff: float, x: np.array, vec: np.array, center=0.0):
    ind = np.absolute( vec - cutoff).argmin()
    return abs(x[ind] - center), ind

class lmbda_param(UserExpression):
    def __init__(self, subdomains, lambda1,lambda2, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.lambda1 = lambda1
        self.lambda2 = lambda2        
    def eval_cell(self, values, x, cell):
        if self.subdomains[cell.index] == 2:
            values[0] = self.lambda2
        else:
            values[0] = self.lambda1

class mu_param(UserExpression):
    def __init__(self, subdomains, mu1, mu2, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.mu1 = mu1
        self.mu2 = mu2        
    def eval_cell(self, values, x, cell):
        if self.subdomains[cell.index] == 2:
            values[0] = self.mu2
        else:
            values[0] = self.mu1

class lmbda_param(UserExpression):
    def __init__(self, subdomains, lambda1,lambda2, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.lambda1 = lambda1
        self.lambda2 = lambda2        
    def eval_cell(self, values, x, cell):
        if self.subdomains[cell.index] == 2:
            values[0] = self.lambda2
        else:
            values[0] = self.lambda1

class Gc_param(UserExpression):
    def __init__(self, subdomains, Gc1, Gc2, **kwargs):
        super().__init__(**kwargs)
        self.subdomains = subdomains
        self.Gc1 = Gc1
        self.Gc2 = Gc2

    def eval_cell(self, values, x, cell):
        if self.subdomains[cell.index] == 2:
            values[0] = self.Gc2
        else:
            values[0] = self.Gc1