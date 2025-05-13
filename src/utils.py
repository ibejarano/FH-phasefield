from dolfin import *
import numpy as np
from tqdm import tqdm
import json

E = 2e7
nu = 0.3
Gc = 2.7 # Energía por unidad de superficie de fractura
Q0 = 1e-4
lmbda = E*nu / ((1+nu)  * (1-2*nu))
mu = E / (2*(1+nu))

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

def get_values_overline_x(x1, x2, y=0.0, vec=None , npoints=100):
    tol = 0.0001 # avoid hitting points outside the domain
    xs = np.linspace(x1 + tol, x2 - tol, npoints)
    points = [(x_, y) for x_ in xs] # 2D points
    values = np.array([vec(point) for point in points])
    return values, xs

def get_distance_from_cutoff(cutoff: float, x: np.array, vec: np.array, center=0.0):
    ind = np.absolute( vec - cutoff).argmin()
    return abs(x[ind] - center), ind

def save_stress(caseDir, step=100):
    mesh = Mesh(f"{caseDir}/malla.xml.gz")
    V = FunctionSpace(mesh, 'CG', 1)
    W = VectorFunctionSpace(mesh, 'CG', 1)
    Vsig = TensorFunctionSpace(mesh, "DG", degree=0)
    phit = Function(V, name="phi")
    ut = Function(W, name="displacement")
    sigt = Function(Vsig, name="stress")
    u_ts = TimeSeries(f"{caseDir}/u_series")
    phi_ts = TimeSeries(f"{caseDir}/phi_series")

    # Definimos las vars para almacenar las tensiones ppales
    Vp = VectorFunctionSpace(mesh, "DG",degree=1)
    svec11 = Function(W, name="S11")
    theta11 = Function(W, name="dir11")
    svec22 = Function(W, name="S22")
    theta22 = Function(W, name="dir22")

    xdmffile = XDMFFile(f"{caseDir}/output_stress.xdmf")
    xdmffile.parameters["flush_output"] = True
    xdmffile.parameters["functions_share_mesh"] = True

    times = u_ts.vector_times()

    pbar = tqdm(total=len(times[::step]))

    for t in times[::step]:
        pbar.update(1)
        u_ts.retrieve(ut.vector(), t)
        phi_ts.retrieve(phit.vector(), t)
        sigt.assign(project( (1-phit)**2 * sigma(ut), Vsig))
        sx = sigt.sub(0)
        sy = sigt.sub(3)
        txy = sigt.sub(1)
        svec11.assign(project(as_vector((sigma11(sx,sy,txy),0.)), Vp))
        svec22.assign(project(as_vector((sigma22(sx,sy,txy),0.)), Vp))

        THETAP =  thethap(sx,sy,txy)

        theta11.assign(project(as_vector((cos(THETAP), sin(THETAP))), Vp))
        theta22.assign(project(as_vector((sin(THETAP), -cos(THETAP))), Vp))

        xdmffile.write(sigt, t)
        xdmffile.write(ut, t)
        xdmffile.write(svec11, t)
        xdmffile.write(svec22, t)
        xdmffile.write(theta11, t)
        xdmffile.write(theta22, t)

    pbar.close()
    return

def sigma11(sigma_x, sigma_y, tau_xy): # magnitude of first principal stress
    return ((sigma_x+sigma_y)/2 + sqrt(((sigma_x-sigma_y)/2)**2 + tau_xy**2))

def sigma22(sigma_x, sigma_y, tau_xy): # magnitude of second principal stress
    return ((sigma_x+sigma_y)/2 - sqrt(((sigma_x-sigma_y)/2)**2 + tau_xy**2))

def thethap(sigma_x, sigma_y, tau_xy):
    return atan((2 * tau_xy) / (sigma_x - sigma_y))


def read_data(fname):
    with open(f"data/{fname}.json", "r") as f:
        data = json.load(f)
    
    E = data["E"]
    nu = data["nu"]
    lmbda = E*nu / ((1+nu)  * (1-2*nu))
    mu = E / (2*(1+nu))
    data.update({"mu": mu})
    data.update({"lmbda": lmbda})
    return data

# ("stress_0"+"stress_4")/2 + sqrt((("stress_0"-"stress_4")/2)^2 + "stress_1"^2)

class CrackDomainAngle(SubDomain):
	def inside(self, x, on_boundary):
		angulo = 38
		alpha = -angulo*(3.1415/180)
		senalpha = np.sin(alpha)
		cosalpha = np.cos(alpha)
		xp = (x[0]*cosalpha - x[1]*senalpha, x[0]*senalpha + x[1]*cosalpha)
		return abs(xp[0]) <= 0.5 and abs(xp[1]) <= 0.015

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


class SimulationData:
    def __init__(self, E=None, nu=None, Gc=None, Q0=None, h_elem=None, l_h=None):
        self.E = E
        self.nu = nu
        self.Gc = Gc
        self.Q0 = Q0
        self.h_elem = h_elem
        self.l_h = l_h

    def set_tolerances(self, tol_vol=0.001, tol_phi=1e-3):
        self.tol_vol = tol_vol
        self.tol_phi = tol_phi
    
    def set_simulation_parameters(self, dt=None, t_final=None):
        self.dt = dt
        self.t_final = t_final

    def set_initial_conditions(self, pinit=None, l0=None, w0=None):
        self.pinit = pinit
        self.l0 = l0
        self.w0 = w0

if __name__ == "__main__":
    #save_stress("results/curving_H_20", step=2)
    import json

    with open("results/test/simulation_output.txt", "r") as f:
        for l in f:
            if l[0] == "{":
                a = json.loads(l)
                print(a)
                print(a["Eyoung"])
    print("Finished")