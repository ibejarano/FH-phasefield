import numpy as np
from sys import argv
from funcs import lame_lambda, lame_mu
# Parámetros
meshName = "lab"
caseDir = "S8"
saveSteps = 10
aspect_hl = 4 # Rel aspecto long. transicion vs altura minima de elemento
l = aspect_hl*7.5e-4 # Longitud de transición

# Material
E = 2450
nu = 0.3
Gc = 0.1 # Energía por unidad de superficie de fractura
g = 9.81
densidad = 1200 * 0.035 # kg/m3 * h = km/m2

# Fzas exteriores
px = 10.
py = 0.

# Condiciones iniciales
w_init = 5e-4
#l_init = 5e-3
l_init = 0.01
p_init = 10
center = [0.0, 0.0]

# Tolerancias
t = 0
Tfinal = 400
#u_r = 0.007
delta_p = 0.1
deltaT = 0.1
tolphi = 1e-3
toldispl = 1e-6
step = 0

## aumentar presión si converge rápido
MAXITE = 5
deltap_counter = 0
DELTAP_AUG = 5

## MESHING ##

mesh = Mesh(meshName+ ".xml")
subdomains = MeshFunction('size_t',mesh,meshName+"_physical_region.xml")
boundary_markers = MeshFunction('size_t',mesh, meshName+"_facet_region.xml")

# 1 . MESH TIENE QUE SER DOMAIN
# --- mesh esta importado, tiene que ser domain
V = FunctionSpace(mesh, 'CG', 1)
W = VectorFunctionSpace(mesh, 'CG', 1)

qo = 3.4 * (1e-6/60)  # 3.4 ml/min

WW = FunctionSpace(mesh, 'DG', 0)
p, q = TrialFunction(V), TestFunction(V)
u, v = TrialFunction(W), TestFunction(W)


lmbda = lame_lambda(E, nu)
mu = lame_mu(E, nu)

class CrackDomainAngle(SubDomain):
  def inside(self, x, on_boundary):
    xc = x.copy()
    xc[0] = x[0] - center[0]
    xc[1] = x[1] - center[1]
    angulo = 360-30
    alpha = -angulo*(3.1415/180)
    senalpha = np.sin(alpha)
    cosalpha = np.cos(alpha)
    xp = (xc[0]*cosalpha - xc[1]*senalpha, xc[0]*senalpha + xc[1]*cosalpha)
    return abs(xp[0]) <= l_init and abs(xp[1]) <= w_init

class CrackDomain(SubDomain):
  def inside(self, x, on_boundary):
    return abs(x[0] - center[0]) <= l_init and abs(x[1] - center[1]) <= w_init

bcright = DirichletBC(W.sub(0), Constant(0.0), boundary_markers, 10)
bciny = DirichletBC(W.sub(0), Constant(0.0), boundary_markers, 40)

# Aplicar desplazamiento en la dirección y positiva
#load = Expression("t", t=0.0, degree=1)
bcbottom = DirichletBC(W.sub(1), Constant(0.0), boundary_markers, 20)

bc_u = [bcright, bcbottom]#, bciny]

# Condicion de borde de la fractura, se marca la entalla inicial con el valor de 1
crack = CrackDomain()
bc_phi = [DirichletBC(V, Constant(1.0), crack)]

unew, uold, ut = Function(W), Function(W), Function(W, name="displacement")
pnew, pold, Hold, phit = Function(V), Function(V), Function(V), Function(V, name="phi")

# Fzas externas
p_value = p_init
pressure = Constant(p_value)
ds = Measure("ds")(subdomain_data=boundary_markers)
n = FacetNormal(mesh)

# Funcional del desplazamiento eq (x)
txx = Constant((-px, 0.0))
tyy = Constant((0.0, py))
# Funcional del desplazamiento eq (x)
E_du = (1-pold)**2*inner(epsilon(v), sigma_0(u))*dx + pressure * inner(v, grad(pold))*dx + inner(txx, v) * ds(30) + inner(tyy, v) * ds(50)

# Funcional phi eq (x)
E_phi = (Gc*l*inner(grad(p), grad(q)) + ((Gc/l)+2.0 *H(uold, unew, Hold))*inner(p,q)-2.0*H(uold, unew,Hold)*q)*dx


p_disp = LinearVariationalProblem(lhs(E_du), rhs(E_du), unew, bc_u)
p_phi = LinearVariationalProblem(lhs(E_phi), rhs(E_phi), pnew, bc_phi)

solver_disp = LinearVariationalSolver(p_disp)
solver_phi = LinearVariationalSolver(p_phi)

#solver_disp.parameters["linear_solver"] = "gmres"
#solver_phi.parameters["linear_solver"] = "gmres"

#solver_disp.parameters["preconditioner"] = "ilu"
#solver_phi.parameters["preconditioner"] = "ilu"

solver_disp.solve()


os.mkdir(caseDir)

file_results = XDMFFile(f"{caseDir}/{caseDir}.xdmf")
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True


fname = open(f"{caseDir}/output.csv", 'w')


while t<=Tfinal:
  step += 1

  deltap_counter += 1 # Si converge sin necesidad de bajar la Presion entonces aceleramos

  if deltap_counter == DELTAP_AUG:
    deltap_counter = 0
    if delta_p < 0.1:
      delta_p = 0.1
      print("Increasing Delta P")
  ## Solve
  t += deltaT
  p_value += delta_p*p_value
  pressure.assign(p_value)

  print(f"Step: {step} Solving for t: {t:.9f} pressure: {p_value:.2f} -----")

  ite = 0
  err_u = 1
  err_phi = 1

  while (err_u > toldispl) or (err_phi > tolphi):
    ite += 1

    if ite % MAXITE == 0:
      p_value *= 0.98
      delta_p = 0.01
      deltap_counter = 0

      pressure.assign(p_value)
      print(f"Disminuir presión a P: {p_value:.2f}")
      # Resetear uold y pold
      uold.assign(ut)
      pold.assign(phit)
      Hold.assign(project(psi(uold), WW))

    solver_disp.solve()
    solver_phi.solve()
    err_u = errornorm(unew, uold, norm_type='l2', mesh=None)
    err_phi = errornorm(pnew, pold, norm_type='l2', mesh=None)
    uold.assign(unew)
    pold.assign(pnew)
    Hold.assign(project(H(uold, unew, Hold), WW))

    print(f"Ite: {ite} -- error phi: {err_phi:.2e} -- error u: {err_u:.2e}" )

  ## Converge
  ut.assign(unew)
  phit.assign(pnew)

  if step%saveSteps == 0:
    file_results.write(ut  , t)
    file_results.write(phit, t)
    print("XMDF Updated")

  fname.write(str(step) + ",")
  fname.write(str(t) + ",")
  fname.write(str(p_value) + "\n")
  #print("Reaccion", reaccion)
  print("Step Converged.")



fname.close()
print("Simulation completed")
