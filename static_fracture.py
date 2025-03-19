from dolfin import *
import numpy as np
import json
from sys import argv
import os
from utils import createSave, line_integral, sigma, epsilon, psi, H, compute_opening_cutoff, compute_opening_grad
from tqdm import tqdm
import matplotlib.pyplot as plt
# Quitar mensajes de compilacion
set_log_active(False)
#set_log_level(LogLevel.INFO)

# Props Material y flujo
E = 2e8
nu = 0.3
Gc = 2.7 # Energía por unidad de superficie de fractura
Q0 = 1e-3

# Mallado
h_elem = 1e-3 # TODO: Cambiar con el tamaño de la malla en zona de fractura
aspect_hl = 3 # aspect_hl = e = l/h
l = aspect_hl*h_elem # Longitud de transición

# Condiciones Iniciales
l0 = 0.025
w0 = h_elem

# Control de simulacion
TOL_PHI = 1e-3 # Tolerancia de phi
T_FINAL = 800
TOL_VOL = 0.001 # 0.1% de tolerancia de volumen inyectado
DT = 0.001
p_init = 100


## MESHING ##
assert len(argv) == 3 , "Case name not found and mesh"
caseDir = os.path.join("./results", argv[1])
meshName = caseDir+"/"+argv[2]
## MESHING ##
mesh = Mesh(meshName+ ".xml")
subdomains = MeshFunction('size_t',mesh,meshName+"_physical_region.xml")
boundary_markers = MeshFunction('size_t',mesh, meshName+"_facet_region.xml")

# 1 . MESH TIENE QUE SER DOMAIN
# --- mesh esta importado, tiene que ser domain
V = FunctionSpace(mesh, 'CG', 1)
W = VectorFunctionSpace(mesh, 'CG', 1)

WW = FunctionSpace(mesh, 'DG', 0)
p, q = TrialFunction(V), TestFunction(V)
u, v = TrialFunction(W), TestFunction(W)


# Parametros de Lame (material isotropo)
lmbda = E*nu / ((1+nu)  * (1-2*nu))
mu = E / (2*(1+nu))

bcright = DirichletBC(W, (0.0, 0.0), boundary_markers, 10)
bcleft  = DirichletBC(W, (0.0, 0.0), boundary_markers, 30)
bc_u = [bcleft, bcright]

class CrackDomain(SubDomain):
	def inside(self, x, on_boundary):
		center = [0, 0.0]
		return abs(x[0] - center[0]) <= l0 and abs(x[1] - center[1]) <= w0
	
crack = CrackDomain()

bc_phi = [DirichletBC(V, Constant(1.0), crack)]

unew, uold, ut = Function(W), Function(W), Function(W, name="displacement")
pnew, pold, Hold, phit = Function(V), Function(V), Function(V), Function(V, name="phi")

# Funcional del desplazamiento eq (13)
pressure = Constant(p_init)

E_du = (1-pold)**2*inner(epsilon(v), sigma(u))*dx + pressure * inner(v, grad(pold))*dx

# Funcional de la variable phi eq(14)
E_phi = (Gc*l*inner(grad(p), grad(q)) + ((Gc/l)+2.0 *H(uold, unew, Hold))*inner(p,q)-2.0*H(uold, unew,Hold)*q)*dx


p_disp = LinearVariationalProblem(lhs(E_du), rhs(E_du), unew, bc_u)
p_phi = LinearVariationalProblem(lhs(E_phi), rhs(E_phi), pnew, bc_phi)

solver_disp = LinearVariationalSolver(p_disp)
solver_phi = LinearVariationalSolver(p_phi)

#solver_disp.parameters["linear_solver"] = "gmres"
solver_phi.parameters["linear_solver"] = "gmres"

#solver_disp.parameters["preconditioner"] = "jacobi"
solver_phi.parameters["preconditioner"] = "ilu"

solver_disp.solve()
solver_phi.solve()

out_xml, u_ts, phi_ts = createSave(mesh, caseDir, "xml")
fname = open(f"./{caseDir}/output.csv", 'w')

t = 0
pn = p_init
pressure.assign(pn)

outfile = open(f"./{caseDir}/simulation_output.txt", 'w')
outfile.write(" -- Algoritmo --- \n")
outfile.write(" Algoritmo con control de volumen \n")
#outfile.write(json.dumps(caseData))
outfile.close()

step = 0
solver_phi.solve()
pold.assign(pnew)

for pn in [2000, 10000]:
	pressure.assign(pn)

	err_phi = 1
	while err_phi > TOL_PHI:
		solver_disp.solve()
		uold.assign(unew)
		Hold.assign(project(psi(unew), WW))
		solver_phi.solve()
		err_phi = errornorm(pnew, pold, norm_type='l2', mesh=mesh)
		pold.assign(pnew)

	ut.assign(unew)
	phit.assign(pnew)

	l0c = l0*1.07
	lf = l0*2.5

	gradphit = Function(W)
	gradphit.assign(project( grad(phit) , W))

	xsgrad_t, wxgrad_t = compute_opening_grad(ut, phit, gradphit, h_elem/20, lf)

	xs_t , wx_t = compute_opening_cutoff(ut, phit, h_elem/2, lf)
	Eprime = E/(1-nu**2)
	wx_an = 2*pn/(Eprime) *l0c  * np.sqrt(1 - (xsgrad_t/l0c)**2)
	out_xml.write(ut, pn)
	out_xml.write(phit, pn)
	#plt.plot(xs_t, wx_t/10, label=f"pressure = {pn:.2f}")
	plt.plot(xsgrad_t, wxgrad_t/10, label=f"$P_f$ = {pn:.2f} Pa")
	plt.plot(xsgrad_t, wx_an, "k--", label=f"Analytic sol.")

plt.ylabel("Opening (m)")
plt.xlabel("Fracture long. [m]")
plt.legend()
plt.show()

print("Simulation completed")