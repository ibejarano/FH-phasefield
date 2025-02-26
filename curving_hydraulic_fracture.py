from dolfin import *
import numpy as np
from sys import argv
import os
from utils import createSave
from tqdm import tqdm
#import matplotlib.pyplot as plt
# Quitar mensajes de compilacion
set_log_active(False)
#set_log_level(LogLevel.INFO)

# Props Material y flujo
E = 2e7
nu = 0.3
Gc = 2.7 # Energía por unidad de superficie de fractura
Q0 = 1e-4

# Mallado
h_elem = 2e-3 # TODO: Cambiar con el tamaño de la malla en zona de fractura
aspect_hl = 3 # aspect_hl = e = l/h
l = aspect_hl*h_elem # Longitud de transición

# Condiciones Iniciales
l0 = 0.05/2
w0 = h_elem

# Control de simulacion
DT = 0.01
TOL_PHI = 1e-3 # Tolerancia de phi
TOL_VOL = 0.001 # 0.1% de tolerancia de volumen inyectado
T_FINAL = 0.63
p_init = 10

# directorios
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

def epsilon(u):
  return sym(grad(u))

def sigma(u):
  return 2.0*mu*epsilon(u)+lmbda*tr(epsilon(u))*Identity(len(u))

# Densidad de energia elastica eq (12)
def psi(u):
  return 0.5*(lmbda+mu)*(0.5*(tr(epsilon(u)) + abs(tr(epsilon(u)))))**2 + mu*inner(dev(epsilon(u)), dev(epsilon(u)))

# Parámetro de trayectoria eq (11)
def H(uold, unew, Hold):
  return conditional(lt(psi(uold), psi(unew)), psi(unew), Hold)

class CrackDomainAngle(SubDomain):
	def inside(self, x, on_boundary):
		angulo = 38
		alpha = -angulo*(3.1415/180)
		senalpha = np.sin(alpha)
		cosalpha = np.cos(alpha)
		xp = (x[0]*cosalpha - x[1]*senalpha, x[0]*senalpha + x[1]*cosalpha)
		return abs(xp[0]) <= 0.5 and abs(xp[1]) <= 0.015

class CrackDomain(SubDomain):
	def inside(self, x, on_boundary):
		center = [0, 0.0]
		return abs(x[0] - center[0]) <= l0 and abs(x[1] - center[1]) <= w0


bcright = DirichletBC(W, (0.0, 0.0), boundary_markers, 10)
bcleft  = DirichletBC(W, (0.0, 0.0), boundary_markers, 30)
#bcup = DirichletBC(W.sub(1), 0.0, boundary_markers, 40)
#bcbottom  = DirichletBC(W.sub(1), 0.0, boundary_markers, 20)
bc_u = [bcleft, bcright]

# Condicion de borde de la fractura, se marca la entalla inicial con el valor de 1
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
outfile.write("Mod young E: " + str(E) + "\n")
outfile.write("nu: " + str(nu) + "\n")
outfile.write("Flow rate: " + str(Q0) + "\n")
outfile.write("Gc: " + str(Gc) + "\n")
outfile.write("helem: " + str(h_elem) + "\n")
outfile.write(" -- Condiciones iniciales -- \n")
outfile.write("Longitud inicial: " + str(l0) + "\n")
outfile.write("Opening inicial: " + str(w0) + "\n")
outfile.write(" -- Parametros y tolerancias -- \n")
outfile.write("delta T: " + str(DT) + "\n")
outfile.write("Tol volumen: " + str(TOL_VOL) + "\n")
outfile.write("Tol phi: " + str(TOL_PHI) + "\n")
outfile.write("Tiempo final: " + str(T_FINAL) + "\n")
outfile.write("Pr inicial: " + str(p_init) + "\n")
outfile.close()

step = 0
solver_phi.solve()
pold.assign(pnew)

pbar = tqdm(total=int(T_FINAL/DT))

while t <= T_FINAL:
	pbar.update(1)
	step += 1
	print(f"Step: {step}")
	ite = 0
	V0 = assemble( inner(grad(phit), -ut) * dx )
	errDV = 1
	errDV1 = None
	errDV2 = None
	DV0 = DT * Q0
	err_phi = 1
	while err_phi > TOL_PHI:
		while abs(errDV)/DV0 > TOL_VOL:
			ite += 1
			DV0 = DT * Q0
			try:
				pn = pn1 - errDV1 * (pn1 - pn2)/(errDV1 - errDV2)
			except:
				pn *= 1.001

			pressure.assign(pn)
			solver_disp.solve()

			VK = assemble( inner(grad(pold), -unew) * dx )
			errDV = DV0 - (VK - V0) # Delta Vol - ( Delta vol )
			uold.assign(unew)
			Hold.assign(project(psi(unew), WW))
			try:
				pn2 = pn1
				errDV2 = errDV1
				pn1 = pn
				errDV1 = errDV
			except:
				pn1 = pn
				errDV1 = errDV
		
		solver_phi.solve()
		err_phi = errornorm(pnew, pold, norm_type='l2', mesh=mesh)
		pold.assign(pnew)
		errDV = 1
		if ite > 20:
			print("Simulation finished by iterations")
			break

	if ite > 20:
		print(" too much iterations ")
		break

	ut.assign(unew)
	phit.assign(pnew)

	vol_frac = assemble( inner(grad(phit), -ut) * dx )
	t += DT
	fname.write(str(t) + ",")
	fname.write(str(pn) + ",")
	fname.write(str(vol_frac) + "\n")

	print(f"Converge t: {t:.4f} dt: {DT:.2e} --- Iteraciones: {ite}")
	# Save files
	if step % 2 == 0:
		out_xml.write(ut, t)
		out_xml.write(phit, t)
		print("Saving VTK")
	u_ts.store(ut.vector(), t)
	phi_ts.store(phit.vector(), t)

fname.close()
