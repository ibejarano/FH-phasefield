from dolfin import *
import json
from sys import argv
import os
from utils import createSave, read_data
from mesh_setup import setup_mesh_and_spaces
from material_model import epsilon, select_sigma, select_psi, H

# Quitar mensajes de compilacion
set_log_active(False)
#set_log_level(LogLevel.INFO)
data = read_data("lab")
# Props Material y flujo
E = data["E"]
nu = data["nu"]
Gc = data["Gc"] # Energía por unidad de superficie de fractura
Q0 = data["Qo"]
pxx = data["px"]
# Mallado
h_elem = data["h"] # TODO: Cambiar con el tamaño de la malla en zona de fractura
aspect_hl = data["aspect_hl"] # aspect_hl = e = l/h
l = aspect_hl*h_elem # Longitud de transición

# Condiciones Iniciales
l0 = data["linit"]
w0 = h_elem
p_init = 100

# Control de simulacion
TOL_PHI = 1e-3 # Tolerancia de phi
TOL_VOL = 0.001 # 0.1% de tolerancia de volumen inyectado
DT = 0.0002
T_FINAL = DT * 10000

assert len(argv) == 3 , "Case name not found and mesh"
caseDir = os.path.join("./results", argv[1])
meshName = caseDir+"/"+argv[2]
## MESHING ##
mesh, subdomains, boundary_markers, V, W, WW = setup_mesh_and_spaces(meshName)
 
# 1 . MESH TIENE QUE SER DOMAIN
# --- mesh esta importado, tiene que ser domain
p, q = TrialFunction(V), TestFunction(V)
u, v = TrialFunction(W), TestFunction(W)

# Parametros de Lame (material isotropo)
lmbda = E*nu / ((1+nu)  * (1-2*nu))
mu = E / (2*(1+nu))
data.update({"mu": mu})
data.update({"lmbda": lmbda})
psi_model = data.get("psi_model", "linear")
psi = select_psi(psi_model)
sigma = select_sigma("linear")  # o "hyperelastic"
print(f"Usando modelo de energía: {psi_model}")

#bcright = DirichletBC(W, (0.0, 0.0), boundary_markers, 10)
#bcleft  = DirichletBC(W, (0.0, 0.0), boundary_markers, 30)
#bcup = DirichletBC(W.sub(1), 0.0, boundary_markers, 40)
bcbottom  = DirichletBC(W, (0.0, 0.0), boundary_markers, 20)
bc_u = [bcbottom]

# Condicion de borde de la fractura, se marca la entalla inicial con el valor de 1
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
ds = Measure("ds", subdomain_data=boundary_markers)

px_vec = Constant((pxx, 0.0))

E_du = (1-pold)**2*inner(epsilon(v), sigma(u, mu, lmbda))*dx + pressure * inner(v, grad(pold))*dx + dot(px_vec, v)*ds(10) - dot(px_vec, v)*ds(30) 

# Funcional de la variable phi eq(14)
E_phi = (Gc*l*inner(grad(p), grad(q)) + ((Gc/l)+2.0 *H(unew, Hold, data, psi))*inner(p,q)-2.0*H(unew, Hold, data, psi)*q)*dx


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
outfile.write(json.dumps(data))
outfile.close()

step = 0
solver_phi.solve()
pold.assign(pnew)

while t <= T_FINAL:
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
			Hold.assign(project(psi(unew, mu, lmbda), WW))
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
