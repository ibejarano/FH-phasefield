from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

L, H = 100, 50
a_fr = 0.1
h_entalla = 0.001
H_sup = 0.2

mesh = Mesh("mesh.xml")
boundaries = MeshFunction("size_t", mesh, "mesh_facet_region.xml")

V = VectorFunctionSpace(mesh, "P", 1)

tol = 1E-14

def bottom(x, on_boundary):
    return on_boundary and near(x[1], -H, tol) and np.less(np.abs(x[0]), L/50)

bc_bottom = DirichletBC(V, Constant((0, 0)), bottom)
# bc_bottom = DirichletBC(V, Constant((0, 0)), boundaries, 1)

def left(x, on_boundary):
    return on_boundary and near(x[0], 0, tol)

# bc_left = DirichletBC(V, Constant((0, 0)), left)
# bc_left = DirichletBC(V, Constant((0, 0)), boundaries, 4)

bcs = [bc_bottom]

E = 2e8
nu = 0.3

mu = E / (2 * (1 + nu))

lmbda = E*nu / ((1 + nu)*(1 - 2*nu))

def epsilon(u):
    return sym(grad(u))

def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(2) + 2*mu*epsilon(u)

p0 = -1e5
p0 = 0.0
p1 =- 1e5

traction = Constant((p0, 0.0))

n = FacetNormal(mesh)
# internal_pressure = Constant((0.0, p1))

class RightBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L/2, tol)

#right_boundary = RightBoundary()

# boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
# boundaries.set_all(0)
# right_boundary.mark(boundaries, 1)

#ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
ds = ds(subdomain_data=boundaries)
u = TrialFunction(V)
v = TestFunction(V)

a = inner(sigma(u), epsilon(v))*dx

p_values = np.geomspace(1e5, 1e5, 1)
u_sup_values = np.zeros(len(p_values))
u_plus_values = np.zeros(len(p_values))
u_minus_values = np.zeros(len(p_values))
sxx_vals = np.zeros(len(p_values))
syy_vals = np.zeros(len(p_values)) 
sxy_vals = np.zeros(len(p_values))
theta_vals = np.zeros(len(p_values))

fig, axs = plt.subplots(3, 1, sharex=True)

for i, p1 in enumerate(p_values):
    print("Internal pressure", p1)
    internal_pressure = Constant((0.0, -p1))
    L_form = dot(traction, v)*ds(2) - dot(traction, v)*ds(4) 
    L_form  += - dot(internal_pressure, v)*ds(5) + dot(internal_pressure, v)*ds(6)
    u_sol = Function(V)
    solve(a == L_form, u_sol, bcs)

    point_disp_1 = u_sol(0.0, H_sup)
    u_plus = u_sol(0.0, h_entalla)
    u_minus = u_sol(0.0, -h_entalla)


    S = sigma(u_sol)

    W = TensorFunctionSpace(mesh, 'P', 1)

    S_proj = project(S, W)

    S_computed = S_proj(a_fr*1.01/2, 0.0)
    sxx = S_computed[0]
    syy = S_computed[3]
    sxy = S_computed[1]

    u_sup_values[i] = point_disp_1[1]
    u_plus_values[i] = u_plus[1]
    u_minus_values[i] = u_minus[1]

    theta_ppal = np.atan(2*sxy/(sxx - syy)) / 2

    print("Theta", np.rad2deg(theta_ppal), "Grados")
    theta_vals[i] = np.rad2deg(theta_ppal)

    num_points = 100
    x_vals = np.linspace(a_fr*1.00/2, a_fr*1.5/2, num_points)
    y_fixed = h_entalla
    sigma_xx_vals = np.zeros(num_points)
    sigma_yy_vals = np.zeros(num_points)
    sigma_xy_vals = np.zeros(num_points)
    PI_vals = np.zeros(num_points)

    for i, x in enumerate(x_vals):
        point = np.array([x, y_fixed])
        S_val = S_proj(point)
        
        sigma_xx_vals[i] = S_val[0]
        sigma_yy_vals[i] = S_val[3]
        sigma_xy_vals[i] = S_val[1]
        PI_vals[i] = 2* S_val[1] / (S_val[0] - S_val[3])

    #axs[0].plot(x_vals, sigma_xx_vals, label=f'{p1} Pa  ')
    axs[1].plot(x_vals, sigma_yy_vals, label=f'{p1} Pa  ')
    axs[2].plot(x_vals, sigma_xy_vals, label=f'{p1} Pa  ')
    axs[0].plot(x_vals, np.rad2deg(np.atan(PI_vals))*0.5)
    #axs[2].set_xlabel('x (m)')
    axs[0].set_ylabel('Tensión XX (Pa)')
    axs[1].set_ylabel('Tensión YY (Pa)')
    axs[2].set_ylabel('Tensión XY (Pa)')
    plt.title(f'Tensiones a lo largo de y = {y_fixed}')


#np.savetxt(f"out_afr_{a_fr*100:.0f}cm.csv", np.array([a_fr*np.ones(len(p_values)), p_values, u_sup_values, u_plus_values, u_minus_values, theta_vals]).T, header="a_fr, pr, u_sup, u_plus, u_minus, theta", delimiter=",",  comments='')
plt.legend()

if False:
    plt.figure()
    plt.plot(p_values, u_plus_values - u_minus_values, label="Opening")
    #plt.plot(-p_values, wdiff_values, label="Opening diff")
    #plt.yscale("log")
    plt.xscale("log")
    plt.legend()


    plt.figure()

    plt.plot(p_values, u_sup_values, label="Top displacement")
    #plt.plot(p_values, u_plus_values, label=r"$u^+$")
    #plt.plot(p_values, np.abs(u_minus_values), label=r"Abs $u^-$")
    plt.yscale("log")
    plt.xscale("log")
    plt.legend()

    plt.figure()


# np.savetxt("a3_results.csv", np.array([p_values, u_sup_values, u_plus_values, u_minus_values]).T, delimiter=",", header="p,u_sup,u_plus,u_minus")


if False:
    plt.figure()
    num_points = 100
    x_vals = np.linspace(a_fr*1.01/2, 1, num_points)
    y_fixed = 0.0
    sigma_xx_vals = np.zeros(num_points)
    sigma_yy_vals = np.zeros(num_points)
    sigma_xy_vals = np.zeros(num_points)

    for i, x in enumerate(x_vals):
        point = np.array([x, y_fixed])
        S_val = S_proj(point)
        
        sigma_xx_vals[i] = S_val[0]
        sigma_yy_vals[i] = S_val[3]
        sigma_xy_vals[i] = S_val[1]



    plt.plot(x_vals, sigma_xx_vals, label=r'$\sigma_{xx}$')
    plt.plot(x_vals, sigma_yy_vals, label=r'$\sigma_{yy}$')
    plt.plot(x_vals, sigma_xy_vals, label=r'$\tau_{xy}$')
    plt.xlabel('x (m)')
    plt.ylabel('Tensión (Pa)')
    plt.title(f'Tensiones a lo largo de y = {y_fixed}')
    plt.legend()

plt.show()

if False:

    import matplotlib.pyplot as plt

    plt.figure()

    plot(mesh, linewidth=0.2)
    c = plot(u_sol*10, mode="displacement", title="Desplazamientos")
    plt.colorbar(c)
    plt.xlim(-a_fr/1.8, a_fr/1.8)
    plt.ylim(-0.02, 0.02)




    plt.figure()
    num_points = 100
    y_vals = np.linspace(h_entalla, H_sup, num_points)
    x_fixed = 0.0
    disply_vals = np.zeros(num_points*2)
    for i, y in enumerate(y_vals):
        point = np.array([x_fixed, y])
        #S_val = S_proj(point)
        disply_vals[i] = u_sol(x_fixed, y)[1]
        #sigma_xx_vals[i] = S_val[0, 0]
        #sigma_yy_vals[i] = S_val[1, 1]
        #sigma_xy_vals[i] = S_val[0, 1]



    plt.plot(y_vals, disply_vals[:num_points], label=r'$u_{y}$')
    plt.xlabel('x (m)')
    plt.ylabel('Tensión (Pa)')
    plt.title(f'Tensiones a lo largo de y = {x_fixed}')
    plt.legend()

    y_vals = np.linspace(-H_sup, -h_entalla, num_points)
    for i, y in enumerate(y_vals):
        point = np.array([x_fixed, y])
        disply_vals[i+num_points] = abs(u_sol(x_fixed, y)[1])

    plt.plot(y_vals, disply_vals[num_points:], label=r'$u_{y}$')
    plt.grid(True)

    #  Componente sigma_xx
    S = sigma(u_sol)

    W = TensorFunctionSpace(mesh, 'P', 1)

    S_proj = project(S, W)
    plt.figure()
    c = plot(S_proj[0, 0], title="Sigma_xx (Tensión normal horizontal)")
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xlim(a_fr/3, a_fr/1.5)
    plt.ylim(-0.02, 0.02)



    # Componente sigma_yy
    plt.figure()
    c = plot(S_proj[1, 1], title="Sigma_yy (Tensión normal vertical)")
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xlim(a_fr/3, a_fr/1.5)
    plt.ylim(-0.02, 0.02)

    # Componente sigma_xy
    plt.figure()
    c = plot(S_proj[0, 1], title="Sigma_xy (Tensión cortante)")
    plt.colorbar(c)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xlim(a_fr/3, a_fr/1.5)
    plt.ylim(-0.02, 0.02)

    plt.show()