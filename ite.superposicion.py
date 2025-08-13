import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import root

def KI_expr(x):
  # x = a/H
  return 1+(1/8)*x + (2/5)*x**2

def KII_expr(x):
  # x = a/H
  return (1/6)*x**2

def compute_y(x, alpha_deg):
  alpha_rad = np.deg2rad(alpha_deg)
  #df = pd.read_csv('curva.csv')

  # 2. Convertir alpha de grados a radianes
  #df['alpha_rad'] = np.deg2rad(df['alpha'])

  # 3. Calcular y usando integración trapezoidal acumulada
  #x = df['x'].values
  #tan_alpha = np.tan(df['alpha_rad'].values)
  tan_alpha = np.tan(alpha_rad)

  # Inicializamos y en cero
  y = np.zeros_like(x)

  # Acumulación paso a paso
  for i in range(1, len(x)):
      dx = x[i] - x[i-1]
      # Trapecio entre tan_alpha[i-1] y tan_alpha[i]
      y[i] = y[i-1] + 0.5 * (tan_alpha[i-1] + tan_alpha[i]) * dx

  return y

def theta_expr(KI, KII):

  R = np.where(np.abs(KI) > np.finfo(float).eps, KII / KI, 0.0)

  # Soluciones para t = tan(θ/2)
  sqrt_term = np.sqrt(1 + 12*R**2)
  t_plus  = (1 + sqrt_term) / (6*R)
  t_minus = (1 - sqrt_term) / (6*R)

  # Ángulos candidatos
  theta_plus  =  2*np.arctan(t_plus)
  theta_minus =  2*np.arctan(t_minus)

  # Elegir la raíz que cumple θ₀ * R < 0
  theta0 = np.where(theta_plus * R < 0, theta_plus, theta_minus)

  # Casos de modo I puro (KII ≈ 0): θ₀ = 0
  theta0 = np.where(np.abs(R) < np.finfo(float).eps, 0.0, theta0)

  return -np.rad2deg(theta0)

def get_fractura_coords(df: pd.DataFrame, H=1):
  xs = df['aH']*H
  KI = df['KI']
  KII = df['KII']
  theta = theta_expr(KI, KII)
  ys = compute_y(xs, theta)

  ys = ys[ ys < H ]
  xs = xs[ : len(ys) ]

  return xs, ys

def compute_from_ajuste(H):
  xs = np.linspace(1e-5, 5, 100)
  KI = KI_expr(xs, beta=0)
  KII = KII_expr(xs, beta=0)
  theta = theta_expr(KI, KII)
  ys = compute_y(xs, theta)
  ys = ys[ ys < H ]
  xs = xs[ : len(ys) ]

  return xs, ys

sim = pd.read_csv('/home/ignacio/repos/FH-phasefield/examples/lefm/shallow_notch/beta_0/ratio_0.csv')
H = 1
cte = np.sqrt(np.pi * H * sim['aH'])
KI =  sim['KI']
KII =  sim['KII']
xs1, ys0 = get_fractura_coords(sim)
theta = theta_expr(KI, KII)
# plt.scatter(sim['aH'], theta, label="Modelo 1")
plt.scatter(xs1, ys0, label="Modelo 1")

theta_rad = np.deg2rad(theta)

alpha = 0

# Oocion 3 (iterativo)
n = len(KI)
KI_3 = np.zeros(n)
KII_3 = np.zeros(n)
theta_3 = np.zeros(n)

def f1(K01, theta, ratio):
    # EJEMPLO: reemplaza por tu f1 real
    return K01 + ratio * np.sin(theta)**2

def f2(K02, theta, ratio):
    # EJEMPLO
    return K02 - ratio * np.cos(theta)*np.sin(theta)

def f3(K1, K2):
    # EJEMPLO
    alpha = theta_expr(K1, K2)
    return np.deg2rad(alpha)

def solve_fixed_point(K01, K02, theta0, ratio=0, alpha=0.01, tol=0.00001, maxit=5000):
    theta = theta0.copy()
    for k in range(maxit):
        K1 = f1(K01, theta, ratio)
        K2 = f2(K02, theta, ratio)
        theta_new = (1-alpha)*theta + alpha*f3(K1, K2)

        if np.linalg.norm(theta_new - theta, ord=np.inf) < tol:
            return K1, K2, np.rad2deg(theta_new), k+1, True
        theta = theta_new
    return K1, K2, theta, maxit, False



# Uso:
for ratio in [0, 0.5, 1, 2, 5, 10, 20]:
  K1, K2, theta_res, iters, ok = solve_fixed_point(KI, KII, theta_3, ratio=ratio)
  print("Convergió:", ok, "en iters:", iters)
  label = r" $\alpha$ =" + str(ratio)
  if ok:
    #plt.plot(sim['aH'], theta_res, label=label)
    ys3 = compute_y(xs1, theta_res)

    # plt.xscale("log")
    plt.xlabel("Coordenada X (mts)")
    plt.ylabel("Coordenada Y (mts)")
    #plt.ylabel(r"Angulo $\theta$ (grad)")
    #plt.show()

    plt.plot(xs1, ys3, label=label)

    #plt.legend()
  else:
    print("No convirgio", f"ratio {ratio}")

plt.legend()
plt.show()
