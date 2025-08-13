import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
plt.scatter(sim['aH'], theta, label="Original")

theta_rad = np.deg2rad(theta)

alpha = 100

# Opcion 1
KI_mod = KI * (np.cos(theta_rad)**2 + alpha * np.sin(theta_rad)**2)
KII_mod =  KII * ( (1-alpha) * np.cos(theta_rad)*np.sin(theta_rad) )
theta_mod = theta_expr(KI_mod, KII_mod)
plt.plot(sim['aH'], theta_mod, label="Modificada Opc. 1")
ys1 = compute_y(xs1, theta_mod)

# Opcion 1B
KI_mod =  (KI * np.cos(theta_rad)**2 + alpha * np.sin(theta_rad)**2)
KII_mod =  ( (KI-alpha) * np.cos(theta_rad)*np.sin(theta_rad) )
theta_mod = theta_expr(KI_mod, KII_mod)
plt.plot(sim['aH'], theta_mod, label="Modificada Opc. 1 B")
ys1b = compute_y(xs1, theta_mod)

# Opcion 2
KI_mod = ( KI + alpha * np.sin(theta_rad)**2)
KII_mod = (KII - alpha * np.cos(theta_rad)*np.sin(theta_rad) )
theta_mod = theta_expr(KI_mod, KII_mod)
ys2 = compute_y(xs1, theta_mod)

plt.plot(sim['aH'], theta_mod, label="Modificada Opc. 2")

# Oocion 3 (incremental)
KI_3 = np.zeros(len(KI))
KII_3 = np.zeros(len(KII))
theta_3 = np.zeros(len(KI))
KI_3[0] = KI[0]
KII_3[0] = KII[0]
theta_i = np.deg2rad(theta_expr(KI, KII)[0])
theta_3[0] = theta_i

for i in range(len(KI_3)):
  KI_3[i] = KI[i] + alpha * np.sin(theta_i)**2
  KII_3[i] = KII[i] - alpha * np.sin(theta_i)* np.cos(theta_i)
  theta_i = theta_expr(KI_3[i], KII_3[i])
  theta_3[i] = theta_i
  theta_i = np.deg2rad(theta_i)


plt.plot(sim['aH'], theta_3, label="Modificada Opc. 3")
ys3 = compute_y(xs1, theta_3)

plt.xscale("log")
plt.legend()
plt.show()

plt.scatter(xs1, ys0, label="Original")
plt.plot(xs1, ys1, label="Mod Opc. 1")
plt.plot(xs1, ys1b, label="Mod Opc. 1B")
plt.plot(xs1, ys2, label="Mod Opc. 2")
plt.plot(xs1, ys3, label="Mod Opc. 3")

plt.legend()
plt.show()