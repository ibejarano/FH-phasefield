import numpy as np
from dolfin import plot

from pre_built_cases.lefm import run_left_notch

Lcrack = 0.025
E = 164.3e9
nu = 0.32
esp = 0.02
p1 = 183e6*esp # Presion interna

r, KI_est, KII_est = run_left_notch(p1, Lcrack, E, nu, geo_name="left_notch.geo", save_vtu=True)



np.savetxt(f"caso_2.csv", np.array([r, KI_est, KII_est]).T, header="r,KI,KII", delimiter=",", comments='')