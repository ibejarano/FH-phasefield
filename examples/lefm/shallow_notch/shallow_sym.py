import numpy as np
from pre_built_cases.lefm import run_shallow_case_symm


# Vamos a editar el geo antes de simular
H = 5
E = 164.3e9
nu = 0.32
esp = 0.02
p1 = 183e6*esp # Presion interna
beta = 0
pxx = 0

ratios_casos = [0]#, 0.2, 0.5, 1, 2, 5]
# casos_H = np.geomspace(0.1, 10, 30)
a_cracks = np.geomspace(H/100, H*5, 40)
a_crack = H / 100

#[KI, KII], [us, up, um] = run_shallow_case_symm(H, p1, pxx, a_crack, E=E, nu=nu, beta=beta, geo_name="sym_shallow.geo", save_vtu=True, mesh=True)


for i, ratio in enumerate(ratios_casos):
    K_calcs = np.zeros((len(a_cracks), 2))
    up_arr = np.zeros(len(a_cracks))
    um_arr = np.zeros(len(a_cracks))
    us_arr = np.zeros(len(a_cracks))
    pxx = ratio*p1
    for j, a_crack in enumerate(a_cracks):
        [KI, KII], [us, up, um] = run_shallow_case_symm(H, p1, pxx, a_crack, E=E, nu=nu, beta=beta, geo_name="sym_shallow.geo", save_vtu=True, mesh=True)
        K_calcs[j] = [KI, KII]
        up_arr[j] = up
        um_arr[j] = um
        us_arr[j] = us


    aH = a_cracks / H
    np.savetxt(f"./sym_beta_{beta}/ratio_{ratio}.csv", np.array([aH, K_calcs[:, 0], K_calcs[:, 1], us_arr, up_arr,um_arr]).T, header="aH,KI,KII,us,up,um", delimiter=",", comments='')
