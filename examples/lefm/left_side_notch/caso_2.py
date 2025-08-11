import numpy as np

from pre_built_cases.lefm import run_left_notch


r, KI_est, KII_est = run_left_notch(geo_name="left_notch.geo")



np.savetxt(f"caso_1.csv", np.array([r, KI_est, KII_est]), header="r,KI,KII", delimiter=",", comments='')