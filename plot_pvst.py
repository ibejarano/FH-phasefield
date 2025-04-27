import numpy as np
import matplotlib.pyplot as plt


caseDir = "h10_2"
data = np.loadtxt(f"./results/{caseDir}/output.csv", delimiter=",")

ts, ps = data[:, 0], data[:, 1]

plt.plot(ts, ps)
plt.show()