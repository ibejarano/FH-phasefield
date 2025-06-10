import numpy as np
import matplotlib.pyplot as plt
import json

def read_data(fname):
    with open(f"{fname}.json", "r") as f:
        data = json.load(f)
    
    E = data["E"]
    nu = data["nu"]
    lmbda = E*nu / ((1+nu)  * (1-2*nu))
    mu = E / (2*(1+nu))
    data.update({"mu": mu})
    data.update({"lmbda": lmbda})
    return data

caseDirs = [f"20250427_q{str(i)}" for i in [1,2,3,4]]
caseDirsh2 = [f"20250428_h15_q{str(i)}" for i in [1,2,3] ]
for c in caseDirsh2:
    caseDirs.append(c)

caseDirs = [caseDirs[i] for i in (1, 5)]

if False:

    for caseDir in caseDirs:
        data = np.loadtxt(f"./results/{caseDir}/output.csv", delimiter=",")
        ts, ps = data[:, 0], data[:, 1]
        plt.plot(ts, ps, label=f" Q = {caseDir[-1]}")

    plt.yscale("log")
    plt.xscale("log")
    plt.legend()
    plt.show()
    #plt.savefig("pvsq.png")

if True:
    avgpoints = 8
    plt.figure(figsize=(10, 6))

    for caseDir in caseDirs:
        props = read_data(f"./results/{caseDir}/lab")

        data = np.loadtxt(f"./results/{caseDir}/{caseDir}_ordered_avg{avgpoints}_Points.csv", delimiter="\t", skiprows=1)
        longs, ws = data[:, 0], data[:, 1]
        
        plt.scatter(longs, ws, marker='.',  label=f"H = {props["H"]} cm")
        plt.xlabel("X coord.")
        plt.ylabel("Y coord.")
        plt.title(f"Trayectoria fractura - Q={props["Qo"]:.1e}")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()


    plt.figure(figsize=(10, 6))

    for caseDir in caseDirs:
        props = read_data(f"./results/{caseDir}/lab")

        data = np.loadtxt(f"./results/{caseDir}/{caseDir}_ordered_avg{avgpoints}_Whewell.csv", delimiter="\t", skiprows=1)
        longs, ws = data[:, 0], data[:, 1]
        
        plt.plot(longs, ws, marker='.', linestyle='-', label=f"H = {props["H"]} cm")
        plt.xlabel("Longitud de Arco (S)")
        plt.ylabel("Ángulo de Whewell (Theta) [rad]")
        plt.title(f"Ecuación de Whewell - Q={props["Qo"]:.1e}")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

    plt.figure(figsize=(10, 6))
    for caseDir in caseDirs:
        props = read_data(f"./results/{caseDir}/lab")

        data = np.loadtxt(f"./results/{caseDir}/{caseDir}_ordered_avg{avgpoints}_Cesaro.csv", delimiter="\t", skiprows=1)
        longs, cs = data[:, 0], data[:, 1]
        plt.plot(longs, cs, marker='.', linestyle='-',  label=f"H = {props["H"]} cm")
        plt.xlabel("Longitud de Arco (S)")
        plt.ylabel("Curvatura")
        plt.title(f"Ecuación de Cesaro - Q={props["Qo"]:.1e}")
        plt.grid(True)
        plt.legend()
        # plt.ylim(0, None) # Comentado por si la curvatura puede ser negativa
        plt.tight_layout()
    plt.show()