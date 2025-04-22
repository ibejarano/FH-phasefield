

import csv
import os
import matplotlib.pyplot as plt
from collections import defaultdict

input_csv = "./fracturepaths/curvatura_vs_propiedades.csv"
output_dir = "./fracturepaths"

# Leer los datos del CSV
data = []
with open(input_csv, newline='') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            row['E'] = float(row['E'])
            row['H'] = float(row['H'])
            row['Curvatura_a'] = float(row['Curvatura_a'])
            data.append(row)
        except ValueError:
            continue

def export_plot(var1, var2):
    """
    var1: str - hue variable
    var2: str - xlabel variable
    """
    print("Plotting: ", var2, "hue:", var1)
    por_var1 = defaultdict(list)
    for row in data:
        por_var1[row[var1]].append((row[var2], row['Curvatura_a']))

    plt.figure()
    for var1_val, pares in por_var1.items():
        pares.sort()
        var2_vals, curv_vals = zip(*pares)
        plt.plot(var2_vals, curv_vals, marker='o', label=f"{var1} = {var1_val}")

    plt.xlabel(var2)
    plt.ylabel("Curvatura (coef a)")
    plt.title(f"Curvatura vs {var2}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"curvatura_vs_{var2}_hue{var1}.png"))
    plt.close()


#export_plot("E", "Gc")
#export_plot("Gc", "E")
export_plot("E", "Qo")
export_plot("Qo", "E")