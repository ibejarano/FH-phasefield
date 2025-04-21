

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

# === Agrupar por E y graficar Curvatura vs h ===
por_E = defaultdict(list)
for row in data:
    por_E[row['E']].append((row['H'], row['Curvatura_a']))

for E_val, pares in por_E.items():
    pares.sort()
    h_vals, curv_vals = zip(*pares)
    plt.figure()
    plt.plot(h_vals, curv_vals, marker='o')
    plt.xlabel("H altura")
    plt.ylabel("Curvatura (coef a)")
    plt.title(f"Curvatura vs h para E = {E_val}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"curvatura_vs_h_E{E_val}.png"))
    plt.close()

# === Agrupar por h y graficar Curvatura vs E ===
por_h = defaultdict(list)
for row in data:
    por_h[row['H']].append((row['E'], row['Curvatura_a']))

for h_val, pares in por_h.items():
    pares.sort()
    E_vals, curv_vals = zip(*pares)
    plt.figure()
    plt.plot(E_vals, curv_vals, marker='o')
    plt.xlabel("E")
    plt.ylabel("Curvatura (coef a)")
    plt.title(f"Curvatura vs E para H = {h_val}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"curvatura_vs_E_h{h_val}.png"))
    plt.close()