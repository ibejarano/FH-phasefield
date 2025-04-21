import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.interpolate import make_interp_spline
from collections import defaultdict
import sys

output_dir = "./fracturepaths"

if len(sys.argv) < 2:
    print("Uso: python post_process_curves.py <nombre_propiedad>")
    sys.exit(1)

propiedad_clave = sys.argv[1]

# Directorio con los resultados

# Buscar todos los archivos CSV que no son de propiedades
files = [f for f in os.listdir(output_dir) if f.endswith(".csv") and not f.startswith("propiedades_")]

# Agrupar archivos por valor de la propiedad deseada
grupos = {}
if False:

    for file in files:
        case_name = file.replace(".csv", "")
        props_path = os.path.join(output_dir, f"propiedades_{case_name}.csv")

        if os.path.exists(props_path):
            with open(props_path, newline='') as f:
                reader = csv.DictReader(f)
                props = next(reader)
                valor_clave = props.get(propiedad_clave, None)

                if valor_clave is not None:
                    if valor_clave not in grupos:
                        grupos[valor_clave] = []
                    grupos[valor_clave].append(file)

    # Graficar por grupo
    for valor, archivos in grupos.items():
        plt.figure(figsize=(10, 6))
        for file in archivos:
            path = os.path.join(output_dir, file)
            x = []
            y = []
            with open(path, newline='') as csvfile:
                reader = csv.reader(csvfile)
                next(reader)
                raw_points = [(float(row[0]), float(row[1])) for row in reader]

            x_groups = defaultdict(list)
            for xi, yi in raw_points:
                x_groups[xi].append(yi)

            x_unique = []
            y_avg = []
            for xi in sorted(x_groups.keys()):
                x_unique.append(xi)
                y_avg.append(np.mean(x_groups[xi]))

            x_sorted = np.array(x_unique)
            y_sorted = np.array(y_avg)

            if len(x_sorted) > 3:
                xnew = np.linspace(x_sorted.min(), x_sorted.max(), 300)
                coef = np.polyfit(x_sorted, y_sorted, deg=2)
                poly = np.poly1d(coef)
                ynew = poly(xnew)
                plt.plot(xnew, ynew, label=file.replace(".csv", ""))
            else:
                plt.plot(x_sorted, y_sorted, label=file.replace(".csv", ""))

        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title(f"Curvas agrupadas por {propiedad_clave} = {valor}")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"curvas_{propiedad_clave}_{valor}.png")

        plt.close()

# === Análisis de sensibilidad: exportar curvatura y propiedades combinadas ===
resultado_global = []
for file in files:
    case_name = file.replace(".csv", "")
    props_path = os.path.join(output_dir, f"propiedades_{case_name}.csv")
    csv_path = os.path.join(output_dir, file)

    if not os.path.exists(props_path):
        continue

    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        raw_points = [(float(row[0]), float(row[1])) for row in reader]

    x_groups = defaultdict(list)
    for xi, yi in raw_points:
        x_groups[xi].append(yi)

    x_unique = []
    y_avg = []
    for xi in sorted(x_groups.keys()):
        x_unique.append(xi)
        y_avg.append(np.mean(x_groups[xi]))

    x_sorted = np.array(x_unique)
    y_sorted = np.array(y_avg)

    if len(x_sorted) > 3:
        coef = np.polyfit(x_sorted, y_sorted, deg=2)
        a = coef[0]
    else:
        a = None

    with open(props_path, newline='') as f:
        reader = csv.DictReader(f)
        props = next(reader)
        props["Caso"] = case_name
        props["Curvatura_a"] = a
        resultado_global.append(props)

# Guardar resultados
if resultado_global:
    keys = resultado_global[0].keys()
    with open(os.path.join(output_dir, "curvatura_vs_propiedades.csv"), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(resultado_global)

print("Archivo curvatura_vs_propiedades.csv generado.")

# === Gráficos de dispersión: propiedades vs curvatura ===
import matplotlib.pyplot as plt

propiedades_numericas = [k for k in resultado_global[0].keys() if k not in ("Caso", "Curvatura_a")]

for prop in propiedades_numericas:
    x_vals = []
    y_vals = []
    for fila in resultado_global:
        try:
            x_vals.append(float(fila[prop]))
            y_vals.append(float(fila["Curvatura_a"]))
        except (ValueError, TypeError):
            continue

    if x_vals and y_vals:
        plt.figure()
        plt.scatter(x_vals, y_vals)
        plt.xlabel(prop)
        plt.ylabel("Curvatura (coef a)")
        plt.title(f"Curvatura vs {prop}")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"curvatura_vs_{prop}.png"))
        plt.close()
