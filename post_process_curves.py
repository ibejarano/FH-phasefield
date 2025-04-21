import os
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy.interpolate import make_interp_spline
from collections import defaultdict
import sys
if len(sys.argv) < 2:
    print("Uso: python post_process_curves.py <nombre_propiedad>")
    sys.exit(1)

propiedad_clave = sys.argv[1]

# Directorio con los resultados
output_dir = "./fracturepaths"

# Buscar todos los archivos CSV que no son de propiedades
files = [f for f in os.listdir(output_dir) if f.endswith(".csv") and not f.startswith("propiedades_")]

# Agrupar archivos por valor de la propiedad deseada
grupos = {}

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
