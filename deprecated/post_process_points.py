import os
import pandas as pd
import numpy as np
import re
import json
import argparse # Para argumentos de línea de comando
import matplotlib.pyplot as plt # Para plotear
import sys # Para salir del script

# --- Configuración ---
RESULTS_BASE_DIR = "results"
PARAMS_FILENAME = "lab.json"
X_COL_NAME = 'Points_0' # Columna X en el CSV de puntos
Y_COL_NAME = 'Points_1' # Columna Y en el CSV de puntos
POLYNOMIAL_DEGREE = 2
OUTPUT_FILENAME = "simulations_summary_with_coeffs.csv"

# --- Función para parsear el archivo de parámetros (sin cambios) ---
def parse_params_file(filepath):
    params = None
    try:
        with open(filepath, 'r') as f:
            all_lines = f.readlines()
            if len(all_lines) > 2:
                json_string_content = "".join(all_lines[2:])
                params = json.loads(json_string_content)
                if not isinstance(params, dict):
                    print(f"    Error: El contenido JSON en {filepath} no es un diccionario. Tipo: {type(params)}")
                    return None
            else:
                print(f"    Advertencia: {filepath} tiene <= 2 líneas. No hay JSON.")
                return None
    except FileNotFoundError:
        print(f"    Error: Archivo de parámetros no encontrado: {filepath}")
        return None
    except json.JSONDecodeError as e:
        print(f"    Error: Fallo al decodificar JSON en {filepath}. Error: {e}")
        return None
    except Exception as e:
        print(f"    Error leyendo o parseando {filepath}: {e}")
        return None
    return params

# --- Función para leer CSV y ajustar polinomio (Refactorizada) ---
def fit_polynomial_for_simulation(points_csv_path, x_col, y_col, degree):
    """
    Lee un archivo CSV de puntos, realiza un ajuste polinómico y devuelve los datos y coeficientes.

    Args:
        points_csv_path (str): Ruta al archivo CSV.
        x_col (str): Nombre de la columna X.
        y_col (str): Nombre de la columna Y.
        degree (int): Grado del polinomio.

    Returns:
        tuple: (x_data, y_data, coeffs) si tiene éxito, (None, None, None) si falla.
               x_data, y_data son Series de Pandas. coeffs es un array de numpy.
    """
    try:
        points_df = pd.read_csv(points_csv_path)

        if x_col not in points_df.columns or y_col not in points_df.columns:
            print(f"  Error: Columnas '{x_col}' o '{y_col}' no encontradas en {points_csv_path}.")
            return None, None, None

        x_data = points_df[x_col].dropna()
        y_data = points_df[y_col].dropna()

        if len(x_data) != len(y_data):
            print(f"  Advertencia: Discrepancia en longitud X ({len(x_data)}) e Y ({len(y_data)}) en {points_csv_path}.")
            return None, None, None

        if len(x_data) <= degree:
            print(f"  Advertencia: No hay suficientes puntos ({len(x_data)}) para ajustar grado {degree} en {points_csv_path}.")
            return None, None, None

        coeffs = np.polyfit(x_data, y_data, degree)
        print(f"  Ajuste polinómico realizado. Coeficientes: {coeffs}")
        return x_data, y_data, coeffs

    except FileNotFoundError:
        print(f"  Error: Archivo CSV no encontrado: {points_csv_path}")
        return None, None, None
    except pd.errors.EmptyDataError:
        print(f"  Error: Archivo CSV está vacío: {points_csv_path}")
        return None, None, None
    except Exception as e:
        print(f"  Error procesando CSV o ajustando polinomio para {points_csv_path}: {e}")
        return None, None, None

# --- Función para plotear el ajuste ---
def plot_fit(x_data, y_data, coeffs, sim_name):
    """Genera y muestra un gráfico del ajuste polinómico."""
    plt.figure(figsize=(10, 6))
    font_settings = {
        'font.size': 24,          # Tamaño base
        'axes.titlesize': 18,     # Tamaño del título del gráfico
        'axes.labelsize': 16,     # Tamaño de las etiquetas de los ejes (X, Y)
        'xtick.labelsize': 14,    # Tamaño de las etiquetas de las marcas en el eje X
        'ytick.labelsize': 14,    # Tamaño de las etiquetas de las marcas en el eje Y
        'legend.fontsize': 14,    # Tamaño del texto de la leyenda
        'legend.title_fontsize': 15 # Tamaño del título de la leyenda (si aplica)
    }
    plt.rcParams.update(font_settings)
    
    plt.scatter(x_data, y_data, label='Datos Originales', s=10, alpha=0.7) # s=tamaño, alpha=transparencia

    # Crear puntos para la curva ajustada
    x_fit = np.linspace(x_data.min(), x_data.max(), 200) # Más puntos para curva suave
    y_fit = np.polyval(coeffs, x_fit) # Evaluar el polinomio en los puntos x_fit

    plt.plot(x_fit, y_fit, color='red', label=f'Ajuste Polinomico (Grado {len(coeffs)-1})')

    plt.title(f'Ajuste Polinómico para Simulación: {sim_name}')
    plt.xlabel("Coordenada X [mts]")
    plt.ylabel("Coordenada Y [mts]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout() # Ajusta el layout para que no se solapen etiquetas
    print(f"\nMostrando gráfico para {sim_name}... Cierra la ventana del gráfico para continuar.")
    plt.show()

# --- Procesamiento principal ---
def main():
    # Configurar el parser de argumentos
    parser = argparse.ArgumentParser(description="Procesa resultados de simulación, ajusta polinomios y opcionalmente plotea un caso específico.")
    parser.add_argument("--plot", metavar="SIM_NAME", type=str,
                        help="Nombre de la simulación específica a plotear. Si se usa, solo se procesa y plotea esta simulación.")

    args = parser.parse_args()

    # --- MODO PLOTEO ---
    if args.plot:
        sim_to_plot = args.plot
        print(f"--- Modo Ploteo Activado para: {sim_to_plot} ---")
        sim_folder_path = os.path.join(RESULTS_BASE_DIR, sim_to_plot)

        if not os.path.isdir(sim_folder_path):
            print(f"Error: La carpeta de simulación '{sim_to_plot}' no existe en '{RESULTS_BASE_DIR}'.")
            sys.exit(1) # Salir con código de error

        points_csv_path = os.path.join(sim_folder_path, f"{sim_to_plot}.csv")

        # Realizar el ajuste para esta simulación
        x_data, y_data, coeffs = fit_polynomial_for_simulation(points_csv_path, X_COL_NAME, Y_COL_NAME, POLYNOMIAL_DEGREE)

        if coeffs is not None:
            # Si el ajuste fue exitoso, plotear
            plot_fit(x_data, y_data, coeffs, sim_to_plot)
        else:
            print(f"No se pudo generar el gráfico para '{sim_to_plot}' debido a errores previos.")
            sys.exit(1) # Salir con código de error

        print("--- Modo Ploteo Finalizado ---")
        sys.exit(0) # Salir exitosamente tras plotear

    # --- MODO NORMAL (Generar CSV Resumen) ---
    print(f"--- Modo Normal: Generando Resumen '{OUTPUT_FILENAME}' ---")
    all_results_data = []
    print(f"Iniciando procesamiento en la carpeta: '{RESULTS_BASE_DIR}'")

    for sim_folder_name in os.listdir(RESULTS_BASE_DIR):
        sim_folder_path = os.path.join(RESULTS_BASE_DIR, sim_folder_name)

        if os.path.isdir(sim_folder_path):
            print(f"\nProcesando simulación: {sim_folder_name}")

            params_file_path = os.path.join(sim_folder_path, PARAMS_FILENAME)
            points_csv_path = os.path.join(sim_folder_path, f"{sim_folder_name}.csv")

            if not os.path.exists(params_file_path):
                print(f"  Advertencia: Archivo '{PARAMS_FILENAME}' no encontrado. Saltando.")
                continue
            if not os.path.exists(points_csv_path):
                print(f"  Advertencia: Archivo '{sim_folder_name}.csv' no encontrado. Saltando.")
                continue

            simulation_params = parse_params_file(params_file_path)
            if simulation_params is None:
                print(f"  Error: No se pudieron leer parámetros. Saltando.")
                continue

            # Ajustar polinomio usando la función refactorizada
            x_data, y_data, coeffs = fit_polynomial_for_simulation(points_csv_path, X_COL_NAME, Y_COL_NAME, POLYNOMIAL_DEGREE)

            polynomial_coeffs_dict = {}
            if coeffs is not None:
                polynomial_coeffs_dict = {f'coeff_{chr(ord("a")+i)}': coeff for i, coeff in enumerate(coeffs)}
            else:
                # Añadir NaNs si el ajuste falló pero queremos incluir los parámetros
                for i in range(POLYNOMIAL_DEGREE + 1):
                    polynomial_coeffs_dict[f'coeff_{chr(ord("a")+i)}'] = np.nan

            # Combinar parámetros y coeficientes
            current_sim_data = simulation_params.copy()
            current_sim_data['simulation_id'] = sim_folder_name
            current_sim_data.update(polynomial_coeffs_dict)

            all_results_data.append(current_sim_data)

        else:
            print(f"Omitiendo '{sim_folder_name}', no es una carpeta.")

    # Crear y guardar el DataFrame final
    if all_results_data:
        final_df = pd.DataFrame(all_results_data)
        cols = final_df.columns.tolist()
        param_cols = [c for c in cols if not c.startswith('coeff_') and c != 'simulation_id']
        coeff_cols = sorted([c for c in cols if c.startswith('coeff_')])
        # Asegurarse de que 'simulation_id' esté presente
        id_col = ['simulation_id'] if 'simulation_id' in cols else []
        final_df = final_df[id_col + param_cols + coeff_cols]

        final_df['K1c'] = np.sqrt(final_df['Gc']*final_df['E']/(1 - final_df['nu']**2))
        final_df['chi'] = final_df['px'] * np.sqrt(final_df['H']) / final_df['K1c'] 

        try:
            final_df.to_csv(OUTPUT_FILENAME, index=False)
            print(f"\n¡Proceso completado! Resultados guardados en: '{OUTPUT_FILENAME}'")
            print("\nVista previa del DataFrame:")
            print(final_df.head())
        except Exception as e:
            print(f"\nError al guardar el archivo CSV '{OUTPUT_FILENAME}': {e}")
    else:
        print("\nNo se procesaron datos de ninguna simulación.")

    print("--- Modo Normal Finalizado ---")

# Asegurarse de que el código principal se ejecute solo cuando se corre el script directamente
if __name__ == "__main__":
    main()
