# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt # Opcional: para visualización
import argparse
import sys
import os
from analisis_curvatura import calcular_datos_intrinsecos
# Import functions and boundary condition type from the other script

def calcular_datos_intrinsecos(x, y):
    """
    Calcula puntos discretos para las ecuaciones intrínsecas de una curva:
    - Ecuación de Cesàro: (longitud de arco s, curvatura kappa)
    - Ecuación de Whewell: (longitud de arco s, ángulo tangencial phi)

    Args:
        x (list or np.ndarray): Coordenadas x de los puntos de la curva.
        y (list or np.ndarray): Coordenadas y de los puntos de la curva.

    Returns:
        tuple: Una tupla conteniendo:
            - s (np.ndarray): Longitudes de arco correspondientes a los puntos calculados.
            - kappa (np.ndarray): Curvaturas calculadas en los puntos.
            - phi (np.ndarray): Ángulos tangenciales (desenvueltos) calculados en los puntos.
            Retorna (None, None, None) si la entrada es inválida.

    Notas:
        - Este método usa np.gradient y es principalmente para comparación en los ejemplos.
        - El resultado omite los puntos extremos donde los cálculos
          son menos fiables con este método.
        - El ángulo phi se devuelve "desenvuelto" para evitar saltos de 2*pi.
    """
    # --- Validación de Entrada ---
    if len(x) != len(y):
        print("Error: Las listas x e y deben tener la misma longitud.", file=sys.stderr)
        return None, None, None
    if len(x) < 3:
        print("Error: Se necesitan al menos 3 puntos para los cálculos.", file=sys.stderr)
        return None, None, None

    # --- Conversión a NumPy Arrays ---
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    # --- 1. Calcular Longitud de Arco (s) ---
    dx = np.diff(x)
    dy = np.diff(y)
    segment_lengths = np.sqrt(dx**2 + dy**2)
    s_full = np.concatenate(([0], np.cumsum(segment_lengths)))

    # --- 2. Calcular Derivadas Aproximadas ---
    xp = np.gradient(x) # dx/dt (aprox)
    yp = np.gradient(y) # dy/dt (aprox)
    xpp = np.gradient(xp) # d^2x/dt^2 (aprox)
    ypp = np.gradient(yp) # d^2y/dt^2 (aprox)

    # --- 3. Calcular Curvatura (kappa) ---
    numerator_kappa = np.abs(xp * ypp - yp * xpp)
    denominator_kappa = (xp**2 + yp**2)**1.5
    epsilon = 1e-12
    kappa_full = np.divide(numerator_kappa, denominator_kappa + epsilon,
                           out=np.zeros_like(numerator_kappa),
                           where=denominator_kappa > epsilon)

    # --- 4. Calcular Ángulo Tangencial (phi) ---
    phi_raw = np.arctan2(yp, xp)
    phi_full = np.unwrap(phi_raw)

    # --- 5. Alinear los Datos y Devolver Resultados ---
    # Omitir los puntos extremos (primero y último) donde los gradientes son menos precisos.
    s_inner = s_full[1:-1]
    kappa_inner = kappa_full[1:-1]
    phi_inner = phi_full[1:-1]

    # Filtrar posibles valores no finitos
    valid_indices_kappa = np.isfinite(kappa_inner)
    valid_indices_phi = np.isfinite(phi_inner)
    valid_indices = valid_indices_kappa & valid_indices_phi

    s_inner = s_inner[valid_indices]
    kappa_inner = kappa_inner[valid_indices]
    phi_inner = phi_inner[valid_indices]

    return s_inner, kappa_inner, phi_inner

def plot_results(x, y, s, kappa, phi, title_prefix, s_comp=None, kappa_comp=None, phi_comp=None, comp_label=None):
    """Genera gráficos de la curva original, Cesàro y Whewell."""
    try:
        plt.figure(figsize=(18, 5)) # Figura más ancha para 3 subplots

        # Subplot 1: Curva Original
        plt.subplot(1, 3, 1)
        plt.plot(x, y, 'b.-', label='Curva Original')
        plt.title(f'Curva Original ({title_prefix})')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.grid(True)
        plt.legend()

        # Subplot 2: Cesàro (kappa vs s)
        plt.subplot(1, 3, 2)
        plt.plot(s, kappa, 'r.-', label='Calculado κ(s)')
        if s_comp is not None and kappa_comp is not None:
            plt.scatter(s_comp, kappa_comp, marker='.', label=comp_label if comp_label else 'Comparación κ(s)')
        # Add theoretical line for circle example
        if title_prefix.lower() == "círculo":
            radius = np.mean(np.sqrt(x**2 + y**2)) # Estimate radius
            if radius > 1e-6:
                plt.axhline(1/radius, color='g', linestyle='--', label=f'Teórico κ≈{1/radius:.3f}')

        plt.title(f'Cesàro (κ vs s) - {title_prefix}')
        plt.xlabel('Longitud de Arco (s)')
        plt.ylabel('Curvatura (κ)')
        plt.ylim(bottom=min(0, np.min(kappa) if len(kappa)>0 else 0)) # Allow negative curvature if present
        plt.grid(True)
        plt.legend()

        # Subplot 3: Whewell (phi vs s)
        plt.subplot(1, 3, 3)
        plt.plot(s, phi, 'm.-', label='Calculado φ(s)')
        if s_comp is not None and phi_comp is not None:
            plt.scatter(s_comp, phi_comp, marker='.', label=comp_label if comp_label else 'Comparación φ(s)')

        plt.title(f'Whewell (φ vs s) - {title_prefix}')
        plt.xlabel('Longitud de Arco (s)')
        plt.ylabel('Ángulo Tangencial φ (rad)')
        plt.grid(True)
        plt.legend()

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print("\nMatplotlib no está instalado o ocurrió un error al graficar.", file=sys.stderr)
        print(f"Error: {e}", file=sys.stderr)
        print("Puedes instalarlo con: pip install matplotlib")

def run_examples():
    """Genera y procesa los ejemplos de círculo y parábola."""
    print("--- Ejecutando Ejemplos Internos ---")

    # --- Ejemplo 1: Círculo ---
    print("\n--- Ejemplo: Círculo ---")
    radio = 5.0
    num_puntos_circulo = 100
    # endpoint=False for better closure with NATURAL BC in analisis_curvatura
    theta = np.linspace(0, 2 * np.pi, num_puntos_circulo, endpoint=False)
    x_circulo = radio * np.cos(theta)
    y_circulo = radio * np.sin(theta)

    # Método 1: np.gradient (calcular_datos_intrinsecos) - For comparison
    s_circ_grad, kappa_circ_grad, phi_circ_grad = calcular_datos_intrinsecos(x_circulo, y_circulo)

    if s_circ_grad is not None:
        print("Resultados (Círculo):")
        print(f"  Curvatura teórica esperada: {1/radio:.4f}")
        print(f"  Curvatura media (np.gradient): {np.mean(kappa_circ_grad):.4f}")
        print(f"  Rango Ángulo (np.gradient): [{phi_circ_grad[0]:.4f}, {phi_circ_grad[-1]:.4f}] rad")
        plot_results(x_circulo, y_circulo, s_circ_grad, kappa_circ_grad, phi_circ_grad, "círculo")


    # --- Ejemplo 2: Parábola ---
    print("\n--- Ejemplo: Parábola ---")
    num_puntos_parabola = 50 # Fewer points to see differences more clearly
    x_parabola = np.linspace(-3, 3, num_puntos_parabola)
    y_parabola = x_parabola**2

    # Método 1: np.gradient
    s_par_grad, kappa_par_grad, phi_par_grad = calcular_datos_intrinsecos(x_parabola, y_parabola)



    if s_par_grad is not None:
        print("Resultados (Parábola):")
        print(f"  Curvatura máx (np.gradient): {np.max(kappa_par_grad):.4f}")
        print(f"  Rango Ángulo (np.gradient): [{phi_par_grad[0]:.4f}, {phi_par_grad[-1]:.4f}] rad")
        plot_results(x_parabola, y_parabola, s_par_grad, kappa_par_grad, phi_par_grad, "Parabola")


def process_csv_file(filename, skiprows, delimiter):
    """Lee un archivo CSV, calcula ecuaciones intrínsecas y plotea."""
    print(f"--- Procesando archivo CSV: {filename} ---")
    try:
        print(f"Leyendo datos de: {filename}")
        # Load only first two columns (X, Y)
        data = np.loadtxt(filename, comments='#', delimiter=delimiter, skiprows=skiprows, usecols=(0,1))
        if data.ndim != 2 or data.shape[1] != 2:
            raise ValueError("El archivo debe contener exactamente dos columnas (X, Y) después de saltar filas.")
        x = data[:, 0]
        y = data[:, 1]
        N = len(x)
        print(f"Puntos leídos: {N}")
        if N < 3:
             print("Error: Se necesitan al menos 3 puntos para los cálculos.", file=sys.stderr)
             return

    except FileNotFoundError:
        print(f"Error: Archivo no encontrado '{filename}'", file=sys.stderr)
        return
    except ValueError as e:
         print(f"Error leyendo el archivo '{filename}'. ¿Es correcto el delimitador ('{delimiter}') y el número de filas a saltar ({skiprows})?", file=sys.stderr)
         print(f"Detalle del error: {e}", file=sys.stderr)
         return
    except Exception as e:
        print(f"Error inesperado leyendo el archivo '{filename}': {e}", file=sys.stderr)
        return

    # Calcular usando analisis_curvatura (método preferido para datos externos)
    print("Calculando ecuaciones intrínsecas usando analisis_curvatura...")
    # Usamos NATURAL como condición de contorno por defecto para curvas abiertas/cerradas de CSV
    curvatura, tu_seg1, tu_seg2, ds = calcular_curvatura(x, y, N, NATURAL)

    if curvatura is None:
        print("No se pudo calcular la curvatura. Abortando.", file=sys.stderr)
        return


def main():
    parser = argparse.ArgumentParser(
        description="Calcula y grafica las ecuaciones intrínsecas (Cesàro, Whewell) de una curva. "
                    "Puede usar ejemplos internos o leer datos de un archivo CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Ejemplos de uso:\n"
               "  python cesaro.py --examples          # Ejecuta los ejemplos internos (círculo, parábola)\n"
               "  python cesaro.py datos_curva.csv     # Procesa el archivo CSV con configuración por defecto\n"
               "  python cesaro.py datos.txt --delimiter=' ' --skiprows=2 # Procesa archivo con espacio como delimitador y salta 2 filas"
    )
    parser.add_argument(
        'input_file',
        nargs='?', # Makes the argument optional
        type=str,
        help="Archivo CSV de entrada con dos columnas (x y). Requerido si no se usa --examples."
    )
    parser.add_argument(
        '--examples',
        action='store_true',
        help="Ejecuta los ejemplos internos (círculo y parábola) con comparación, en lugar de leer un archivo."
    )
    parser.add_argument(
        "--skiprows",
        type=int,
        default=0, # Default to skipping 0 header rows for general CSV
        help="Número de filas a saltar al inicio del archivo CSV (encabezado). Predeterminado: 0"
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        default=",", # Default to comma delimiter
        help="Delimitador de columnas en el archivo CSV. Predeterminado: ','"
    )

    args = parser.parse_args()

    if args.examples:
        if args.input_file:
            print("Advertencia: Se especificó un archivo de entrada pero también --examples. Ejecutando solo los ejemplos.", file=sys.stderr)
        run_examples()
    elif args.input_file:
        # Check if the file exists before processing
        if not os.path.isfile(args.input_file):
             print(f"Error: Archivo de entrada no encontrado '{args.input_file}'", file=sys.stderr)
             sys.exit(1)
        process_csv_file(args.input_file, args.skiprows, args.delimiter)
    else:
        # Neither --examples nor input_file was provided
        parser.print_help()
        print("\nError: Debes especificar un archivo de entrada o usar la opción --examples.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
