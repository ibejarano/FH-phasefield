import numpy as np
import argparse
import sys
import os
import matplotlib.pyplot as plt # Import matplotlib


def calcular_datos_intrinsecos(x, y, filter_negative_x=True):
    

    if len(x) != len(y):
        print("Error: Las listas x e y deben tener la misma longitud.", file=sys.stderr)
        return
    if len(x) < 3:
        print("Error: Se necesitan al menos 3 puntos para los cálculos.", file=sys.stderr)
        return

    # --- Conversión a NumPy Arrays ---
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)

    # --- NUEVO: Filtrar puntos con x < 0 si se solicita ---
    if filter_negative_x:
        print("Filtrando puntos con x < 0...")
        positive_x_indices = np.where(x >= 0)[0]
        if len(positive_x_indices) < 3:
            print("Error: Después de filtrar x < 0, quedan menos de 3 puntos.", file=sys.stderr)
            return None, None, None
        x = x[positive_x_indices]
        y = y[positive_x_indices]
        print(f"Puntos restantes después del filtrado: {len(x)}")
    # ------------------------------------------------------

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


# --- Función Principal ---
def main():
    parser = argparse.ArgumentParser(
        description="Calcula las ecuaciones de Cesaro (s, curvatura) y Whewell (s, theta) "
                    "a partir de un archivo de datos XY."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Archivo de entrada con dos columnas (x y) separadas por espacio/tab/coma. Opcionalmente con encabezado."
    )
    parser.add_argument(
        "--plot",
        action="store_true", # Makes it a flag, args.plot will be True if present
        help="Muestra gráficos de las ecuaciones de Whewell y Cesaro además de guardar los archivos."
    )
    parser.add_argument(
        "--skiprows",
        type=int,
        default=1, # Default to skipping 1 header row
        help="Número de filas a saltar al inicio del archivo de entrada (encabezado). Predeterminado: 1"
    )
    parser.add_argument(
        "--delimiter",
        type=str,
        default=",", # Default to comma delimiter
        help="Delimitador de columnas en el archivo de entrada. Predeterminado: ','"
    )
    parser.add_argument(
        "--order",
        action="store_true",
        help="Ordenar los puntos de entrada por la coordenada X ascendente antes de calcular. ¡Esto cambia la curva analizada!"
    )
    # --- NUEVO ARGUMENTO PARA PROMEDIAR ---
    parser.add_argument(
        "--average-every",
        type=int,
        default=0, # Default 0 means no averaging
        metavar='N',
        help="Promediar cada N puntos después de leer (y opcionalmente ordenar). Si N > 1, se realiza el promedio. Los puntos restantes se descartan. Predeterminado: 0 (sin promedio)."
    )
    # ------------------------------------

    args = parser.parse_args()
    input_filename = args.input_file
    average_n = args.average_every # Tamaño del grupo para promediar

    # 1. Leer datos del archivo de entrada
    try:
        print(f"Leyendo datos de: {input_filename}")
        # Use specified skiprows and delimiter
        data = np.loadtxt(input_filename, comments='#', delimiter="\t", skiprows=args.skiprows, usecols=(0,1))
        if data.ndim != 2 or data.shape[1] != 2:
            raise ValueError("El archivo debe contener exactamente dos columnas (X, Y) después de saltar filas.")
        U = data[:, 0]
        V = data[:, 1]
        N = len(U)
        print(f"Puntos leídos: {N}")
        if N < 3:
             print("Error: Se necesitan al menos 3 puntos inicialmente.", file=sys.stderr)
             sys.exit(1)

    except FileNotFoundError:
        print(f"Error: Archivo no encontrado '{input_filename}'", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
         print(f"Error leyendo el archivo '{input_filename}'. ¿Es correcto el delimitador ('{args.delimiter}') y el número de filas a saltar ({args.skiprows})?", file=sys.stderr)
         print(f"Detalle del error: {e}", file=sys.stderr)
         sys.exit(1)
    except Exception as e:
        print(f"Error inesperado leyendo el archivo '{input_filename}': {e}", file=sys.stderr)
        sys.exit(1)

    # --- ORDENAR SI SE SOLICITÓ ---
    if args.order:
        print("Ordenando puntos por coordenada X ascendente...")
        # Obtener los índices que ordenarían U
        sorted_indices = np.argsort(U)
        # Reordenar U y V usando esos índices
        U = U[sorted_indices]
        V = V[sorted_indices]
        print("Puntos ordenados.")
        # Advertencia: Esto cambia fundamentalmente la secuencia de puntos y, por lo tanto, la curva analizada.
        # La curvatura y el ángulo se calcularán para esta nueva secuencia ordenada por X.
    # -----------------------------

    # --- PROMEDIAR SI SE SOLICITÓ ---
    if average_n > 1:
        print(f"Promediando cada {average_n} puntos...")
        num_groups = N // average_n
        if num_groups == 0:
            print(f"Advertencia: No hay suficientes puntos ({N}) para formar ni un grupo de {average_n}. No se realizará el promedio.", file=sys.stderr)
        else:
            num_points_to_average = num_groups * average_n
            remainder = N % average_n
            if remainder > 0:
                print(f"Advertencia: Se descartarán los últimos {remainder} puntos que no forman un grupo completo de {average_n}.", file=sys.stderr)

            # Tomar solo los puntos que forman grupos completos
            U_to_avg = U[:num_points_to_average]
            V_to_avg = V[:num_points_to_average]

            # Reorganizar en grupos y calcular la media
            U_avg = U_to_avg.reshape(num_groups, average_n).mean(axis=1)
            V_avg = V_to_avg.reshape(num_groups, average_n).mean(axis=1)

            # Actualizar U, V y N
            U = U_avg
            V = V_avg
            N = num_groups # El nuevo número de puntos es el número de grupos
            print(f"Número de puntos después de promediar: {N}")

            # Volver a verificar si hay suficientes puntos después de promediar
            if N < 3:
                print(f"Error: Después de promediar, quedan {N} puntos, que son menos de los 3 necesarios para calcular la curvatura.", file=sys.stderr)
                sys.exit(1)
    # --------------------------------

    # 2. Calcular curvatura, tangentes y longitudes de segmento
    print("Calculando curvatura...")
    # Pass N and cc_type correctly
    s_in, kappa_in, phi_in = calcular_datos_intrinsecos(U, V)


    # Verificar si el cálculo de curvatura falló (e.g., N < 3)
    if kappa_in is None:
        print("No se pudo calcular la curvatura. Abortando.")
        sys.exit(1)


    # 4. Preparar nombres de archivo de salida
    base_name = os.path.splitext(input_filename)[0]
    # Añadir sufijos si se ordenó o promedió
    order_suffix = "_ordered" if args.order else ""
    avg_suffix = f"_avg{average_n}" if average_n > 1 else ""
    output_cesaro_whewell = f"{base_name}{order_suffix}{avg_suffix}_Curvatura.csv"
    output_points = f"{base_name}{order_suffix}{avg_suffix}_Points.csv"

    # 5. Guardar resultados
    try:

        print(f"Guardando ecuación de Cesaro y Whewell en: {output_cesaro_whewell}")
        np.savetxt(output_cesaro_whewell, np.vstack((s_in, kappa_in, phi_in)).T, fmt='%.6e', delimiter='\t',
                   header='LongitudArco_S\tCurvatura\tAngulo_Theta_rad', comments='')
        
        print(f"Guardando puntos promediados en: {output_points}")
        np.savetxt(output_points, np.vstack((U, V)).T, fmt='%.6e', delimiter='\t',
                   header='Point_X\tPoint_Y', comments='')

    except Exception as e:
        print(f"Error guardando los archivos de salida: {e}", file=sys.stderr)
        # Continue to plotting if requested, even if saving fails
        # sys.exit(1) # Optionally exit here

    # 6. Plotear si se solicitó
    if args.plot:
        print("Generando gráficos...")
        try:
            plot_title_suffix = ""
            if args.order:
                plot_title_suffix += " (Ordenado por X)"
            if average_n > 1:
                plot_title_suffix += f" (Promediado cada {average_n})"

            # --- Gráfico de Whewell ---
            plt.figure(figsize=(10, 6))
            plt.plot(S, Theta, marker='.', linestyle='-', label='Whewell (Calculado)')
            plt.plot(s_in, phi_in, marker='.', linestyle='-', label='Whewell (Intrinseco)')
            plt.xlabel("Longitud de Arco (S)")
            plt.ylabel("Ángulo de Whewell (Theta) [rad]")
            plt.title(f"Ecuación de Whewell{plot_title_suffix}\nArchivo: {os.path.basename(input_filename)}")
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{base_name}{order_suffix}{avg_suffix}_Whewell.png")

            # --- Gráfico de Cesaro ---
            plt.figure(figsize=(10, 6))
            plt.plot(S, curvatura, marker='.', linestyle='-', label='Cesaro (Calculado)')
            plt.plot(s_in, kappa_in, marker='.', linestyle='-', label='Cesaro (Intrinseco)')
            plt.xlabel("Longitud de Arco (S)")
            plt.ylabel("Curvatura")
            plt.title(f"Ecuación de Cesaro{plot_title_suffix}\nArchivo: {os.path.basename(input_filename)}")
            plt.grid(True)
            plt.legend()
            # plt.ylim(0, None) # Comentado por si la curvatura puede ser negativa
            plt.tight_layout()
            plt.savefig(f"{base_name}{order_suffix}{avg_suffix}_Cesaro.png")

            # --- Gráfico Original XY (para verificar) ---
            # Este gráfico mostrará los puntos en el orden en que fueron procesados
            # (ordenados y/o promediados si se usaron las opciones)
            plt.figure(figsize=(8, 8))
            plt.plot(U, V, marker='.', linestyle='-', label=f'Datos XY Procesados{plot_title_suffix}')
            if N > 0: # Check if there are points left to plot
                plt.scatter(U[0], V[0], color='red', s=50, zorder=5, label='Punto Inicial (post-procesado)') # Mark start point
            plt.xlabel("Coordenada X")
            plt.ylabel("Coordenada Y")
            plt.title(f"Datos XY Procesados{plot_title_suffix}\nArchivo: {os.path.basename(input_filename)}")
            plt.grid(True)
            plt.axis('equal') # Ensure aspect ratio is equal
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{base_name}{order_suffix}{avg_suffix}_points.png")


            print("Mostrando gráficos... Cierra las ventanas para finalizar.")
            plt.show() # Display all created figures

        except Exception as e:
            print(f"Error generando los gráficos: {e}", file=sys.stderr)
            print("Asegúrate de tener matplotlib instalado (`pip install matplotlib`)", file=sys.stderr)
            # sys.exit(1) # Optionally exit if plotting fails

    print("Proceso completado.")

# --- Punto de Entrada ---
if __name__ == "__main__":
    main()
