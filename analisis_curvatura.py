import numpy as np
import argparse
import sys
import os
import matplotlib.pyplot as plt # Import matplotlib

# --- Constantes (equivalentes a #define) ---
SI = 1
NO = 0

PERIODICA = 0
NATURAL = 1
LIBRE = 2

# --- Configuración (ajustar según sea necesario) ---
ES_SOLIDO_DE_REVOLUCION = NO
CCONTORNO = NATURAL # Elegir entre PERIODICA, NATURAL, LIBRE
LONG_C = 1000.0 # Solo relevante si CCONTORNO == PERIODICA

def calcular_curvatura(U, V, N, cc_type, long_c=None):
    """
    Calcula la longitud de los segmentos, vectores tangentes,
    y la curvatura en cada punto.

    Args:
        U (np.array): Coordenadas X (tamaño N).
        V (np.array): Coordenadas Y (tamaño N).
        N (int): Número de puntos.
        cc_type (int): Tipo de condición de contorno (PERIODICA, NATURAL, LIBRE).
        long_c (float, optional): Longitud característica para contorno periódico.

    Returns:
        tuple: (curvatura, tu_seg1, tu_seg2, ds)
               curvatura (np.array): Curvatura en cada punto (tamaño N).
               tu_seg1, tu_seg2 (np.array): Componentes del vector tangente
                                            en cada segmento (tamaño N).
               ds (np.array): Longitudes de los segmentos (tamaño N).
    """
    if N < 3: # Necesitamos al menos 3 puntos para calcular diferencias centradas
        print(f"Error: Se necesitan al menos 3 puntos para calcular la curvatura (se tienen {N}).", file=sys.stderr)
        return None, None, None, None # Devuelve None para indicar error

    if cc_type == LIBRE:
        print("Error: Condición de contorno LIBRE aún no implementada.", file=sys.stderr)
        sys.exit(1) # O manejar de otra forma si prefieres no salir

    # Crear copias extendidas para manejar condiciones de contorno fácilmente
    # U_ext contendrá [U[0], U[1], ..., U[N-1], U_boundary]
    U_ext = np.zeros(N + 1)
    V_ext = np.zeros(N + 1)
    U_ext[:N] = U
    V_ext[:N] = V

    # Aplicar condición de contorno al último punto del array extendido
    if cc_type == PERIODICA:
        if long_c is None:
             print("Error: Se requiere LongC para condición de contorno PERIODICA.", file=sys.stderr)
             sys.exit(1)
        U_ext[N] = long_c + U[0]
        V_ext[N] = V[0]
    elif cc_type == NATURAL:
        # For natural boundary conditions, the curve connects back to the start
        # We need the segment between the last point (N-1) and the first point (0)
        # The extended array helps calculate diff_u[N-1] and diff_v[N-1] correctly
        U_ext[N] = U[0]
        V_ext[N] = V[0]


    # Calcular diferencias entre puntos consecutivos (vectorizado)
    # diff_u[k] = U_ext[k+1] - U_ext[k] for k=0 to N-1
    diff_u = np.diff(U_ext)
    diff_v = np.diff(V_ext)

    # Calcular longitud de cada segmento (Ds en C)
    # ds[k] = longitud del segmento entre el punto k y k+1 (for k=0 to N-1)
    ds = np.sqrt(diff_u**2 + diff_v**2)

    # Evitar división por cero si hay puntos duplicados
    zero_ds_indices = np.where(ds == 0)[0]
    if len(zero_ds_indices) > 0:
        print(f"Advertencia: Se detectaron {len(zero_ds_indices)} puntos duplicados (longitud de segmento cero).", file=sys.stderr)
        # Replace ds=0 with a small value to avoid NaN/Inf in tangent calculation
        ds[zero_ds_indices] = 1e-15 # Valor pequeño

    # Calcular componentes del vector tangente unitario en cada segmento
    # tu_seg1[k] = componente x de la tangente en el segmento k (k=0 to N-1)
    tu_seg1 = diff_u / ds
    tu_seg2 = diff_v / ds

    # Calcular curvatura en cada punto usando diferencia finita centrada
    # Necesitamos TUSeg[k] y TUSeg[k-1], Ds[k] y Ds[k-1]
    # Usamos np.roll para obtener los valores k-1 fácilmente (maneja el wrap-around)
    # For point k, we use segment k (ending at k+1) and segment k-1 (ending at k)
    tu_seg1_prev = np.roll(tu_seg1, 1) # TUSeg1[k-1] (for k=0, this is TUSeg1[N-1])
    tu_seg2_prev = np.roll(tu_seg2, 1) # TUSeg2[k-1] (for k=0, this is TUSeg2[N-1])
    ds_prev = np.roll(ds, 1)          # Ds[k-1]     (for k=0, this is Ds[N-1])

    # Denominador: ds_avg = (Ds[k] + Ds[k-1]) / 2.0
    # ds_avg[k] = (ds[k] + ds_prev[k]) / 2.0
    ds_avg = 0.5 * (ds + ds_prev)

    # Evitar división por cero en ds_avg (can happen if two consecutive segments have length 0)
    zero_ds_avg_indices = np.where(ds_avg == 0)[0]
    if len(zero_ds_avg_indices) > 0:
         print(f"Advertencia: Promedio de longitud de segmento cero detectado en {len(zero_ds_avg_indices)} puntos.", file=sys.stderr)
         ds_avg[zero_ds_avg_indices] = 1e-15

    # Numerador: d(TU)/ds_avg = (TU[k] - TU[k-1]) / ds_avg[k]
    # aux1 = (TUSeg1[k] - TUSeg1[k-1]) / ds_avg[k]
    # aux2 = (TUSeg2[k] - TUSeg2[k-1]) / ds_avg[k]
    # Note: TU[k] refers to the tangent of the segment *starting* at point k.
    # The curvature at point k is related to the change in tangent *around* point k.
    # It uses the tangent of the segment *ending* at k (tu_seg_prev) and the
    # tangent of the segment *starting* at k (tu_seg).
    aux1 = (tu_seg1 - tu_seg1_prev) / ds_avg
    aux2 = (tu_seg2 - tu_seg2_prev) / ds_avg

    # Magnitud de la curvatura (sin signo)
    curv_mag = np.sqrt(aux1**2 + aux2**2)

    # Determinar el signo de la curvatura
    # Cross product of tangent vectors (tu_prev x tu) z-component sign
    # sign = tu_seg1_prev * tu_seg2 - tu_seg2_prev * tu_seg1
    # C condition: TUSeg1[k]*TUSeg2[mod(k-1,N)] <= TUSeg2[k]*TUSeg1[mod(k-1,N)]
    # Rearranged: tu_seg2 * tu_seg1_prev - tu_seg1 * tu_seg2_prev >= 0
    # This is equivalent to -(tu_prev x tu)_z >= 0
    # If true, curvature is positive (counter-clockwise turn), else negative (clockwise turn).
    sign_condition = (tu_seg1 * tu_seg2_prev) <= (tu_seg2 * tu_seg1_prev)
    curvatura = np.where(sign_condition, curv_mag, -curv_mag)

    # --- Término opcional para Sólido de Revolución ---
    if ES_SOLIDO_DE_REVOLUCION == SI:
        print("Aplicando corrección por sólido de revolución...")
        # Termino C: -0.5 * (TUSeg1[k] + TUSeg1[mod(k-1,N)]) / V[k]
        # Equivalente: -0.5 * (tu_seg1 + tu_seg1_prev) / V
        radial_term = np.zeros_like(curvatura)
        # Evitar división por cero o valores muy grandes cuando V es pequeño
        valid_v = np.abs(V) > 1e-9 # Umbral pequeño
        # Check if V has valid entries before using it
        if np.any(valid_v):
            radial_term[valid_v] = -0.5 * (tu_seg1[valid_v] + tu_seg1_prev[valid_v]) / V[valid_v]
        else:
            print("Advertencia: Todos los valores de V son cercanos a cero, no se puede aplicar corrección radial.", file=sys.stderr)

        # Añadir término y promediar (como en C)
        curvatura = 0.5 * (curvatura + radial_term)


    # Return the values calculated for the N original points
    # ds corresponds to segments 0 to N-1
    # tu_seg1, tu_seg2 correspond to segments 0 to N-1
    # curvatura corresponds to points 0 to N-1
    return curvatura, tu_seg1, tu_seg2, ds


def calcular_arco_y_theta(N, ds, curvatura, tu_seg1, tu_seg2):
    """
    Calcula la longitud de arco (S) y el ángulo de Whewell (Theta).

    Args:
        N (int): Número de puntos.
        ds (np.array): Longitudes de los segmentos (tamaño N).
        curvatura (np.array): Curvatura en cada punto (tamaño N).
        tu_seg1 (np.array): Componente X del vector tangente en cada segmento (tamaño N).
        tu_seg2 (np.array): Componente Y del vector tangente en cada segmento (tamaño N).


    Returns:
        tuple: (S, Theta)
               S (np.array): Longitud de arco hasta cada punto (tamaño N).
               Theta (np.array): Ángulo de Whewell en cada punto (tamaño N).
    """
    # Calcular longitud de arco S
    # S[k] = suma de ds[0] hasta ds[k-1]
    S = np.zeros(N)
    S[1:] = np.cumsum(ds[:-1]) # Suma acumulada de las longitudes de segmento hasta N-2

    # Calcular ángulo Theta
    Theta = np.zeros(N)

    # Inicializar Theta[0]
    # Use the tangent vector of the segment starting at point 0 (tu_seg1[0], tu_seg2[0])
    # Theta[0] = atan2(tu_seg2[0], tu_seg1[0])
    # The C code averages the Y component of the first and last *segment* tangents.
    # Let's use atan2 of the average tangent vector around point 0.
    # Average tangent vector components at point 0:
    avg_tx0 = 0.5 * (tu_seg1[0] + np.roll(tu_seg1, 1)[0]) # (tu_seg1[0] + tu_seg1[N-1]) / 2
    avg_ty0 = 0.5 * (tu_seg2[0] + np.roll(tu_seg2, 1)[0]) # (tu_seg2[0] + tu_seg2[N-1]) / 2
    Theta[0] = np.arctan2(avg_ty0, avg_tx0)


    # Calcular Theta[k] integrando la curvatura (Regla del Trapecio)
    # Theta[k] = Theta[k-1] + 0.5 * (Curvatura[k] + Curvatura[k-1]) * ds[k-1]
    # Necesitamos Curvatura[k-1], que obtenemos con np.roll
    # curvatura_prev = np.roll(curvatura, 1) # No se usa directamente en el bucle
    # Necesitamos ds[k-1], que también obtenemos con np.roll
    # ds_prev = np.roll(ds, 1) # No se usa directamente en el bucle

    # Iterar para calcular Theta
    for k in range(1, N):
        ds_step = ds[k-1] # Longitud del segmento entre k-1 y k
        avg_curv = 0.5 * (curvatura[k] + curvatura[k-1]) # Curvatura promedio entre k-1 y k
        Theta[k] = Theta[k-1] + avg_curv * ds_step

    # Desenvolver el ángulo Theta para evitar saltos de 2*pi
    Theta = np.unwrap(Theta)

    return S, Theta

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
        data = np.loadtxt(input_filename, comments='#', delimiter=args.delimiter, skiprows=args.skiprows, usecols=(0,1))
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
    curvatura, tu_seg1, tu_seg2, ds = calcular_curvatura(U, V, N, CCONTORNO, LONG_C if CCONTORNO == PERIODICA else None)

    # Verificar si el cálculo de curvatura falló (e.g., N < 3)
    if curvatura is None:
        print("No se pudo calcular la curvatura. Abortando.")
        sys.exit(1)

    # 3. Calcular longitud de arco (S) y ángulo de Whewell (Theta)
    print("Calculando longitud de arco y ángulo de Whewell...")
    S, Theta = calcular_arco_y_theta(N, ds, curvatura, tu_seg1, tu_seg2)

    # 4. Preparar nombres de archivo de salida
    base_name = os.path.splitext(input_filename)[0]
    # Añadir sufijos si se ordenó o promedió
    order_suffix = "_ordered" if args.order else ""
    avg_suffix = f"_avg{average_n}" if average_n > 1 else ""
    output_whewell = f"{base_name}{order_suffix}{avg_suffix}_Whewell.csv"
    output_cesaro = f"{base_name}{order_suffix}{avg_suffix}_Cesaro.csv"
    output_points = f"{base_name}{order_suffix}{avg_suffix}_Points.csv"

    # 5. Guardar resultados
    try:
        print(f"Guardando ecuación de Whewell en: {output_whewell}")
        np.savetxt(output_whewell, np.vstack((S, Theta)).T, fmt='%.6e', delimiter='\t',
                   header='LongitudArco_S\tAngulo_Theta_rad', comments='')

        print(f"Guardando ecuación de Cesaro en: {output_cesaro}")
        np.savetxt(output_cesaro, np.vstack((S, curvatura)).T, fmt='%.6e', delimiter='\t',
                   header='LongitudArco_S\tCurvatura', comments='')
        
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
