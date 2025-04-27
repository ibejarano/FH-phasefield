import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import os

# --- Configuración ---
INPUT_CSV = "simulations_summary_with_coeffs.csv"
Y_VAR_COLUMN = 'coeff_a' # La columna Y está fija a 'coeff_a'

# --- Argument Parser ---
def setup_arg_parser():
    """Configura y devuelve el parser de argumentos de línea de comando."""
    parser = argparse.ArgumentParser(
        description=f"Lee '{INPUT_CSV}', agrupa los datos y plotea una variable vs '{Y_VAR_COLUMN}'."
    )
    parser.add_argument(
        "--group-by",
        required=True,
        type=str,
        metavar='COLUMN_NAME',
        help="Nombre de la columna por la cual agrupar los datos para el plot."
    )
    parser.add_argument(
        "--x-var",
        required=True,
        type=str,
        metavar='COLUMN_NAME',
        help=f"Nombre de la columna a plotear en el eje X (contra '{Y_VAR_COLUMN}' en el eje Y)."
    )
    parser.add_argument(
        "--output-plot",
        type=str,
        metavar='FILENAME.png',
        help="Opcional: Ruta para guardar el gráfico como archivo PNG en lugar de mostrarlo."
    )
    return parser

# --- Función Principal ---
def main():
    parser = setup_arg_parser()
    args = parser.parse_args()

    group_by_col = args.group_by
    x_var_col = args.x_var
    output_plot_path = args.output_plot

    # 1. Leer el archivo CSV
    try:
        print(f"Leyendo el archivo: {INPUT_CSV}")
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print(f"Error: El archivo '{INPUT_CSV}' no fue encontrado en el directorio actual.")
        sys.exit(1)
    except Exception as e:
        print(f"Error al leer el archivo CSV '{INPUT_CSV}': {e}")
        sys.exit(1)

    # 2. Validar que las columnas existan
    required_cols = [group_by_col, x_var_col, Y_VAR_COLUMN]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Las siguientes columnas requeridas no se encontraron en '{INPUT_CSV}': {', '.join(missing_cols)}")
        print(f"Columnas disponibles: {', '.join(df.columns)}")
        sys.exit(1)

    # 3. Preparar datos (opcional: manejar NaNs si es necesario)
    df_cleaned = df.dropna(subset=required_cols)
    if len(df_cleaned) < len(df):
        print(f"Advertencia: Se eliminaron {len(df) - len(df_cleaned)} filas con valores faltantes en las columnas: {', '.join(required_cols)}.")

    if df_cleaned.empty:
        print("Error: No quedan datos válidos para plotear después de limpiar NaNs.")
        sys.exit(1)

    # --- AJUSTE DE TAMAÑO DE FUENTE PARA PRESENTACIÓN ---
    # Actualizar rcParams para fuentes más grandes
    # Puedes ajustar estos valores según tus preferencias
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
    # ----------------------------------------------------

    # 4. Agrupar y Plotear
    print(f"Agrupando por '{group_by_col}' y ploteando '{x_var_col}' vs '{Y_VAR_COLUMN}'...")
    plt.figure(figsize=(10, 7)) # Tamaño del gráfico (puedes ajustarlo también)

    grouped = df_cleaned.groupby(group_by_col)

    for name, group in grouped:
        group_sorted = group.sort_values(by=x_var_col)
        # Intentar formatear como float con 2 decimales, si falla, usar como string
        try:
            label_text = f"{group_by_col} = {float(name):.2f}"
        except (ValueError, TypeError):
            label_text = f"{group_by_col} = {name}" # Usar el valor tal cual si no es numérico

        plt.plot(
            group_sorted[x_var_col],
            group_sorted[Y_VAR_COLUMN],
            marker='o',
            linestyle='-',
            label=label_text # Usar la etiqueta formateada
        )

    # 5. Formatear el Gráfico
    plt.xlabel(x_var_col) # Usará axes.labelsize de rcParams
    plt.ylabel(Y_VAR_COLUMN) # Usará axes.labelsize de rcParams
    #plt.title(f"{Y_VAR_COLUMN} vs {x_var_col} (Agrupado por {group_by_col})") # Usará axes.titlesize
    #plt.legend(title=f"Grupos ({group_by_col})", bbox_to_anchor=(1.05, 1), loc='upper left') # Usará legend.fontsize y legend.title_fontsize
    plt.legend()
    plt.grid(True)
    # Ajustar layout para dejar espacio a la leyenda que está fuera
    # El rect [left, bottom, right, top] define el área del plot dentro de la figura
    plt.tight_layout(rect=[0, 0, 0.80, 1]) # Ajustado para leyenda potencialmente más ancha

    # 6. Mostrar o Guardar el Gráfico
    if output_plot_path:
        try:
            output_dir = os.path.dirname(output_plot_path)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
            # Guardar con mayor DPI para mejor calidad en presentaciones
            plt.savefig(output_plot_path, dpi=300, bbox_inches='tight')
            print(f"Gráfico guardado en: {output_plot_path} (dpi=300)")
        except Exception as e:
            print(f"Error al guardar el gráfico en '{output_plot_path}': {e}")
            print("Mostrando gráfico en pantalla...")
            plt.show()
    else:
        print("Mostrando gráfico...")
        plt.show()

    # Opcional: Restaurar los valores por defecto de rcParams si se van a generar más gráficos después
    # plt.rcdefaults()

# --- Punto de Entrada ---
if __name__ == "__main__":
    main()
