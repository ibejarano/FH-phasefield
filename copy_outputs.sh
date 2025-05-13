#!/bin/bash

# Script para copiar archivos 'output.csv' de subcarpetas a la carpeta raíz
# y renombrarlos con el nombre de la subcarpeta.

# --- Configuración ---
# Ruta a la carpeta raíz que contiene las subcarpetas.
# Puedes pasarla como argumento al script, o descomentar la línea de abajo
# y establecerla aquí directamente si prefieres.
# ROOT_DIR="/ruta/a/tu/carpeta/principal"
# ---------------------

# --- Validar Argumentos o Configuración ---
# Si no se ha configurado ROOT_DIR en el script, usar el primer argumento.
if [ -z "$ROOT_DIR" ]; then
    if [ -z "$1" ]; then
        echo "Uso: $0 <directorio_raiz>"
        echo "Ejemplo: $0 /home/usuario/mis_datos"
        exit 1
    fi
    ROOT_DIR="$1"
fi

# Verificar si la ruta proporcionada es un directorio válido
if [ ! -d "$ROOT_DIR" ]; then
    echo "Error: '$ROOT_DIR' no es un directorio válido."
    exit 1
fi

# Asegurarse de que ROOT_DIR sea una ruta absoluta si es posible (mejora robustez)
ROOT_DIR=$(cd "$ROOT_DIR" && pwd) || { echo "Error al resolver la ruta: '$ROOT_DIR'"; exit 1; }

echo "Iniciando búsqueda y copia en '$ROOT_DIR'..."

# --- Recorrer Subcarpetas ---
# Itera sobre todos los elementos dentro de ROOT_DIR
for item in "$ROOT_DIR"/*; do
    # Verifica si el elemento actual es un directorio
    if [ -d "$item" ]; then
        # Obtiene solo el nombre de la subcarpeta (sin la ruta completa)
        SUBFOLDER_NAME=$(basename "$item")

        # Define la ruta completa al archivo 'output.csv' esperado dentro de la subcarpeta
        SOURCE_FILE="$item/output.csv"

        # Define la ruta de destino y el nuevo nombre para el archivo en la carpeta raíz
        # Se añade .csv al final del nombre de la subcarpeta
        DEST_FILE="$ROOT_DIR/$SUBFOLDER_NAME-output.csv"

        # --- Procesar el Archivo ---
        # Verifica si el archivo 'output.csv' existe en la subcarpeta actual
        if [ -f "$SOURCE_FILE" ]; then
            echo "  Encontrado '$SOURCE_FILE'."
            echo "  Copiando a '$DEST_FILE'..."

            # Copia el archivo, sobrescribiendo si ya existe en el destino
            cp "$SOURCE_FILE" "$DEST_FILE"

            # Verifica el código de salida del comando cp
            if [ $? -eq 0 ]; then
                echo "  Copia exitosa."
            else
                echo "  Error al copiar '$SOURCE_FILE'."
            fi
        #else
            # Opcional: descomenta la siguiente línea si quieres un mensaje
            # para las subcarpetas que no tienen output.csv
            # echo "  '$SOURCE_FILE' no encontrado en '$SUBFOLDER_NAME'. Saltando."
        fi
    fi
done

echo "Proceso completado."