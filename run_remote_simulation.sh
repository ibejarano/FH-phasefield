#!/bin/bash
# --- Configuración ---
REMOTE_USER="ignacio"
REMOTE_HOST="192.168.0.12"
REMOTE_PROJECT_DIR="simulaciones/FH-phasefield"
LOCAL_PROJECT_DIR="." # Directorio actual
MAIN_SIMULATION_SCRIPT="main.py" # El script principal a ejecutar remotamente
# Archivos o directorios adicionales a copiar, separados por espacios
# Ejemplo: FILES_TO_COPY="datos_entrada/ utils/"
FILES_TO_COPY="main.py src/"
RSYNC_OPTS="-avz --delete --exclude results/ --exclude __pycache__/ --exclude *.ipynb --exclude *.csv"
PYTHON_SCRIPT_NAME="main.py"
REMOTE_CONDA_INIT_SCRIPT_PATH=""
CONDA_ENV_NAME="fenicsproject"
REMOTE_CONDA_INIT_SCRIPT_PATH="/home/ignacio/miniconda3/etc/profile.d/conda.sh"

# Casos de datos
SIM_NAME=$1
RESULTS_DIR="results/$1"
MESH_NAME="curving_H"

echo ""
echo "Paso 1: Sincronizando con '${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PROJECT_DIR}/'..."

# La barra al final de LOCAL_DIR es importante para rsync:
# Copia el *contenido* de LOCAL_DIR a REMOTE_FULL_SIM_DIR.
# Si REMOTE_FULL_SIM_DIR no existe, rsync lo creará (si el directorio padre existe y hay permisos).
# Primero, asegurémonos de que el directorio base en el remoto exista (rsync crea el último componente, no toda la ruta)
ssh "${REMOTE_USER}@${REMOTE_HOST}" "mkdir -p ${REMOTE_PROJECT_DIR}"
if [ $? -ne 0 ]; then
  echo "Error: No se pudo crear o asegurar la existencia del directorio base '${REMOTE_PROJECT_DIR}' en el servidor remoto."
fi


rsync ${RSYNC_OPTS} "." "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PROJECT_DIR}/"
if [ $? -eq 0 ]; then
  echo "Sincronización completada exitosamente."
else
  echo "Error durante la sincronización con rsync. Código de salida: $?"
  echo "Por favor, revisa los mensajes de error de rsync."
  exit 1
fi


# --- Paso 2: Ejecutar el script de Python en la máquina remota ---
MESHGEO_DIR=meshes/$MESH_NAME.geo
MESHMSH_DIR=meshes/$MESH_NAME.msh
XML_DIR=$RESULTS_DIR/$MESH_NAME.xml

REMOTE_COMMAND="
echo 'Intentando inicializar Conda...' &&
source \"${REMOTE_CONDA_INIT_SCRIPT_PATH}\" &&
echo 'Conda inicializado. Intentando activar entorno ${CONDA_ENV_NAME}...' &&
conda activate \"${CONDA_ENV_NAME}\" &&
echo 'Entorno Conda activo: \$CONDA_DEFAULT_ENV (Esperado: ${CONDA_ENV_NAME})' &&
echo 'Cambiando al directorio ${REMOTE_PROJECT_DIR}...' &&
cd \"${REMOTE_PROJECT_DIR}\" &&
rm -rf \"${RESULTS_DIR}\" &&
mkdir  -p \"${RESULTS_DIR}\" &&
gmsh -2 -format msh2 \"${MESHGEO_DIR}\" -o \"${MESHMSH_DIR}\" &&
dolfin-convert \"${MESHMSH_DIR}\" \"${XML_DIR}\" &&
echo 'Ejecutando python3 ${PYTHON_SCRIPT_NAME}...' &&
python3 \"${PYTHON_SCRIPT_NAME}\" \"${SIM_NAME}\" &&
echo 'Desactivando entorno Conda...' &&
conda deactivate
"

ssh "${REMOTE_USER}@${REMOTE_HOST}" "${REMOTE_COMMAND}"

# Comprobar el código de salida de la ejecución remota
if [ $? -eq 0 ]; then
  echo "Script de Python ejecutado exitosamente en el servidor remoto."
else
  echo "Error durante la ejecución del script de Python en el servidor remoto. Código de salida: $?"
  echo "Por favor, revisa los mensajes de error del script."
  exit 1
fi

echo ""
echo "Proceso completado."
exit 0