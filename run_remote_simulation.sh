#!/bin/bash
# --- Configuración ---
REMOTE_USER="ignacio"
REMOTE_HOST="192.168.0.13"
REMOTE_PROJECT_DIR="simulaciones/FH-phasefield"
LOCAL_PROJECT_DIR="." # Directorio actual
RSYNC_OPTS="-avz --delete --exclude results/ --exclude __pycache__/ --exclude *.ipynb --exclude *.csv"
PYTHON_SCRIPT_NAME="main.py"
REMOTE_CONDA_INIT_SCRIPT_PATH="/home/ignacio/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_NAME="fenicsproject"

# Argumentos: nombre del caso y archivo de configuración
CONFIG_FILE=$1

if [ -z "$CONFIG_FILE" ]; then
  echo "Uso: $0 <nombre_del_caso> <archivo_configuracion>"
  exit 1
fi

echo ""
echo "Sincronizando proyecto con '${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PROJECT_DIR}/'..."
ssh "${REMOTE_USER}@${REMOTE_HOST}" "mkdir -p ${REMOTE_PROJECT_DIR}"
rsync ${RSYNC_OPTS} "${LOCAL_PROJECT_DIR}/" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PROJECT_DIR}/"
if [ $? -ne 0 ]; then
  echo "Error durante la sincronización con rsync."
  exit 1
fi

# --- Ejecutar main.py en el servidor remoto ---
REMOTE_COMMAND="
source \"${REMOTE_CONDA_INIT_SCRIPT_PATH}\" &&
conda activate \"${CONDA_ENV_NAME}\" &&
cd \"${REMOTE_PROJECT_DIR}\" &&
python3 \"${PYTHON_SCRIPT_NAME}\" \"${CONFIG_FILE}\" &&
conda deactivate
"

echo "Ejecutando simulación en el servidor remoto..."
ssh "${REMOTE_USER}@${REMOTE_HOST}" "${REMOTE_COMMAND}"


REMOTE_RESULTS_DIR="${REMOTE_PROJECT_DIR}/results"
LOCAL_RESULTS_DIR="./results"
echo "Copiando archivos de resultados desde el servidor (excluyendo .h5 y .xdmf)..."
rsync -avz --exclude '*.h5' --exclude '*.xdmf' "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/" "${LOCAL_RESULTS_DIR}/"

if [ $? -eq 0 ]; then
  echo "Simulación ejecutada exitosamente en el servidor remoto."
  echo "Copiando archivos de resultados desde el servidor (excluyendo .h5 y .xdmf)..."
  rsync -avz --exclude '*.h5' --exclude '*.xdmf' "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/" "${LOCAL_RESULTS_DIR}/"
else
  echo "Error durante la ejecución remota."
  exit 1
fi

echo "Proceso completado."
exit 0