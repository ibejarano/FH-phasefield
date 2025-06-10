#!/bin/bash
# --- Configuración ---
REMOTE_USER="ignacio"
REMOTE_HOST="192.168.0.19"
REMOTE_PROJECT_DIR="repos/FH-phasefield"
LOCAL_PROJECT_DIR="." # Directorio actual
RSYNC_OPTS="-avz --delete --exclude results/ --exclude __pycache__/ --exclude *.ipynb --exclude *.csv"
PYTHON_SCRIPT_NAME="main.py"
REMOTE_CONDA_INIT_SCRIPT_PATH="/home/ignacio/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV_NAME="fenicsproject"


REMOTE_RESULTS_DIR="${REMOTE_PROJECT_DIR}/results"
LOCAL_RESULTS_DIR="./results"
echo "Copiando archivos de resultados desde el servidor (excluyendo .h5 y .xdmf)..."
rsync -avz --exclude '*.h5' --exclude '*.xdmf' "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/" "${LOCAL_RESULTS_DIR}/"

if [ $? -eq 0 ]; then
  echo "Simulación ejecutada exitosamente en el servidor remoto."
  echo "Copiando archivos de resultados desde el servidor (excluyendo .h5 y .xdmf)..."
  rsync -avz --exclude '*.h5' --exclude '*.xdmf'  --exclude '*.xml' "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_RESULTS_DIR}/" "${LOCAL_RESULTS_DIR}/"
else
  echo "Error durante la ejecución remota."
  exit 1
fi

echo "Proceso completado."
exit 0
