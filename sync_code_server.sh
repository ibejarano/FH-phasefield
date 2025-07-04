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


echo ""
echo "Sincronizando proyecto con '${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PROJECT_DIR}/'..."
rsync ${RSYNC_OPTS} "${LOCAL_PROJECT_DIR}/" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PROJECT_DIR}/"
if [ $? -ne 0 ]; then
  echo "Error durante la sincronización con rsync."
  exit 1
fi


echo "Proceso completado."
exit 0