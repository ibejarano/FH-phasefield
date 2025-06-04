import argparse
import os
import subprocess
from src.utils import read_from_json
from src.simulation_controller import Simulation
from src.utils import update_geo
from mpi4py import MPI
import logging

# Argumentos de línea de comandos
parser = argparse.ArgumentParser(description="Simulación Phase-Field FEniCSx")
parser.add_argument("--config", required=True, help="Archivo de configuración JSON")
parser.add_argument("--log", default="INFO", choices=["INFO", "DEBUG"], help="Nivel de logging")
args = parser.parse_args()

# Configuración del logger
log_level = logging.DEBUG if args.log == "DEBUG" else logging.INFO
logging.basicConfig(
    level=log_level,
    format="%(asctime)s [%(levelname)s]  %(funcName)s: %(message)s",
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("simulation.log"),
    ]
)
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        config_file = args.config
        config_data = read_from_json(config_file)
        if config_data is None:
            logger.info(f"Error al leer el archivo de configuración: {config_file}")
            exit(1)

        case_name = config_data.get("name", None)
        if case_name is None:
            logger.info("El archivo de configuración debe contener un campo 'name'.")
            exit(1)

        caseDir = os.path.join("./results", case_name)
        config_data["caseDir"] = caseDir
        mesh_name = config_data.get("mesh_data", {}).get("file_name", "mesh")

        geo_path = f"meshes/{mesh_name}.geo"
        h = config_data.get("h", None)
        h_coarse = config_data.get("h_coarse", None)
        H = config_data.get("H", None)
        if h is not None and h_coarse is not None:
            update_geo(geo_path, h, h_coarse, H)
        else:
            logger.warning("Advertencia: No se encontraron los parámetros 'h' y 'h_coarse' en el archivo de configuración.")

        #if os.path.isdir(caseDir):
            #confirm = input(f"El directorio '{caseDir}' ya existe. ¿Continuar y sobrescribir? [y/N]: ").lower()
            #if confirm != "y":
                #logger.info("Abortando.")
                #exit(1)
            #subprocess.run(["rm", "-rf", caseDir])

        os.makedirs(caseDir, exist_ok=True)
    else:
        config_data = None

    # Broadcast a todos los procesos
    config_data = comm.bcast(config_data, root=0)

    # Leer configuración y ejecutar simulación
    simulation = Simulation(config_data)
    for handler in logger.handlers:
        if hasattr(handler, 'flush'):
            handler.flush()
    simulation.run()