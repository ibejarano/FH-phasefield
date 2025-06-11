from sys import argv
import os
from src.utils import read_data
from src.simulation_controller import Simulation
from dolfin import LogLevel, set_log_level
set_log_level(LogLevel.ERROR)  # Configura el nivel de log de DOLFIN a WARNING

import logging

# Configuración básica de logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s]  %(funcName)s: %(message)s",
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("simulation.log"),
    ]
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    if len(argv) < 2:
        logger.info("Uso: python main.py <archivo_configuracion> [clave=valor ...]")
        exit(1)

    config_file = argv[1]
    config_data = read_data(config_file, overrrides=argv[2:])
    if config_data is None:
        logger.info(f"Error al leer el archivo de configuración: {config_file}")
        exit(1)

    case_name = config_data.get("name", None)
    if case_name is None:
        logger.info("El archivo de configuración debe contener un campo 'name'.")
        exit(1)

    caseDir = os.path.join("./results", case_name)
    config_data["caseDir"] = caseDir

    simulation = Simulation(config_data)
    for handler in logger.handlers:
        if hasattr(handler, 'flush'):
            handler.flush()
    simulation.run()