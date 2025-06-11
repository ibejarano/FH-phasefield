from sys import argv
import os
import subprocess
from src.utils import read_data

import logging

# Configuración básica de logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s]  %(funcName)s: %(message)s",
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("simulation.log"),  # Descomenta para guardar en archivo
    ]
)

logger = logging.getLogger(__name__)

def actualizar_geo_con_parametros(geo_path, h, h_coarse, H, l_max):
    """
    Reemplaza los valores de gridsize y ref_gridsize en el archivo .geo.
    """
    with open(geo_path, "r") as f:
        lines = f.readlines()

    nuevas_lineas = []
    for line in lines:
        if line.strip().startswith("gridsize"):
            nuevas_lineas.append(f"gridsize = {h_coarse};\n")
        elif line.strip().startswith("ref_gridsize"):
            nuevas_lineas.append(f"ref_gridsize = {h};\n")
        elif line.strip().startswith("H_sup"):
            nuevas_lineas.append(f"H_sup = {H};\n")
        elif line.strip().startswith("dx"):
            nuevas_lineas.append(f"dx = {l_max/2};\n")
        else:
            nuevas_lineas.append(line)

    with open(geo_path, "w") as f:
        f.writelines(nuevas_lineas)

if __name__ == "__main__":

    if len(argv) < 2:
        logger.info("Uso: python main.py <archivo_configuracion>")
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
    mesh_name = config_data.get("mesh_data", {}).get("file_name", "mesh")

    geo_path = f"meshes/{mesh_name}.geo"
    h = config_data.get("h", None)
    h_coarse = config_data.get("h_coarse", None)
    H = config_data.get("H", None)
    l_max = config_data.get("l_max", None)
    if h is not None and h_coarse is not None:
        actualizar_geo_con_parametros(geo_path, h, h_coarse, H, l_max)
    else:
        logger.warning("Advertencia: No se encontraron los parámetros 'h' y 'h_coarse' en el archivo de configuración.")

    subprocess.run(["rm", "-rf", caseDir])
    os.makedirs(caseDir, exist_ok=True)

    # Generar la malla con GMSH
    gmsh_cmd = [
        "gmsh", "-2", "-format", "msh2",
        f"meshes/{mesh_name}.geo",
        "-o", f"{caseDir}/{mesh_name}.msh"
    ]

    subprocess.run(gmsh_cmd, check=True)

    # Convertir la malla a XML con dolfin-convert
    dolfin_cmd = [
        "dolfin-convert",
        f"{caseDir}/{mesh_name}.msh",
        f"{caseDir}/{mesh_name}.xml"
    ]

    subprocess.run(dolfin_cmd, check=True)
