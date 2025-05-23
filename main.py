from sys import argv
import os
import subprocess
from src.utils import read_data
from src.simulation_controller import Simulation
from dolfin import set_log_level, LogLevel

if __name__ == "__main__":
    if len(argv) < 2:
        print("Uso: python main.py <archivo_configuracion>")
        exit(1)

    config_file = argv[1]
    config_data = read_data(config_file)
    if config_data is None:
        print(f"Error al leer el archivo de configuración: {config_file}")
        exit(1)

    case_name = config_data.get("name", None)
    if case_name is None:
        print("El archivo de configuración debe contener un campo 'name'.")
        exit(1)

    caseDir = os.path.join("./results", case_name)
    config_data["caseDir"] = caseDir
    mesh_name = config_data.get("mesh_data", {}).get("file_name", "mesh")

    # Si el directorio existe, preguntar confirmación
    if os.path.isdir(caseDir):
        confirm = input(f"El directorio '{caseDir}' ya existe. ¿Continuar y sobrescribir? [y/N]: ").lower()
        if confirm != "y":
            print("Abortando.")
            exit(1)
        # Eliminar el directorio si se confirma
        subprocess.run(["rm", "-rf", caseDir])

    os.makedirs(caseDir, exist_ok=True)

    # Generar la malla con GMSH
    gmsh_cmd = [
        "gmsh", "-2", "-format", "msh2",
        f"meshes/{mesh_name}.geo",
        "-o", f"{caseDir}/{mesh_name}.msh", "-v","0"
    ]
    subprocess.run(gmsh_cmd, check=True)

    # Convertir la malla a XML con dolfin-convert
    dolfin_cmd = [
        "dolfin-convert",
        f"{caseDir}/{mesh_name}.msh",
        f"{caseDir}/{mesh_name}.xml"
    ]
    subprocess.run(dolfin_cmd, check=True)

    # Set FEniCS log level
    set_log_level(LogLevel.ERROR)

    # Leer configuración y ejecutar simulación
    simulation = Simulation(config_data)
    simulation.run()