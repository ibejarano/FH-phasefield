from sys import argv
import os
import subprocess
from src.utils import read_data
from src.simulation_controller import Simulation
from dolfin import set_log_level, LogLevel

if __name__ == "__main__":
    if len(argv) < 3:
        print("Uso: python main.py <nombre_del_caso> <archivo_configuracion>")
        exit(1)

    case_name = argv[1]
    config_file = argv[2]
    caseDir = os.path.join("./results", case_name)
    mesh_name = "curving_H"

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

    # Set FEniCS log level
    set_log_level(LogLevel.ERROR)

    # Leer configuración y ejecutar simulación
    config_data = read_data(config_file)
    simulation = Simulation(config_data)
    simulation.run()