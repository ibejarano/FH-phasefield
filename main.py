from sys import argv
import os
from src.utils import read_data
from src.simulation_controller import Simulation
from dolfin import set_log_level, LogLevel

if __name__ == "__main__":
    caseDir = os.path.join("./results", argv[1])
    # Ensure results directory exists
    os.makedirs(caseDir, exist_ok=True)

    # Set FEniCS log level (optional)
    set_log_level(LogLevel.WARNING) # Or INFO, ERROR, etc. 50 is CRITICAL

    # Read configuration data
    config_data = read_data("lab_gmsh") # Assuming 'lab_gmsh.json' in 'data/' folder
    # Create and run the simulation
    simulation = Simulation(config_data)
    simulation.run()