import argparse
import logging
from pathlib import Path

from dolfin import set_log_level, LogLevel

from src.config import Config
from src.problem import Problem
from src.solver_strategy import StaggeredSolver

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.StreamHandler()
        ]
    )

def main():
    parser = argparse.ArgumentParser(description="Run a phase-field simulation.")

    parser.add_argument(
        "--case_dir",
        type=str,
        default=None,
        help="Directorio de resultados (caso) sobreescribe archivo config.json"
        )

    parser.add_argument(
        "--config",
        type=str,
        default="data/gmsh_1.json",
        help="Path to the configuration JSON file.",
    )
    args = parser.parse_args()

    setup_logging()
    set_log_level(LogLevel.ERROR)

    try:
        # 1. Load Configuration
        config_path = Path(args.config)
        config = Config(config_path, case_dir=args.case_dir)

        # 2. Set up the physical problem
        problem = Problem(config)

        # 3. Instantiate the solver and run the simulation
        solver = StaggeredSolver(problem, config)
        solver.run()

    except FileNotFoundError as e:
        logging.error(e)
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}", exc_info=True)

if __name__ == "__main__":
    main()