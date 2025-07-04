import json
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

class Config:
    def __init__(self, config_path: Path, case_dir: str):
        if not isinstance(config_path, Path):
            config_path = Path(config_path)

        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found at {config_path}")

        logger.info(f"Loading configuration from: {config_path}")
        with open(config_path, 'r') as f:
            self.params = json.load(f)
        self.case_dir_str = case_dir

    def get(self, key, default=None):
        return self.params.get(key, default)

    def __getitem__(self, key):
        return self.params[key]

    @property
    def case_dir(self) -> Path:
        return Path(f"results/{self.case_dir_str}")