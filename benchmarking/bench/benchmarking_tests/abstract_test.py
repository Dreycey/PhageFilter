""" Abstract class for each benchmarking test.

Description:
    This file contains the abstract base class used by each
    benchmarking test. 
"""
from pathlib import Path
from abc import ABC, abstractmethod
import yaml
from typing import Dict
import os
import shutil
import csv
# in-house
from bench.tools.tool_template import ToolOp
import bench.utils as utils



class BenchmarkStrategy(ABC):

    def load_configuration(self, config_path: Path) -> dict:
        """Load configuration from a YAML file."""
        with open(config_path, "r") as config_file:
            return yaml.safe_load(config_file)

    def append_to_csv(self, filepath, row_dict):
        """Add row and header to output CSV."""
        file_exists = Path(filepath).exists()
        
        with open(filepath, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=row_dict.keys())
            if not file_exists:
                writer.writeheader()
            writer.writerow(row_dict)

    def delete_file_if_exists(self, filepath: str):
        """ Deletes a file if it exists. """
        if os.path.exists(filepath):
            os.remove(filepath)

    @abstractmethod
    def run(self, **kwargs):
        pass
