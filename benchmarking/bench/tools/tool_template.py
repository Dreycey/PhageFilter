"""
Description:
    This script provides classes for working with various bioinformatics tools, 
    such as PhageFilter, Kraken2, and FastViromeExplorer. Each class inherits 
    from the abstract base class ToolOp, which defines a common interface for
    parsing output, building the tool's database, and running the tool. The use 
    of an abstract base class enables the implementation of the Template Pattern, 
    which promotes code reusability and modularity.

Classes:
    ToolOp: Abstract base class for bioinformatics tools.
        parse_output: Parse the output of PhageFilter and return a dictionary of NCBI ID to read count.
        build: Build the tool's database.
        run: Run the tool and output command-line arguments.
    PhageFilter: Class for working with the PhageFilter tool.
    Kraken2: Class for working with the Kraken2 tool.
    FastViromeExplorer: Class for working with the FastViromeExplorer tool.
"""

# standard libraries
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Tuple, Dict
import os




class ToolOp(ABC):
    """
    Abstract base class that defines a tool and its sub-operations.
    """

    @abstractmethod
    def parse_output():
        """
        Abstract method to parse an output file/directory (depends on the tool).
        """
        pass

    @abstractmethod
    def build():
        """
        Abstract method to run the tool, based on input arguments, it outputs a command-line array.
        """
        pass

    @abstractmethod
    def run():
        """
        Abstract method to run the tool, based on input arguments, it outputs a command-line array.
        """
        pass