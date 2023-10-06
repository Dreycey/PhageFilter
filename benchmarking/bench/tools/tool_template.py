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

    @staticmethod
    def get_taxid2ncbi(genomes_path: Path) -> Dict[str, int]:
        """_summary_
        This function is used for mapping taxonomy IDs to NCBI accessions. This is useful when
        comparing the output of Kraken2 to PhageFilter.

        Args:
            genomes_path (Path): Path to the directory of genomes used to build the Kraken2 DB

        Returns:
            Dict[str, int]: An output dictionary mapping NCBI taxomony IDs to NCBI IDs
        """
        taxid2ncbi = {}
        for genome in os.listdir(genomes_path):
            with open(os.path.join(genomes_path, genome), 'r') as f:
                line = f.readline()
                while line:
                    if line.startswith(">"):
                        try:
                            taxid = line.strip(">").strip("\n").split("|kraken:taxid|")[1].strip()
                            ncbi = line.strip(">").strip("\n").split(" ")[0].strip()
                        except:
                            print(f"line: {line}")
                            break # assumes non-multifasta
                        if ncbi in taxid2ncbi:
                            taxid2ncbi[taxid].append(ncbi)
                        else:
                            taxid2ncbi[taxid] = [ncbi]
                    line = f.readline()
        return taxid2ncbi