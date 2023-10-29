"""
This module maps tools to their respective classes.

Design Pattern: Factory pattern.

Make sure to update the ToolName enum when adding a new tool.
"""
from enum import Enum
from typing import List, Dict
import os
import shutil
from pathlib import Path
# in-house
from bench.tools import (ToolOp, BioBloomTools, Clark, FACS, PhageFilter, Kraken2)
from bench.utils import BenchmarkResult
from bench import utils



class ToolName(Enum):
    PHAGEFILTER = "PhageFilter"
    BIOBLOOMTOOLS = "BioBloomTools"
    FACS = "FACS"
    CLARK = "Clark"
    KRACKEN = "Kraken2"

    @classmethod
    def create_tools(cls, tools_to_instantiate: List['ToolName'], configuration: dict):
        """Instantiate tools from list of tool names.

        Args:
            tools_to_instantiate (List[ToolName]): List of tool names to instantiate.
            configuration (dict): Dictionary containing tool configuration.

        Returns:
            dict: Dictionary mapping tool names to instantiated tool objects.
        """
        toolname2class = {
            cls.PHAGEFILTER.value: PhageFilter,
            cls.BIOBLOOMTOOLS.value: BioBloomTools,
            cls.FACS.value: FACS,
            cls.CLARK.value: Clark,
            cls.KRACKEN.value: Kraken2
        }

        # create tools mapping
        tools = {}
        for tool_name in tools_to_instantiate:
            if tool_name.value not in configuration:
                raise ValueError(f"Configuration for {tool_name.value} not found.")

            print(tool_name.value, configuration[tool_name.value])
            tool_class = toolname2class[tool_name.value]
            tools[tool_name.value] = tool_class(**configuration[tool_name.value])

        return tools

    @staticmethod
    def build_dbs(tools: Dict[str, ToolOp], configuration, genome_directory: Path) -> Dict[str, BenchmarkResult]:
        """Given a map from name to each tool wrapper, create a db."""
        # build DBs
        toolname2buildresult = {}
        for toolname, tool in tools.items():
            print(f"Building DB for {toolname}...")
            tool_DB = configuration[toolname].get("database_name", None)
            if os.path.exists(tool_DB):
                if os.path.isfile(tool_DB):
                    os.remove(tool_DB)
                else:
                    shutil.rmtree(tool_DB) # delete DB if exists
            tool_build_cmd = tool.build(tool_DB, genome_directory)
            toolname2buildresult[toolname]: BenchmarkResult = utils.run_command(tool_build_cmd)
            
        return toolname2buildresult