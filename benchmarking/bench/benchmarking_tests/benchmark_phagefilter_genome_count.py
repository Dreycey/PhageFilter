"""
This benchmark tests how PhageFilter performs when the number of genomes is changed.
"""
# std
from pathlib import Path
import os
from typing import Dict

# in-house
from bench.benchmarking_tests.abstract_test import BenchmarkStrategy
from bench.tools import ToolName
from bench.utils import BenchmarkResult, Experiment




class GenomeCountBenchmark(BenchmarkStrategy):
    
    def __init__(self):
        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.PHAGEFILTER]

    def run(self, pos_genome_path: Path, neg_genome_path: Path, config, result_csv: Path, variation_count: int = 10):
        """
        Performs benchmarking of PhageFilter to obtain an estimate on the impact on time and memory usage
        of having N genomes in the database.

        Args:
            phagefilter (PhageFilter): instance of PhageFilter
            genome_path (Path): Path to the genome directory
            phagefilter_db (Path): Path to the DB for benchmarking (will rewrite for each combination)
            result_csv (Path): Path to desired output file.
            variation_count (int): Essentially equivalent to step size for the number of genomes tested per buildtime benchmark. 
                                    We want to test intervals between 10 and N genomes for built time, where the number 
                                    between 10 and N is given by variation count.
        """
        # open config
        configuration = self.load_configuration(config)
        
        # Create tools used in benchmarking and build DBs
        tools = ToolName.create_tools(self.tools_to_instantiate, configuration)
        
        # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
        number_of_genomes = len(os.listdir(pos_genome_path))
        step_size = int(number_of_genomes/variation_count)+1

        # delete results file if it already exists.
        self.delete_file_if_exists(result_csv)

        # benchmarking
        genomecount2Result: Dict[int, BenchmarkResult] = {} # benchmark impact of number of genomes on build.
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, pos_genome_path) as exp:
                for tool_name, tool in tools.items():
                    assert(tool_name == ToolName.PHAGEFILTER.value) # should only be phagefilter for this test.
                    for cache_size in range(step_size, genome_count+1, step_size):
                        tool.cache_size = cache_size

                        # run tool on tmp build directory
                        dict2result = ToolName.build_dbs(tools, configuration, exp.genome_dir())
                        genomecount2Result[genome_count] = dict2result[ToolName.PHAGEFILTER.value]

                        # save results
                        result_dict = {
                            "cache size": cache_size,
                            "genome count": genome_count,
                            "time (ns)": genomecount2Result[genome_count].elapsed_time,
                            "memory (bytes)": genomecount2Result[genome_count].max_memory
                        }
                        self.append_to_csv(result_csv, result_dict)