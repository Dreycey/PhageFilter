"""
Testing the impact of using a tunable depth search with PhageFilter.
"""
# std
from pathlib import Path
import os
from typing import Dict
import math

# in-house
from bench.benchmarking_tests.abstract_test import BenchmarkStrategy
from bench.tools import ToolName
from bench.utils import BenchmarkResult, Experiment
import bench.simulate_reads as simulate_reads
from bench import utils




class TunableDepthBenchmark(BenchmarkStrategy):
    
    def __init__(self):
        self.number_of_replicates = 3
        self.error_rates = [0.0, 0.02, 0.04]
        
        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.PHAGEFILTER]
    
    def get_expected_tree_depth(self, genome_count: int):
        """

        Args:
            genome_count (int): The number of genomes in the database.

        Returns:
            int: The expected tree depth for a given number of genomes.
        """
        return int(math.log(genome_count, 2)) + 2
    
    def run(self, pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, read_count=100000, contamination_fraction=0.5, fraction2sample = 0.20):
        """_summary_
        This function performs the testing for the impact of changing depth on both time and filtering performance.

        Args:
            pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
            neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
            config (Path): Path to a tool configuration file (default in provided directory)
            result_csv (Path): Path to desired output file.
            read_count (int): The number of reads used in the simulated reads file.
            contamination_fraction (int): the percentage of unknown genomes to add to the tested read file. (sampled from neg_genome_path)
            fraction2sample (float): the fraction of true genomes, from pos_genome_path, to use within a simulated reads file.
        """
        # open config
        configuration = self.load_configuration(config)

        # Create tools used in benchmarking and build DBs
        tools = ToolName.create_tools(self.tools_to_instantiate, configuration)
        ToolName.build_dbs(tools, configuration, pos_genome_path)

        # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
        number_of_genomes = len(os.listdir(pos_genome_path))
        number_of_negative_genomes = len(os.listdir(neg_genome_path))

        # benchmark on test files.
        for tool_name, tool in tools.items():
            assert(tool_name == ToolName.PHAGEFILTER.value) # should only be phagefilter for this test.
            for replicates in range(0, self.number_of_replicates):
                for error_rate in self.error_rates:
                    # create test name
                    test_name = f"depth_test_{replicates}_{error_rate}"

                    # get read counts from contamination fraction.
                    neg_read_count = contamination_fraction * read_count
                    pos_read_count = (1 - contamination_fraction) * read_count

                    # use pos and neg genome set to simulate contaminated reads.
                    number_of_genomes_to_sample = int(number_of_genomes * fraction2sample)
                    pos_reads_path = simulate_reads.multi_simulate(pos_genome_path, max(number_of_genomes_to_sample, 1), pos_read_count, "posreads", error_rate=error_rate)
                    neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, "negreads", error_rate=error_rate)
                    combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                    for depth in range(0, self.get_expected_tree_depth(number_of_genomes)):
                        # get true maps
                        truth_map = utils.get_true_maps(pos_reads_path)

                        # run tool.
                        output_path = f"{test_name}"
                        run_cmd = tool.run(combined_test_path, output_path, filter_reads=True, depth=depth)
                        run_result: BenchmarkResult = utils.run_command(run_cmd)

                        # benchmark filtering
                        result_map = tool.parse_output(output_path, genomes_path=combined_test_path, filter_reads=True)
                        filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                        # save results
                        result_dict = {
                            "replicate": replicates,
                            "tree depth": depth,
                            "genome count": number_of_genomes,
                            "error": error_rate,
                            "time": run_result.elapsed_time,
                            "filter recall": filter_recall,
                            "filter precision": filter_precision,
                        }
                        self.append_to_csv(result_csv, result_dict)

                    # remove output
                    utils.delete_files_with_string(test_name)