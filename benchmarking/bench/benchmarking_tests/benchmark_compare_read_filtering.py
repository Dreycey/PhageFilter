"""
This benchmark performs a filtering performance testing for several tools.
"""
# std
from pathlib import Path
import os
import yaml
import shutil
from enum import Enum
from typing import List, Dict

# in-house
from bench.benchmarking_tests.abstract_test import BenchmarkStrategy
from bench.utils import BenchmarkResult, Experiment
import bench.simulate_reads as simulate_reads
import bench.utils as utils

# tools
from bench.tools import ToolName




class FilterBenchmark(BenchmarkStrategy):
    
    def __init__(self) -> None:
        self.error_rates = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]
        self.contamination_percentages = [0, 50, 95, 99, 99.5, 99.9]
    
        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.BIOBLOOMTOOLS, ToolName.FACS, ToolName.CLARK, ToolName.PHAGEFILTER]
            
    def run(self, pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, read_count=100000, fraction2sample = 0.10):
        """_summary_
        This function performs the testing of differing filter tools.

        In essence, it tests how different filtering tools perform on differing contamination rates.

        Args:
            pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
            neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
            config (Path): Path to a tool configuration file (default in provided directory)
            result_csv (Path): Path to desired output file.
            read_count (int): The number of reads used in the simulated reads file.
            fraction2sample (float): the fraction of true genomes, from pos_genome_path, to use within a simulated reads file.
        """
        # open config
        configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

        # Create tools used in benchmarking and build DBs
        tools = ToolName.create_tools(self.tools_to_instantiate, configuration)
        ToolName.build_dbs(tools, configuration, pos_genome_path)

        # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
        number_of_genomes = len(os.listdir(pos_genome_path))
        number_of_negative_genomes = len(os.listdir(neg_genome_path))

        # delete results file if it already exists.
        self.delete_file_if_exists(result_csv)
        
        # If the directory exists, delete it (saving testing information)
        test_results_path = Path(result_csv).parent / Path(result_csv).stem
        if os.path.exists(test_results_path):
            shutil.rmtree(test_results_path)
        os.mkdir(test_results_path)

        print("starting benchmarking")
        # benchmark on test files.
        for error_rate in self.error_rates:
            # test differing contamination percentages on tools
            for contamination_percentage in self.contamination_percentages:
                # get read counts from contamination fraction.
                contamination_fraction = contamination_percentage / 100

                # get read counts from contamination fraction.
                neg_read_count = int(contamination_fraction * read_count)
                pos_read_count = read_count - neg_read_count

                # create test name (used in all output for easy cleanup!)
                test_name = f"filter_test_{error_rate}_{contamination_percentage}"

                # use pos and neg genome set to simulate contaminated reads.
                number_of_genomes_to_sample = int(number_of_genomes * fraction2sample)
                pos_reads_path = simulate_reads.multi_simulate(pos_genome_path, max(number_of_genomes_to_sample, 1), pos_read_count, f"posreads_{test_name}", error_rate=error_rate)
                neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, f"negreads_{test_name}", error_rate=error_rate)
                combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                # get true maps
                truth_map = utils.get_true_maps(pos_reads_path)

                # save information about test files.
                utils.save_truth_information(outdir=test_results_path, 
                                            error_rate=error_rate, 
                                            viral_mapping_dictionary=utils.get_true_maps(pos_reads_path), 
                                            bacterial_mapping_dictionary=utils.get_true_maps(neg_reads_path))

                # iterate through tools
                for tool_name, tool in tools.items():
                    print(f"{tool_name}, {test_name}")

                    # output name
                    output_path = f"{tool_name}_{test_name}"

                    # get query command from each tool, then run
                    run_cmd = tool.run(combined_test_path, output_path, filter_reads=True)
                    run_result: BenchmarkResult = utils.run_command(run_cmd)

                    # benchmark
                    result_map = tool.parse_output(output_path, genomes_path=pos_genome_path, filter_reads=True)
                    recall, precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)
                    utils.get_filter_metrics(true_map=truth_map, out_map=result_map)           
                    
                    # save results
                    result_dict = {
                        "tool name": tool_name,
                        "contamination percentage": contamination_percentage,
                        "error rate": error_rate,
                        "time": run_result.elapsed_time,
                        "memory": run_result.max_memory,
                        "recall": recall,
                        "precision": precision,
                    }
                    self.append_to_csv(result_csv, result_dict)

                # remove output
                utils.delete_files_with_string(test_name)