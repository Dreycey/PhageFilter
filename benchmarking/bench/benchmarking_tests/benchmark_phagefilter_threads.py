"""
Testing the impact of using multiple threads with PhageFilter.
"""
# std
from pathlib import Path
import os
import shutil
import numpy as np

# in-house
from bench.benchmarking_tests.abstract_test import BenchmarkStrategy
from bench.utils import BenchmarkResult, Experiment
import bench.simulate_reads as simulate_reads
import bench.utils as utils

# tools
from bench.tools import ToolName




class ThreadsBenchmarking(BenchmarkStrategy):
    
    def __init__(self):
        self.number_of_replicates = 3
        self.thread_counts = [1, 2, 3, 4]
        
        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.PHAGEFILTER]
    
    def run(self, pos_genome_path: Path, neg_genome_path: Path, config, result_csv: Path, variation_count: int = 3, contamination_fraction=0.2, read_count=100000, fraction2sample = 0.10):
        """_summary_
        Performs benchmarking of PhageFilter, testing how threads changes timing.

        Args:
            pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
            neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
            config (Path): Path to a tool configuration file (default in provided directory)
            result_csv (Path): Path to desired output file.
            variation_count (int): Essentially equivalent to step size for the number of genomes tested per buildtime benchmark. 
                                    We want to test intervals between 10 and N genomes for built time, where the number 
                                    between 10 and N is given by variation count.
            contamination_fraction (int): the percentage of unknown genomes to add to the tested read file. (sampled from neg_genome_path)
            fraction2sample (float): the fraction of true genomes, from pos_genome_path, to use within a simulated reads file.
        """
        # open config
        configuration = self.load_configuration(config)

        # Create tools used in benchmarking and build DBs
        tools = ToolName.create_tools(self.tools_to_instantiate, configuration)

        # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
        number_of_genomes = len(os.listdir(pos_genome_path))

        # get read counts from contamination fraction.
        neg_read_count = int(contamination_fraction * read_count)
        pos_read_count = read_count - neg_read_count

        # benchmarking
        for tool_name, tool in tools.items():
            assert(tool_name == ToolName.PHAGEFILTER.value) # should only be phagefilter for this test.
            for replicate in range(0, self.number_of_replicates):
                for thread_count in self.thread_counts:
                    # update threads
                    tool.threads = thread_count

                    # run tool on tmp build directory
                    toolname2result = ToolName.build_dbs(tools, configuration, Path(pos_genome_path))
                    build_result = toolname2result[ToolName.PHAGEFILTER.value]

                    # make test name
                    test_name = f"threads_test_count{read_count}_genomes{number_of_genomes}"

                    # use pos and neg genome set to simulate contaminated reads.
                    number_of_genomes_to_sample = int(number_of_genomes * fraction2sample)
                    pos_reads_path = simulate_reads.multi_simulate(Path(pos_genome_path), max(number_of_genomes_to_sample, 1), pos_read_count, f"posreads_{test_name}")
                    neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, 1, neg_read_count, f"negreads_{test_name}")
                    combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                    # truth should only contain positive results
                    truth_map = utils.get_true_maps(pos_reads_path)

                    # make output path
                    output_path = f"PhageFilter_{test_name}"

                    # run tool
                    run_cmd = tool.run(combined_test_path, output_path, filter_reads=True)
                    run_result: BenchmarkResult = utils.run_command(run_cmd)

                    # benchmark
                    result_map = tool.parse_output(output_path)
                    classification_recall, classification_precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)
                    result_map = tool.parse_output(output_path, filter_reads=True)
                    filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                    # save results
                    result_dict = {
                        "replicate": replicate,
                        "threads": thread_count,
                        "read count": read_count,
                        "build time (ns)": build_result.elapsed_time,
                        "build memory (bytes)": build_result.max_memory,
                        "query time (ns)": run_result.elapsed_time,
                        "query memory (bytes)": run_result.max_memory,
                    }
                    self.append_to_csv(result_csv, result_dict)

                    # remove output
                    utils.delete_files_with_string(test_name)
                    utils.delete_files_with_string(output_path)