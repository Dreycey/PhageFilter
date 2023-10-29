
""" Benchmarking taxonomic classification.

Description:
    This benchmarking method tests several tools on taxonomic classification 
    benchmarking.
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



class TaxonomicBenchmark(BenchmarkStrategy):

    def __init__(self):
        # parameters used during testing
        self.error_rates = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]

        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.KRACKEN, ToolName.BIOBLOOMTOOLS, ToolName.CLARK, ToolName.PHAGEFILTER]

    def run(self, pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, read_count=100000, contamination_fraction=0.01, variation_count = 10):
        """
        This function performs the relative performance benchmarking of PhageFilter,
        showing how accuracy and precision compare to other tools for taxonomic classification.

        Args:
            pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
            neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
            config (Path): Path to a tool configuration file (default in provided directory)
            result_csv (Path): Path to desired output file.
            read_count (int): The number of reads used in the simulated reads file.
            contamination_fraction (int): the percentage of unknown genomes to add to the tested read file. (sampled from neg_genome_path)
        """
        # open config file
        configuration = self.load_configuration(config)

        # Create tools used in benchmarking and build DBs
        tools = ToolName.create_tools(self.tools_to_instantiate, configuration)
        ToolName.build_dbs(tools, configuration, pos_genome_path)

        # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
        number_of_genomes = len(os.listdir(pos_genome_path))
        number_of_negative_genomes = len(os.listdir(neg_genome_path))
        step_size = int(number_of_genomes/variation_count)+1

        # get read counts from contamination fraction.
        neg_read_count = int(contamination_fraction * read_count)
        pos_read_count = read_count - neg_read_count

        # delete results file if it already exists.
        self.delete_file_if_exists(result_csv)

        # If the directory exists, delete it (saving testing information)
        test_results_path = Path(result_csv).parent / Path(result_csv).stem
        if os.path.exists(test_results_path):
            shutil.rmtree(test_results_path)
        os.mkdir(test_results_path)

        # iterate through different genome numbers (number of genomes from DB in metagenomic mixture)
        for genome_count in range(step_size, number_of_genomes, step_size):
            for error_rate in self.error_rates:
                # create test name (used in all output for easy cleanup!)
                test_name = f"relative_performance_benchmarking_{error_rate}_{genome_count}"

                # simulate reads: use pos and neg genome set to simulate contaminated reads.
                pos_reads_path = simulate_reads.multi_simulate(pos_genome_path, genome_count, pos_read_count, f"posreads_{test_name}", error_rate=error_rate)
                neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, f"negreads_{test_name}", error_rate=error_rate)
                combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                # get truth map. The only true genomes should come from the genome dir used to build the tool DBs
                truth_map = utils.get_true_maps(pos_reads_path)

                # save information about test files.
                utils.save_truth_information(outdir=test_results_path, 
                                                error_rate=error_rate, 
                                                viral_mapping_dictionary=utils.get_true_maps(pos_reads_path), 
                                                bacterial_mapping_dictionary=utils.get_true_maps(neg_reads_path))

                # test all tools
                for tool_name, tool in tools.items():
                    output_path = f"{tool_name}_{test_name}"

                    # run the tool (on the contaminated file)
                    run_cmd = tool.run(combined_test_path, output_path)
                    run_result: BenchmarkResult = utils.run_command(run_cmd)

                    # benchmark
                    result_map = tool.parse_output(output_path, genomes_path=pos_genome_path)
                    recall, precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)
                    read_count_error = utils.get_readcount_metrics(true_map=truth_map, out_map=result_map)

                    # save results
                    result_dict = {
                        "tool name": tool_name,
                        "genome count in reads": genome_count,
                        "error rate": error_rate,
                        "time": run_result.elapsed_time,
                        "memory": run_result.max_memory,
                        "recall": recall,
                        "precision": precision,
                        "read count error": np.average(read_count_error)
                    }
                    self.append_to_csv(result_csv, result_dict)

                # remove output
                utils.delete_files_with_string(test_name)