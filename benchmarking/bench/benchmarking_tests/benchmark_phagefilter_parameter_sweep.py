""" This benchmark performs a parameter sweep for PhageFilter.
"""
# std
from pathlib import Path
import os
import numpy as np

# in-house
from bench.benchmarking_tests.abstract_test import BenchmarkStrategy
from bench.tools import ToolName, PhageFilter
from bench.utils import BenchmarkResult, Experiment
import bench.simulate_reads as simulate_reads
import bench.utils as utils



class ParameterBenchmark(BenchmarkStrategy):

    def __init__(self):
        # parameters used during testing
        self.error_rates = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]
        self.kmer_sizes = range(10, 50, 5)
        self.genome_sizes = [1000, 10000, 100000, 1000000]
        self.theta_range = [theta / 10 for theta in range(10)]
        self.fraction_of_bacteria = [0, 0.5]
        self.cuttoff_range = [0, 0.001, 0.1] # fraction of reads a genome must have mapped compared to max.
        self.number_of_replicates = 1
        self.fraction_of_pos_genomes = 0.1 # fraction of genomes to use per genomes within a DB.
        
        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.PHAGEFILTER]

    def run(self, pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, read_count=1000, variation_count: int = 3):
        """
        parameter suite for manual testing.

        Args:
            pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
            neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
            config (Path): Path to a tool configuration file (default in provided directory)
            result_csv (Path): Path to desired output file.
            read_count (int): The number of reads used in the simulated reads file.
            variation_count (int): Essentially equivalent to step size for the number of genomes tested per buildtime benchmark. 
                                    We want to test intervals between 10 and N genomes for built time, where the number 
                                    between 10 and N is given by variation count.
        """
        configuration = self.load_configuration(config)

        # create tools
        tools = ToolName.create_tools(self.tools_to_instantiate, configuration)

        # calculate number of genomes to test
        num_genomes = len(os.listdir(pos_genome_path))
        num_neg_genomes = len(os.listdir(neg_genome_path))
        step_size = int(num_genomes / variation_count)

        # delete results file if it already exists.
        self.delete_file_if_exists(result_csv)

        for tool_name, tool in tools.items():
            assert(tool_name == ToolName.PHAGEFILTER.value) # should only be phagefilter for this test.
            for replicate in range(self.number_of_replicates):
                for contamination_fraction in self.fraction_of_bacteria:
                    # update read counts based on fraction of bacterial reads
                    neg_read_count = contamination_fraction * read_count
                    pos_read_count = read_count - neg_read_count

                    for genome_count in range(step_size, num_genomes+1, step_size):
                        with Experiment(genome_count, pos_genome_path) as genome_subset:
                            for error_rate in self.error_rates:

                                # create simulated reads
                                number_of_pos_genomes_to_sample = max(int(self.fraction_of_pos_genomes * genome_count), 1)
                                sim_reads_name = f"param_sweep_{replicate}_{genome_count}"
                                pos_reads_path = simulate_reads.multi_simulate(genome_subset.genome_dir(), number_of_pos_genomes_to_sample, pos_read_count, 
                                                                               f"posreads_{sim_reads_name}", error_rate=error_rate, out_dir="pos/")
                                neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, num_neg_genomes, neg_read_count, f"negreads_{sim_reads_name}", error_rate=error_rate, out_dir="neg/")
                                combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"combined_{sim_reads_name}.fastq")

                                for kmer_size in self.kmer_sizes:
                                    for genome_size in self.genome_sizes:
                                        # update kmer and genome size, then rebuild db
                                        tool.k = kmer_size
                                        tool.largest_genome = genome_size
                                        ToolName.build_dbs(tools, configuration, genome_subset.genome_dir())

                                        for theta in self.theta_range:
                                            tool.theta = theta
                                            for cutoff in self.cuttoff_range:
                                                # run phagefilter
                                                output_file = f"phagefilter_{sim_reads_name}_{genome_size}_{kmer_size}_{theta}_{error_rate}"
                                                pf_run_cmd = tool.run(combined_test_path, output_file, filter_reads=True)
                                                run_result: BenchmarkResult = utils.run_command(pf_run_cmd)

                                                # classification performance
                                                truth_map = utils.get_true_maps(pos_reads_path)
                                                result_map = tool.parse_output(output_file, cutoff=cutoff)
                                                classification_recall, classification_precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)
                                                read_count_error = utils.get_readcount_metrics(true_map=truth_map, out_map=result_map)
                                                
                                                # filtering performance
                                                result_map = tool.parse_output(output_file, filter_reads=True)
                                                filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                                                # save results
                                                result_dict = {
                                                    "replicate": replicate,
                                                    "kmer size": kmer_size,
                                                    "theta": theta,
                                                    "cutoff percentage": cutoff,
                                                    "genome size": genome_size,
                                                    "error rate": error_rate,
                                                    "number of genomes": genome_count,
                                                    "number of genomes in sample": number_of_pos_genomes_to_sample,
                                                    "fraction of bacteria": contamination_fraction,
                                                    "read count": read_count,
                                                    "time": run_result.elapsed_time,
                                                    "memory": run_result.max_memory,
                                                    "classification recall": classification_recall,
                                                    "classification precision": classification_precision,
                                                    "filter recall": filter_recall,
                                                    "filter precision": filter_precision,
                                                    "avg read count error": np.average(read_count_error)
                                                }
                                                self.append_to_csv(result_csv, result_dict)
                                                
                                                # clean up
                                                utils.delete_files_with_string(output_file)
                                utils.delete_files_with_string(sim_reads_name)
