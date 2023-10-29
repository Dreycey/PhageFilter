"""
PhageFilter testing with read length variation.
"""
# std
from pathlib import Path
import os
from typing import Dict

# in-house
from bench.benchmarking_tests.abstract_test import BenchmarkStrategy
from bench.tools import ToolName
from bench.utils import BenchmarkResult, Experiment
import bench.simulate_reads as simulate_reads
from bench import utils

class ReadLengthBenchmark(BenchmarkStrategy):
    
    def __init__(self):
        
        self.read_lengths = [100, 150, 200, 1000, 10000]
        self.read_counts = [1000, 10000]
        
        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.PHAGEFILTER]

    def run(self, pos_genome_path: Path, neg_genome_path: Path, config, result_csv: Path, variation_count: int = 5, contamination_fraction =0.5, fraction2sample = 0.10):
        """_summary_
        Performs benchmarking of PhageFilter, testing how the read length impacts accuracy.

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
        number_of_negative_genomes = len(os.listdir(neg_genome_path))
        step_size = int(number_of_genomes/variation_count)+1

        # delete results file if it already exists.
        self.delete_file_if_exists(result_csv)

        # benchmark impact of number of genomes on build.
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, pos_genome_path) as genome_subset:
                # build DBs
                ToolName.build_dbs(tools, configuration, genome_subset.genome_dir())
                for read_count in self.read_counts:
                    # get read counts from contamination fraction.
                    neg_read_count = contamination_fraction * read_count
                    pos_read_count = (1 - contamination_fraction) * read_count
                    for read_length in self.read_lengths:
                        for tool_name, tool in tools.items():
                            assert(tool_name == ToolName.PHAGEFILTER.value) # should only be phagefilter for this test.
                            
                            # make test name
                            test_name = f"readlen_test_{read_length}_count{read_count}_genomes{number_of_genomes}"

                            # use pos and neg genome set to simulate contaminated reads.
                            number_of_genomes_to_sample = int(number_of_genomes * fraction2sample)
                            pos_reads_path = simulate_reads.multi_simulate(genome_subset.genome_dir(), max(number_of_genomes_to_sample, 1), pos_read_count, f"posreads_{test_name}", error_rate=0.05, readlength=read_length)
                            neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, f"negreads_{test_name}", error_rate=0, readlength=read_length)
                            combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                            # truth should only contain positive results
                            truth_map = utils.get_true_maps(pos_reads_path)

                            # make output path 
                            output_path = f"PhageFilter_{test_name}"

                            # run tool
                            run_cmd = tool.run(combined_test_path, output_path, filter_reads=True)
                            run_result: BenchmarkResult = utils.run_command(run_cmd)

                            # classification benchmarking
                            result_map = tool.parse_output(output_path)
                            classification_recall, classification_precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)

                            # filtering benchmarking
                            result_map = tool.parse_output(output_path, filter_reads=True)
                            filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                            # save results
                            result_dict = {
                                "genome count": genome_count,
                                "read count": read_count,
                                "read length": read_length,
                                "query time (ns)": run_result.elapsed_time,
                                "query memory (bytes)": run_result.max_memory,
                                "classification recall": classification_recall,
                                "classification precision": classification_precision,
                                "filter recall": filter_recall,
                                "filter precision": filter_precision
                            }
                            self.append_to_csv(result_csv, result_dict)
                            
                            # remove output
                            utils.delete_files_with_string(output_path)
                            utils.delete_files_with_string(test_name)
                            