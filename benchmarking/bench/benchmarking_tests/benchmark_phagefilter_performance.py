"""
Benchmarking the performance of PhageFilter. Tests PhageFilter on time, memory, and accuracy of classification and filtering.
"""
# std
from pathlib import Path
import os
import yaml
import numpy as np

# in-house
from bench.benchmarking_tests.abstract_test import BenchmarkStrategy
from bench.utils import BenchmarkResult, Experiment
import bench.simulate_reads as simulate_reads
import bench.utils as utils

# tools
from bench.tools import ToolName


class PerformanceBenchmark(BenchmarkStrategy):
    
    def __init__(self) -> None:
        # parameters used during testing
        self.read_counts = [100, 1000, 10000, 100000]

        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.PHAGEFILTER]
        
        
    def run(self, pos_genome_path: Path, neg_genome_path: Path, config, result_csv: Path, variation_count: int = 10, contamination_fraction=0.2, fraction2sample = 0.10):
        """_summary_
        Performs benchmarking of PhageFilter, testing time and memory consumption for the build process AND query process.
        For the query process, differing number of reads and DB sizes will tested in combination.

        Args:
            pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
            neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
            config (Path): Path to a tool configuration file (default in provided directory)
            result_csv (Path): Path to desired output file.
            variation_count (int): Essentially equivalent to step size for the number of genomes tested per buildtime benchmark. 
                                    We want to test intervals between 10 and N genomes for built time, where the number 
                                    between 10 and N is given by variation count.
        """
        # open config
        configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

        # Create tools used in benchmarking and build DBs
        tools = ToolName.create_tools(self.tools_to_instantiate, configuration)
        
        # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
        number_of_genomes = len(os.listdir(pos_genome_path))
        number_of_negative_genomes = len(os.listdir(neg_genome_path))
        step_size = int(number_of_genomes/variation_count)

        # delete results file if it already exists.
        self.delete_file_if_exists(result_csv)

        # benchmark impact of number of genomes on build.
        for genome_count in range(step_size, number_of_genomes+1, step_size):
            with Experiment(genome_count, pos_genome_path) as genome_exp:
                assert(len(os.listdir(genome_exp.genome_dir())) == genome_count)
                # run tool on tmp build directory
                dict2result = ToolName.build_dbs(tools, configuration, genome_exp.genome_dir())
                build_result = dict2result[ToolName.PHAGEFILTER.value]
                
                for read_count in self.read_counts:
                    for tool_name, tool in tools.items():
                        assert(tool_name == ToolName.PHAGEFILTER.value) # should only be phagefilter for this test.

                        # make output path and run tool.
                        test_name = f"performanceTest_{genome_count}genomes_{read_count}reads_"
                        output_path = f"{tool_name}_{test_name}"

                        # get read counts from contamination fraction.
                        neg_read_count = int(contamination_fraction * read_count)
                        pos_read_count = read_count - neg_read_count

                        # use pos and neg genome set to simulate contaminated reads.
                        number_of_genomes_to_sample = int(genome_count * fraction2sample)
                        pos_reads_path = simulate_reads.multi_simulate(genome_exp.genome_dir(), max(number_of_genomes_to_sample, 1), pos_read_count, f"posreads_{test_name}")
                        neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, f"negreads_{test_name}")
                        combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                        # truth should only contain positive results
                        truth_map = utils.get_true_maps(pos_reads_path)

                        # run PhageFilter
                        run_cmd = tool.run(combined_test_path, output_path, filter_reads=True)
                        run_result: BenchmarkResult = utils.run_command(run_cmd)

                        # classification benchmark
                        result_map = tool.parse_output(output_path, genomes_path=genome_exp.genome_dir())
                        classification_recall, classification_precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)
                        read_count_error = utils.get_readcount_metrics(true_map=truth_map, out_map=result_map)

                        # filter benchmark
                        result_map = tool.parse_output(output_path, genomes_path=genome_exp.genome_dir(), filter_reads=True)
                        filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                        # save results
                        result_dict = {
                            "genome count": genome_count,
                            "number of reads": read_count,
                            "build time (ns)": build_result.elapsed_time,
                            "build mem (bytes)": run_result.elapsed_time,
                            "query time (ns)": run_result.max_memory,
                            "query memory (bytes)": build_result.max_memory,
                            "classification recall": classification_recall,
                            "classification precision": classification_precision,
                            "filter recall": filter_recall,
                            "filter precision": filter_precision,
                            "avg. read count error": np.average(read_count_error),
                        }
                        self.append_to_csv(result_csv, result_dict)

                        # remove output
                        utils.delete_files_with_string(test_name)