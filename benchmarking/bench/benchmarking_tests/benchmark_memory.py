"""
Benchmarking test for measuring memory consumption between different tools.
"""
# std
from pathlib import Path
import os

# in-house
from bench.benchmarking_tests.abstract_test import BenchmarkStrategy
from bench.utils import BenchmarkResult, Experiment
import bench.simulate_reads as simulate_reads
import bench.utils as utils
from bench.tools import ToolName




class MemoryBenchmark(BenchmarkStrategy):
    
    def __init__(self):
        # parameters used during testing
        self.read_counts = [10000, 100000]
        
        # tools used in this benchmarking
        self.tools_to_instantiate = [ToolName.KRACKEN, ToolName.BIOBLOOMTOOLS, ToolName.CLARK, ToolName.PHAGEFILTER, ToolName.FACS]
    
    def run(self, pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, variation_count: int = 3, contamination_fraction = 0.2):
        """
        This function performs the testing of differing filter tools for time
        and memory consumption.

        Args:
            pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
            neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
            config (Path): Path to a tool configuration file (default in provided directory)
            result_csv (Path): Path to desired output file.
            variation_count (int): Essentially equivalent to step size for the number of genomes tested per buildtime benchmark. 
                                    We want to test intervals between 10 and N genomes for built time, where the number 
                                    between 10 and N is given by variation count.
            contamination_fraction (int): the percentage of unknown genomes to add to the tested read file. (sampled from neg_genome_path)
        """
        # open config
        configuration = self.load_configuration(config)

        # Create tools used in benchmarking and build DBs
        tools = ToolName.create_tools(self.tools_to_instantiate, configuration)
        ToolName.build_dbs(tools, configuration, pos_genome_path)

        # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
        number_of_genomes = len(os.listdir(pos_genome_path))
        number_of_negative_genomes = len(os.listdir(neg_genome_path))
        step_size = int(number_of_genomes/variation_count)+1
        contamination_percentage = int(contamination_fraction * 100)

        # iterate through different genome numbers (for DBs)
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, pos_genome_path, tmp_name=f"tmp_genomes_{genome_count}") as genome_subset:
                # make sure number of genomes is correct
                print(f"genome count (set) = {genome_count}")
                print(f"genome count (actual) = {len(os.listdir(genome_subset.genome_dir()))}")
                assert(len(os.listdir(genome_subset.genome_dir())) == genome_count)

                # build DBs, if not exists
                tool_build_result = ToolName.build_dbs(tools, configuration, genome_subset.genome_dir())

                for read_count in self.read_counts:
                    # get read counts from contamination fraction.
                    neg_read_count = contamination_fraction * read_count
                    pos_read_count = (1 - contamination_fraction) * read_count

                    # create test name
                    test_name = f"filter_test_genomes{genome_count}_reads{read_count}_contam{contamination_percentage}"

                    # use pos and neg genome set to simulate contaminated reads.
                    pos_reads_path = simulate_reads.multi_simulate(Path(genome_subset.genome_dir()), genome_count, pos_read_count, f"posreads_{test_name}", error_rate=0.3)
                    neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, f"negreads_{test_name}", error_rate=0)
                    combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                    # get true maps
                    truth_map = utils.get_true_maps(pos_reads_path)

                    # iterate through tools
                    for tool_name, tool in tools.items():
                        assert(len(os.listdir(genome_subset.genome_dir())) == genome_count)
                        print(f"{tool_name}, {test_name}")
                        output_path = f"{tool_name}_{test_name}"
                        run_cmd = tool.run(combined_test_path, output_path, filter_reads=True)
                        run_result: BenchmarkResult = utils.run_command(run_cmd)

                        # remove output
                        utils.delete_files_with_string(output_path)

                        # save results
                        result_dict = {
                            "tool name": tool_name,
                            "genome count": genome_count,
                            "read count": read_count,
                            "build time": tool_build_result[tool_name].elapsed_time,
                            "build memory": tool_build_result[tool_name].max_memory,
                            "contamination percentage": contamination_percentage,
                            "query time": run_result.elapsed_time,
                            "query memory": run_result.max_memory
                        }
                        self.append_to_csv(result_csv, result_dict)
                        
                    # remove output
                    utils.delete_files_with_string(test_name)