"""
Description:
    This module contains functions for benchmarking the performance of PhageFilter. The three functions included are:

1. benchtest_genomecount:
    - Performs benchmarking of PhageFilter to obtain an estimate on the impact on time and memory usage of 
    having N genomes in the database.

2. benchtest_parameter_sweep:
    - Performs a parameter sweep for different combinations of kmer size and theta. For each combination it saves 
    the output to a given result_csv.

3. benchtest_relative_performance:
    - Performs the relative performance benchmarking of PhageFilter, showing how accuracy and precision compare 
    to other tools.
"""
# standard libraries
import os
from typing import List, Tuple, Dict
from pathlib import Path
import re
# third party libraries
import yaml
import numpy as np
import shutil
# custom libraries
from bench.utils import BenchmarkResult, Experiment, run_command, get_classification_metrics, get_readcount_metrics, get_true_maps
from bench.tools import PhageFilter, Kraken2, FastViromeExplorer
import bench.simulate_reads as simulate_reads



def benchtest_performance_testing(phagefilter: PhageFilter, genome_path: Path, phagefilter_db: Path, result_csv: Path, variation_count: int = 5):
    """_summary_
    Performs benchmarking of PhageFilter, testing time and memory consumption for the build process AND query process.
    For the query process, differing number of reads and DB sizes will tested in combination.

    Args:
        phagefilter (PhageFilter): instance of PhageFilter
        genome_path (Path): Path to the genome directory
        phagefilter_db (Path): Path to the DB for benchmarking (will rewrite for each combination)
        result_csv (Path): Path to desired output file.
        variation_count (int): Essentially equivalent to step size for the number of genomes tested per buildtime benchmark. 
                                We want to test intervals between 10 and N genomes for built time, where the number 
                                between 10 and N is given by variation count.
    """
    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(genome_path))
    step_size = int(number_of_genomes/variation_count)+1

    # open the output file and write the header
    with open(result_csv, "w+") as output_file:
        # create output CSV header.
        output_file.write("genome count" + ",")
        output_file.write("number of reads" + ",")
        output_file.write("build time (ns)" + ",")
        output_file.write("build mem (bytes)" + ",")
        output_file.write("query time (ns)" + ",")
        output_file.write("query memory (bytes)" + ",")
        output_file.write("recall" + ",")
        output_file.write("precision" + ",")
        output_file.write("avg. read count error" + "\n")

        # benchmark impact of number of genomes on build.
        genomecount2Result: Dict[int, BenchmarkResult] = {}
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, genome_path) as genome_exp:
                # run tool on tmp build directory
                genome_dir = Path(genome_exp.genome_dir())
                pf_build_cmd = phagefilter.build(phagefilter_db, genome_dir)
                genomecount2Result[genome_count] = run_command(pf_build_cmd)
                print(f"before: {genome_dir}")
                for read_count in [100, 1000, 10000, 100000]:
                    # simulate reads and get true genome counts.
                    test_file_path = simulate_reads.multi_simulate(genome_dir, step_size, read_count, "simreads")
                    truth_map = get_true_maps(test_file_path)

                    # print statements.
                    print(f"\n{truth_map}\n")
                    print(len(os.listdir(genome_dir)))
                    
                    # make output path and run tool.
                    test_name = str(test_file_path).strip('.fq')
                    output_path = f"PhageFilter_{test_name}"
                    run_cmd = phagefilter.run(test_file_path, output_path, filter_reads=True)
                    run_result: BenchmarkResult = run_command(run_cmd)

                    # find number of reads from file name (should be in name of simulated read file)
                    number_of_reads = simulate_reads.get_read_counts(test_name)

                    # benchmark
                    result_map = phagefilter.parse_output(output_path+"/CLASSIFICATION.csv", genomes_path=genome_path)
                    recall, precision = get_classification_metrics(true_map=truth_map, out_map=result_map)
                    read_count_error = get_readcount_metrics(true_map=truth_map, out_map=result_map)

                    # save results to file
                    output_file.write(f"{genome_count}" + ",")
                    output_file.write(f"{number_of_reads}" + ",")
                    output_file.write(f"{genomecount2Result[genome_count].elapsed_time}" + ",")
                    output_file.write(f"{genomecount2Result[genome_count].max_memory}" + ",")
                    output_file.write(f"{run_result.elapsed_time}" + ",")
                    output_file.write(f"{run_result.max_memory}" + ",")
                    output_file.write(f"{recall}" + ",")
                    output_file.write(f"{precision}" + ",")
                    output_file.write(f"{np.average(read_count_error)}\n")

                    # remove tmp fasta file and output directory
                    os.remove(test_file_path)
                    shutil.rmtree(output_path)

def benchtest_genomecount(phagefilter: PhageFilter, genome_path: Path, phagefilter_db: Path, result_csv: Path, variation_count: int = 10):
    """_summary_
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
    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(genome_path))
    step_size = int(number_of_genomes/variation_count)+1

    # open the output file and write the header
    with open(result_csv, "w+") as output_file:
        output_file.write(f"genome count, time (ns), memory (bytes)\n")
        # benchmark impact of number of genomes on build.
        genomecount2Result: Dict[int, BenchmarkResult] = {}
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, genome_path) as exp:
                # run tool on tmp build directory
                pf_build_cmd = phagefilter.build(phagefilter_db, exp.genome_dir())
                genomecount2Result[genome_count] = run_command(pf_build_cmd)
                # save the result to the output file
                output_file.write(f"{genome_count}, {genomecount2Result[genome_count].elapsed_time}, {genomecount2Result[genome_count].max_memory}\n")

def benchtest_parameter_sweep(phagefilter: PhageFilter, test_directory: Path, phagefilter_db: Path, genome_path: Path, result_csv: Path):
    """_summary_
    Performs a parameter sweep for different combinations of
    kmer size and theta. For each combination it saves the output
    to a given result_csv.

    Args:
        phagefilter (PhageFilter): instance of PhageFilter
        test_directory (Path): Path to directory of test reads (Fasta/Fastq)
        phagefilter_db (Path): Path to the DB for benchmarking (will rewrite for each combination)
        genome_path (Path): Path to the genome directory
        result_csv (Path): Path to desired output file.

    Returns:
        N/A
    """
    # perform a parameterization for kmer_size and theta.
    result_file = open(result_csv, "w+")
    result_file.write(
        "kmer size, theta, error rate, number of genomes, read count, time, memory, recall, precision, avg read count error\n")
    for kmer_size in range(15, 51, 5):
        # build a tree for each kmer_size
        phagefilter.k = kmer_size  # update kmer_size
        pf_build_cmd = phagefilter.build(phagefilter_db, genome_path)
        run_command(pf_build_cmd)
        # test tree on kmer size for different thresholds.
        for theta in range(0, 10, 2):
            theta /= 10
            for test_file in os.listdir(test_directory):
                test_file_path = os.path.join(test_directory, test_file)
                output_file = f"phagefilter_{kmer_size}_{theta}_{test_file}.csv"
                # parse output file (NOTE: will fail if simulated reads file name does not contain values below in name)
                print(test_file)
                error_rate = float(test_file.strip(
                    ".fq").split("_")[-1].strip("e"))
                number_of_genomes = int(test_file.strip(
                    ".fq").split("_")[-2].strip("n"))
                number_of_reads = int(test_file.strip(
                    ".fq").split("_")[-3].strip("c"))
                # update theta.
                phagefilter.theta = theta
                # query
                pf_run_cmd = phagefilter.run(test_file_path, output_file)
                run_result: BenchmarkResult = run_command(pf_run_cmd)
                # benchmark
                truth_map = get_true_maps(test_file_path)
                result_map = phagefilter.parse_output(output_file)
                recall, precision = get_classification_metrics(
                    true_map=truth_map, out_map=result_map)
                read_count_error = get_readcount_metrics(
                    true_map=truth_map, out_map=result_map)
                # save to file
                result_file.write(
                    f"{kmer_size}, {theta}, {error_rate}, {number_of_genomes}, {number_of_reads}, {run_result.elapsed_time}, {run_result.max_memory}, {recall}, {precision}, {np.average(read_count_error)}\n")
    # close result file
    result_file.close()

def benchtest_relative_performance(genome_path: Path, config: Path, result_csv: Path, test_directory: Path):
    """_summary_
    This function performs the relative performance benchmarking of PhageFilter,
    showing how accuracy and precision compare to other tools.

    Args:
        genome_path (Path): Path to the genome directory
        config (Path): Path to a tool configuration file (default in provided directory)
        result_csv (Path): Path to desired output file.
        test_directory (Path): Path to directory of test reads (Fasta/Fastq)
    """
    configuration = yaml.load(
        open(f"{config}", "r"), Loader=yaml.SafeLoader)

    # create tool instances
    kraken2 = Kraken2(kmer_size=configuration["Kraken2"]["kmer_size"])
    phagefilter = PhageFilter(
        kmer_size=configuration["PhageFilter"]["kmer_size"], filter_thresh=configuration["PhageFilter"]["theta"])
    fve = FastViromeExplorer(tool_path=configuration["FastViromeExplorer"]["tool_path"], kmer_size=configuration[
                                "FastViromeExplorer"]["kmer_size"], list_file_path=configuration["FastViromeExplorer"]["list_file_path"])
    tools = {"FastViromeExplorer": fve,
                "Kraken2": kraken2, "PhageFilter": phagefilter}

    # build DBs, if not exists
    for toolname, tool in tools.items():
        tool_DB = configuration[toolname]["database_name"]
        if not os.path.exists(tool_DB):
            tool_build_cmd = tool.build(tool_DB, genome_path)
            tool_build_result: BenchmarkResult = run_command(
                tool_build_cmd)
        else:
            tool.db_path = tool_DB

    # benchmark on test files.
    with open(result_csv, "w+") as result_file:
        result_file.write(
            "tool name, test name, time, memory, recall, precision, read count error\n")
        for test_file in os.listdir(test_directory):
            test_file_path = os.path.join(test_directory, test_file)
            truth_map = get_true_maps(test_file_path)
            test_name = test_file.strip('.fq')
            for tool_name, tool in tools.items():
                output_path = f"{tool_name}_{test_name}"
                run_cmd = tool.run(test_file_path, output_path)
                run_result: BenchmarkResult = run_command(run_cmd)
                # benchmark
                result_map = tool.parse_output(
                    output_path, genomes_path=genome_path)
                recall, precision = get_classification_metrics(
                    true_map=truth_map, out_map=result_map)
                read_count_error = get_readcount_metrics(
                    true_map=truth_map, out_map=result_map)
                # save to file
                result_file.write(
                    f"{tool_name}, {test_name}, {run_result.elapsed_time}, {run_result.max_memory}, {recall}, {precision}, {np.average(read_count_error)}\n")
