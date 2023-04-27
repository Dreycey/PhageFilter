"""
Description:
    This module contains functions for benchmarking the performance of PhageFilter. The methods included in 
    this module may be long, as they are made to be called directly to perform a specific benchmarking
    goal.

Functions: (in order)
    benchtest_performance_testing
        - This method performs benchmarking on speed and memory usage of PhageFilter for different combinations
          of genomes and reads used with PhageFilter.
    benchtest_threads_testing
        - The benchmarking method tests the speed improvement for build and query based on the number 
          of threads allocated.
    benchtest_read_length_testing
        - Measures the performance of PhageFilter on different read lengths.
    benchtest_genomecount
        - Performs benchmarking of PhageFilter to obtain an estimate on the impact on time and memory usage of 
          having N genomes in the database
    benchtest_parameter_sweep
        - Performs a parameter sweep for different combinations of kmer size and theta. For each combination it saves 
          the output to a given result_csv.
    benchtest_relative_performance
        - Performs the relative classificationperformance benchmarking of PhageFilter, showing how accuracy and precision 
          compare to other tools.
    benchtest_filter_performance
        - This method performs benchmarking of the filter performance between different tools for read filtering.
    benchtest_memory
        - This method compares all tools for memory usage performance.
    benchtest_depth
        - Performs benchmarking on how search depth in the gSBT impacts recall, precision and time.
"""
# standard libraries
import os
from typing import List, Tuple, Dict
from pathlib import Path
import re
import math
import time
# third party libraries
import yaml
import numpy as np
import shutil
# custom libraries
import bench.utils as utils
from bench.utils import BenchmarkResult, Experiment
import bench.simulate_reads as simulate_reads
from bench.tools.kraken2 import Kraken2
from bench.tools.fastviromeexplorer import FastViromeExplorer
from bench.tools.phage_filter import PhageFilter 
from bench.tools.facs import FACS
from bench.tools.biobloomtools import BioBloomTools




def benchtest_performance_testing(phagefilter: PhageFilter, genome_path: Path, phagefilter_db: Path, result_csv: Path, variation_count: int = 10):
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
                genomecount2Result[genome_count] = utils.run_command(pf_build_cmd)
                print(f"before: {genome_dir}")
                for read_count in [100, 1000, 10000, 100000]:
                    # simulate reads and get true genome counts.
                    test_file_path = simulate_reads.multi_simulate(genome_dir, step_size, read_count, "simreads")
                    truth_map = utils.get_true_maps(test_file_path)

                    # print statements.
                    print(f"\n{truth_map}\n")
                    print(len(os.listdir(genome_dir)))
                    
                    # make output path and run tool.
                    test_name = str(test_file_path).strip('.fq')
                    output_path = f"PhageFilter_{test_name}"
                    run_cmd = phagefilter.run(test_file_path, output_path, filter_reads=True)
                    run_result: BenchmarkResult = utils.run_command(run_cmd)

                    # find number of reads from file name (should be in name of simulated read file)
                    number_of_reads = simulate_reads.SimReadParser.get_read_counts(test_name)

                    # benchmark
                    result_map = phagefilter.parse_output(output_path, genomes_path=genome_path)
                    recall, precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)
                    read_count_error = utils.get_readcount_metrics(true_map=truth_map, out_map=result_map)

                    # save results to file
                    output_file.write(f"{genome_count}" + ",")
                    output_file.write(f"{number_of_reads}" + ",")
                    output_file.write(f"{genomecount2Result[genome_count].elapsed_time}" + ",")
                    output_file.write(f"{genomecount2Result[genome_count].max_memory}" + ",")
                    output_file.write(f"{run_result.elapsed_time}" + ",")
                    output_file.write(f"{run_result.max_memory}" + ",")
                    output_file.write(f"{recall}" + ",")
                    output_file.write(f"{precision}" + ",")
                    output_file.write(f"{np.average(read_count_error)}" + "\n")

                    # remove output
                    utils.delete_files_with_string(test_name)


def benchtest_threads_testing(pos_genome_path: Path, neg_genome_path: Path, config, result_csv: Path, variation_count: int = 3, contamination_fraction=0.2, read_count=100000, fraction2sample = 0.10):
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
    configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

    # instantiate tool wrapper/adapter
    phagefilter = PhageFilter(kmer_size=configuration["PhageFilter"]["kmer_size"],
                              filter_thresh=configuration["PhageFilter"]["theta"])
    # name of output tree
    phagefilter_db = configuration["PhageFilter"]["database_name"]

    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(pos_genome_path))
    step_size = int(number_of_genomes/variation_count)+1

    # get read counts from contamination fraction.
    neg_read_count = contamination_fraction * read_count
    pos_read_count = (1 - contamination_fraction) * read_count

    # open the output file and write the header
    with open(result_csv, "w+") as output_file:
        # create output CSV header.
        output_file.write("threads" + ",")
        output_file.write("read count" + ",")
        output_file.write("build time (ns)" + ",")
        output_file.write("build memory (bytes)" + ",")
        output_file.write("query time (ns)" + ",")
        output_file.write("query memory (bytes)" + "\n")

        # benchmark impact of number of genomes on build.
        for replicate in range(0, 3):
            for thread_count in range(1, 5):
                # update threads
                phagefilter.threads = thread_count

                # run tool on tmp build directory
                pf_build_cmd = phagefilter.build(
                    phagefilter_db, Path(pos_genome_path))
                build_result = utils.run_command(pf_build_cmd)

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
                run_cmd = phagefilter.run(combined_test_path, output_path, filter_reads=True)
                run_result: BenchmarkResult = utils.run_command(run_cmd)

                # benchmark
                result_map = phagefilter.parse_output(output_path)
                classification_recall, classification_precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)
                result_map = phagefilter.parse_output(output_path, filter_reads=True)
                filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                # save results to file
                output_file.write(f"{thread_count}" + ",")
                output_file.write(f"{read_count}" + ",")
                output_file.write(f"{build_result.elapsed_time}" + ",")
                output_file.write(f"{build_result.max_memory}" + ",")
                output_file.write(f"{run_result.elapsed_time}" + ",")
                output_file.write(f"{run_result.max_memory}" + "\n")

                # remove output
                utils.delete_files_with_string(test_name)


def benchtest_read_length_testing(pos_genome_path: Path, neg_genome_path: Path, config, result_csv: Path, variation_count: int = 5, contamination_fraction =0.5, fraction2sample = 0.10):
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
    configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

    # instantiate tool wrapper/adapter
    phagefilter = PhageFilter(kmer_size=configuration["PhageFilter"]["kmer_size"], 
                              filter_thresh=configuration["PhageFilter"]["theta"])
    phagefilter_db = configuration["PhageFilter"]["database_name"] # name of output tree

    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(pos_genome_path))
    step_size = int(number_of_genomes/variation_count)+1

    # open the output file and write the header
    with open(result_csv, "w+") as output_file:
        # create output CSV header.
        output_file.write("genome count" + ",")
        output_file.write("read count" + ",")
        output_file.write("read length" + ",")
        output_file.write("query time (ns)" + ",")
        output_file.write("query memory (bytes)" + ",")
        output_file.write("recall (classification)" + ",")
        output_file.write("precision (classification)" + ",")
        output_file.write("recall (filter)" + ",")
        output_file.write("precision (filter)" + "\n")

        # benchmark impact of number of genomes on build.
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, pos_genome_path) as genome_subset:
                # run tool on tmp build directory
                pf_build_cmd = phagefilter.build(phagefilter_db, genome_subset.genome_dir())
                build_result = utils.run_command(pf_build_cmd)

                for read_count in [1000, 10000]:
                    # get read counts from contamination fraction.
                    neg_read_count = contamination_fraction * read_count
                    pos_read_count = (1 - contamination_fraction) * read_count
                    for read_length in range(100, 10000, 1000):
                        # make test name
                        test_name = f"readlen_test_{read_length}_count{read_count}_genomes{number_of_genomes}"

                        # use pos and neg genome set to simulate contaminated reads.
                        number_of_genomes_to_sample = int(number_of_genomes * fraction2sample)
                        pos_reads_path = simulate_reads.multi_simulate(genome_subset.genome_dir(), max(number_of_genomes_to_sample, 1), pos_read_count, f"posreads_{test_name}", error_rate=0.1, readlength=read_length)
                        neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, 1, neg_read_count, f"negreads_{test_name}", error_rate=0, readlength=read_length)
                        combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                        # truth should only contain positive results
                        truth_map = utils.get_true_maps(pos_reads_path)

                        # make output path 
                        output_path = f"PhageFilter_{test_name}"

                        # run tool
                        run_cmd = phagefilter.run(combined_test_path, output_path, filter_reads=True)
                        run_result: BenchmarkResult = utils.run_command(run_cmd)

                        # classification benchmarking
                        result_map = phagefilter.parse_output(output_path)
                        classification_recall, classification_precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)

                        # filtering benchmarking
                        result_map = phagefilter.parse_output(output_path, filter_reads=True)
                        filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                        # save results to file
                        output_file.write(f"{genome_count}" + ",")
                        output_file.write(f"{read_count}" + ",")
                        output_file.write(f"{read_length}" + ",")
                        output_file.write(f"{run_result.elapsed_time}" + ",")
                        output_file.write(f"{run_result.max_memory}" + ",")
                        output_file.write(f"{classification_recall}" + ",")
                        output_file.write(f"{classification_precision}" + ",")
                        output_file.write(f"{filter_recall}" + ",")
                        output_file.write(f"{filter_precision}" + "\n")

                        # remove output
                        utils.delete_files_with_string(test_name)


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
        # create output file CSV header.
        output_file.write(f"cache size" + ",")
        output_file.write(f"genome count" + ",")
        output_file.write(f"time (ns)" + ",")
        output_file.write(f"memory (bytes)" + ",")
        # benchmark impact of number of genomes on build.
        genomecount2Result: Dict[int, BenchmarkResult] = {}
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, genome_path) as exp:
                # test different cache sizes on directory of differing genome count
                for cache_size in range(step_size, genome_count+1, step_size):
                    # run tool on tmp build directory
                    pf_build_cmd = phagefilter.build(phagefilter_db, exp.genome_dir(), cache_size=cache_size)
                    genomecount2Result[genome_count] = utils.run_command(pf_build_cmd)

                    # save the result to the output file
                    output_file.write(f"{cache_size}" + ",")
                    output_file.write(f"{genome_count}" + ",")
                    output_file.write(f"{genomecount2Result[genome_count].elapsed_time}" + ",")
                    output_file.write(f"{genomecount2Result[genome_count].max_memory}" + "\n")


def benchtest_parameter_sweep(pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, read_count=100000, variation_count: int = 3, contamination_fraction=0.5, fraction2sample = 0.10):
    """_summary_
    Performs a parameter sweep for different combinations of
    kmer size and theta. For each combination it saves the output
    to a given result_csv.

    Args:
        pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
        neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
        config (Path): Path to a tool configuration file (default in provided directory)
        result_csv (Path): Path to desired output file.
        read_count (int): The number of reads used in the simulated reads file.
        variation_count (int): Essentially equivalent to step size for the number of genomes tested per buildtime benchmark. 
                                We want to test intervals between 10 and N genomes for built time, where the number 
                                between 10 and N is given by variation count.
        contamination_fraction (int): the percentage of unknown genomes to add to the tested read file. (sampled from neg_genome_path)
        fraction2sample (float): the fraction of true genomes, from pos_genome_path, to use within a simulated reads file.
    """
    # open config
    configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

    # create phagefilter instances
    phagefilter = PhageFilter(kmer_size=configuration["PhageFilter"]["kmer_size"], filter_thresh=configuration["PhageFilter"]["theta"])
    phagefilter_db = configuration["PhageFilter"]["database_name"]

    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(pos_genome_path))
    number_of_negative_genomes = len(os.listdir(neg_genome_path))

    # get read counts from contamination fraction.
    neg_read_count = contamination_fraction * read_count
    pos_read_count = (1 - contamination_fraction) * read_count
    step_size = int(number_of_genomes/variation_count)

    # perform a parameterization for kmer_size and theta.
    with open(result_csv, "w+") as result_file:
        # write output CSV header
        result_file.write("kmer size" + ",")
        result_file.write("theta" + ",")
        result_file.write("error rate" + ",")
        result_file.write("number of genomes" + ",")
        result_file.write("read count" + ",")
        result_file.write("time" + ",")
        result_file.write("memory" + ",")
        result_file.write("classification recall" + ",")
        result_file.write("classification precision" + ",")
        result_file.write("filter recall" + ",")
        result_file.write("filter precision" + ",")
        result_file.write("avg read count error" + "\n")
        
        # iterate through different genome numbers (for DBs)
        for genome_count in range(step_size, number_of_genomes+1, step_size):
            with Experiment(genome_count, pos_genome_path) as genome_subset:
                for kmer_size in range(15, 51, 5):
                    # build a tree for each kmer_size
                    phagefilter.k = kmer_size  # update kmer_size

                    #utils.delete_files_with_string(phagefilter_db) # delete old DB
                    pf_build_cmd = phagefilter.build(phagefilter_db, genome_subset.genome_dir())
                    utils.run_command(pf_build_cmd)

                    # test different error rates
                    for error_rate in [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]:
                        # test tree on kmer size for different thresholds.
                        theta_set = [theta / 10 for theta in range(0, 11, 2)]
                        for theta in theta_set:
                            # update theta.
                            phagefilter.theta = theta

                            # create output name
                            output_file = f"phagefilter_{kmer_size}_{theta}_{error_rate}"
                            print(f"\n\n  {output_file} \n\n")

                            # simulate reads: use pos and neg genome set to simulate contaminated reads.
                            number_of_genomes_to_sample = int(number_of_genomes * fraction2sample)
                            pos_reads_path = simulate_reads.multi_simulate(genome_subset.genome_dir(), max(number_of_genomes_to_sample, 1), pos_read_count, f"posreads_{output_file}", error_rate=error_rate, out_dir="pos/")
                            neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, f"negreads_{output_file}", error_rate=error_rate, out_dir="neg/")
                            combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"combined_{output_file}.fastq")

                            # query
                            pf_run_cmd = phagefilter.run(combined_test_path, output_file, filter_reads=True)
                            run_result: BenchmarkResult = utils.run_command(pf_run_cmd)

                            # benchmark classification
                            truth_map = utils.get_true_maps(pos_reads_path)
                            result_map = phagefilter.parse_output(output_file)
                            classification_recall, classification_precision = utils.get_classification_metrics(true_map=truth_map, out_map=result_map)
                            read_count_error = utils.get_readcount_metrics(true_map=truth_map, out_map=result_map)

                            # benchmark filtering
                            truth_map = utils.get_true_maps(pos_reads_path) # duplicated for consistency.
                            result_map = phagefilter.parse_output(output_file, filter_reads=True)
                            filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                            # remove output
                            utils.delete_files_with_string(output_file)

                            # save to file
                            result_file.write(f"{kmer_size}" + ",")
                            result_file.write(f"{theta}" + ",")
                            result_file.write(f"{error_rate}" + ",")
                            result_file.write(f"{genome_count}" + ",")
                            result_file.write(f"{read_count}" + ",")
                            result_file.write(f"{run_result.elapsed_time}" + ",")
                            result_file.write(f"{run_result.max_memory}" + ",")
                            result_file.write(f"{classification_recall}" + ",")
                            result_file.write(f"{classification_precision}" + ",")
                            result_file.write(f"{filter_recall}" + ",")
                            result_file.write(f"{filter_precision}" + ",")
                            result_file.write(f"{np.average(read_count_error)}" + "\n")


def benchtest_relative_performance(pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, read_count=100000, contamination_fraction=0.5, fraction2sample = 0.10):
    """_summary_
    This function performs the relative performance benchmarking of PhageFilter,
    showing how accuracy and precision compare to other tools.

    Args:
        pos_genome_path (Path): path to genomes used to build DBs (in the case of classification) and to true positive genomes (in the case of filtering)
        neg_genome_path (Path): path to unkown genomes (in the case of classification and filtering)
        config (Path): Path to a tool configuration file (default in provided directory)
        result_csv (Path): Path to desired output file.
        read_count (int): The number of reads used in the simulated reads file.
        contamination_fraction (int): the percentage of unknown genomes to add to the tested read file. (sampled from neg_genome_path)
        fraction2sample (float): the fraction of true genomes, from pos_genome_path, to use within a simulated reads file.
    """
    # open config file
    configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

    # create tool instances
    kraken2 = Kraken2(kmer_size=configuration["Kraken2"]["kmer_size"])
    phagefilter = PhageFilter(kmer_size=configuration["PhageFilter"]["kmer_size"], 
                              filter_thresh=configuration["PhageFilter"]["theta"])
    fve = FastViromeExplorer(tool_path=configuration["FastViromeExplorer"]["tool_path"], 
                             kmer_size=configuration["FastViromeExplorer"]["kmer_size"], 
                             list_file_path=configuration["FastViromeExplorer"]["list_file_path"])
    biobloomtools = BioBloomTools(kmer_size=configuration["BioBloomTools"]["kmer_size"])

    # map from toolname to tool adapter
    tools = {
             "FastViromeExplorer": fve, 
             "PhageFilter": phagefilter,
             "Kraken2": kraken2,
             "BioBloomTools": biobloomtools,
            }

    # build DBs
    for toolname, tool in tools.items():
        tool_DB = configuration[toolname]["database_name"]
        utils.delete_files_with_string(tool_DB.strip("/")) # delete DB if exists
        tool_build_cmd = tool.build(tool_DB, pos_genome_path)
        tool_build_result: BenchmarkResult = utils.run_command(tool_build_cmd)

    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(pos_genome_path))
    number_of_negative_genomes = len(os.listdir(neg_genome_path))

    # get read counts from contamination fraction.
    neg_read_count = contamination_fraction * read_count
    pos_read_count = (1 - contamination_fraction) * read_count

    # benchmark on test files.
    with open(result_csv, "w+") as result_file:
        result_file.write("tool name" + ",")
        result_file.write("test name" + ",")
        result_file.write("time" + ",")
        result_file.write("memory" + ",")
        result_file.write("recall" + ",")
        result_file.write("precision" + ",")
        result_file.write("read count error" + "\n")

        # iterate through different genome numbers (for DBs)
        for genome_count in range(5, number_of_genomes, 5):
            for error_rate in [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]:
                # create test name (used in all output for easy cleanup!)
                test_name = f"relative_performance_benchmarking_{error_rate}_{genome_count}"

                # simulate reads: use pos and neg genome set to simulate contaminated reads.
                number_of_genomes_to_sample = int(number_of_genomes * fraction2sample)
                pos_reads_path = simulate_reads.multi_simulate(pos_genome_path, max(number_of_genomes_to_sample, 1), pos_read_count, f"posreads_{test_name}", error_rate=error_rate)
                neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, f"negreads_{test_name}", error_rate=error_rate)
                combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                # get truth map. The only true genomes should 
                # come from the genome dir used to build the tool DBs
                truth_map = utils.get_true_maps(pos_reads_path)

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

                    # save to file
                    result_file.write(f"{tool_name}" + ",")
                    result_file.write(f"{test_name}" + ",")
                    result_file.write(f"{run_result.elapsed_time}" + ",")
                    result_file.write(f"{run_result.max_memory}" + ",")
                    result_file.write(f"{recall}" + ",")
                    result_file.write(f"{precision}" + ",")
                    result_file.write(f"{np.average(read_count_error)}" + "\n")

                # remove output
                utils.delete_files_with_string(test_name)


def benchtest_filter_performance(pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, read_count=100000, fraction2sample = 0.10):
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

    # create tool instances
    phagefilter = PhageFilter(kmer_size=configuration["PhageFilter"]["kmer_size"], filter_thresh=configuration["PhageFilter"]["theta"])
    facs = FACS(kmer_size=configuration["FACS"]["kmer_size"], filter_thresh=configuration["FACS"]["theta"])
    biobloomtools = BioBloomTools(kmer_size=configuration["BioBloomTools"]["kmer_size"])

    # map from toolname to tool adapter
    tools = {
             "PhageFilter": phagefilter,
             "FACS": facs,
             "BioBloomTools": biobloomtools,
            }

    # build DBs
    for toolname, tool in tools.items():
        tool_DB = configuration[toolname]["database_name"]
        utils.delete_files_with_string(tool_DB.strip("/")) # delete DB if exists
        tool_build_cmd = tool.build(tool_DB, pos_genome_path)
        tool_build_result: BenchmarkResult = utils.run_command(tool_build_cmd)

    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(pos_genome_path))
    number_of_negative_genomes = len(os.listdir(neg_genome_path))

    # benchmark on test files.
    with open(result_csv, "w+") as result_file:
        # create header for the output CSV
        result_file.write("tool name" + ",")
        result_file.write("contamination percentage" + ",")
        result_file.write("error rate" + ",")
        result_file.write("time" + ",")
        result_file.write("memory" + ",")
        result_file.write("recall" + ",")
        result_file.write("precision" + "\n")

        for error_rate in [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]:
            # test differing contamination percentages on tools
            for contamination_percentage in range(0, 101, 20):
                # get read counts from contamination fraction.
                contamination_fraction = contamination_percentage / 100
                neg_read_count = contamination_fraction * read_count
                pos_read_count = (1 - contamination_fraction) * read_count

                # create test name (used in all output for easy cleanup!)
                test_name = f"filter_test_{error_rate}_{contamination_percentage}"

                # use pos and neg genome set to simulate contaminated reads.
                number_of_genomes_to_sample = int(number_of_genomes * fraction2sample)
                pos_reads_path = simulate_reads.multi_simulate(pos_genome_path, max(number_of_genomes_to_sample, 1), pos_read_count, f"posreads_{test_name}", error_rate=error_rate)
                neg_reads_path = simulate_reads.multi_simulate(neg_genome_path, number_of_negative_genomes, neg_read_count, f"negreads_{test_name}", error_rate=error_rate)
                combined_test_path = simulate_reads.combine_files([pos_reads_path, neg_reads_path], f"{test_name}.fastq")

                # get true maps
                truth_map = utils.get_true_maps(pos_reads_path)

                # iterate through tools
                for tool_name, tool in tools.items():
                    print(f"{tool_name}, {test_name}")

                    # output name
                    output_path = f"{tool_name}_{test_name}"

                    # get query command from each tool, then run
                    run_cmd = tool.run(combined_test_path, output_path, filter_reads=True)
                    run_result: BenchmarkResult = utils.run_command(run_cmd)

                    # benchmark
                    result_map = tool.parse_output(output_path, genomes_path=combined_test_path, filter_reads=True)
                    recall, precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)
                    utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                    # save to file
                    result_file.write(f"{tool_name}" + ",")
                    result_file.write(f"{contamination_percentage}" + ",")
                    result_file.write(f"{error_rate}" + ",")
                    result_file.write(f"{run_result.elapsed_time}" + ",")
                    result_file.write(f"{run_result.max_memory}" + ",")
                    result_file.write(f"{recall}" + ",")
                    result_file.write(f"{precision}" + "\n")

                # remove output
                utils.delete_files_with_string(test_name)


def benchtest_memory(pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, variation_count: int = 3, contamination_fraction = 0.2):
    """_summary_
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
    configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

    # create tool instances
    kraken2 = Kraken2(kmer_size=configuration["Kraken2"]["kmer_size"])
    fve = FastViromeExplorer(tool_path=configuration["FastViromeExplorer"]["tool_path"], 
                             kmer_size=configuration["FastViromeExplorer"]["kmer_size"], 
                             list_file_path=configuration["FastViromeExplorer"]["list_file_path"])
    phagefilter = PhageFilter(kmer_size=configuration["PhageFilter"]["kmer_size"], filter_thresh=configuration["PhageFilter"]["theta"])
    facs = FACS(kmer_size=configuration["FACS"]["kmer_size"], filter_thresh=configuration["FACS"]["theta"])
    biobloomtools = BioBloomTools(kmer_size=configuration["BioBloomTools"]["kmer_size"])

    # map from toolname to tool adapter
    tools = {
            "FastViromeExplorer": fve, 
            "Kraken2": kraken2,
            "BioBloomTools": biobloomtools,
            "FACS": facs,
            "PhageFilter": phagefilter,
            }

    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(pos_genome_path))
    number_of_negative_genomes = len(os.listdir(neg_genome_path))
    step_size = int(number_of_genomes/variation_count)+1
    contamination_percentage = int(contamination_fraction * 100)

    # benchmark on test files.
    with open(result_csv, "w+") as result_file:
        # create header for the output CSV
        result_file.write("tool name" + ",")
        result_file.write("genome count" + ",")
        result_file.write("read count" + ",")
        result_file.write("build time" + ",")
        result_file.write("build memory" + ",")
        result_file.write("contamination percentage" + ",")
        result_file.write("query time" + ",")
        result_file.write("query memory" + "\n")

        # iterate through different genome numbers (for DBs)
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, pos_genome_path, tmp_name=f"tmp_genomes_{genome_count}") as genome_subset:
                # make sure number of genomes is correct
                print(f"genome count (set) = {genome_count}")
                print(f"genome count (actual) = {len(os.listdir(genome_subset.genome_dir()))}")
                assert(len(os.listdir(genome_subset.genome_dir())) == genome_count)

                # build DBs, if not exists
                tool_build_result = {}
                for toolname, tool in tools.items():
                    tool_DB = configuration[toolname]["database_name"] # get DB name
                    utils.delete_files_with_string(tool_DB.strip("/")) # delete DB if exists
                    tool_build_cmd = tool.build(tool_DB, genome_subset.genome_dir()) # GET BUILD CMD
                    tool_build_result[toolname]: Dict[str, BenchmarkResult] = utils.run_command(tool_build_cmd) # save mem and time.
                    assert(len(os.listdir(genome_subset.genome_dir())) == genome_count)

                for read_count in [100, 1000]:
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

                        # save to file
                        result_file.write(f"{tool_name}" + ",")
                        result_file.write(f"{genome_count}" + ",")
                        result_file.write(f"{read_count}" + ",")
                        result_file.write(f"{tool_build_result[tool_name].elapsed_time}" + ",")
                        result_file.write(f"{tool_build_result[tool_name].max_memory}" + ",")
                        result_file.write(f"{contamination_percentage}" + ",")
                        result_file.write(f"{run_result.elapsed_time}" + ",")
                        result_file.write(f"{run_result.max_memory}" + "\n")
                    # remove output
                    utils.delete_files_with_string(test_name)


def benchtest_depth(pos_genome_path: Path, neg_genome_path: Path, config: Path, result_csv: Path, read_count=100000, contamination_fraction=0.5, fraction2sample = 0.10):
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
    configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

    # create phagefilter instances
    phagefilter = PhageFilter(kmer_size=configuration["PhageFilter"]["kmer_size"], filter_thresh=configuration["PhageFilter"]["theta"])

    # build DB
    tool_DB = configuration["PhageFilter"]["database_name"]
    tool_build_cmd = phagefilter.build(tool_DB, pos_genome_path)
    tool_build_result: BenchmarkResult = utils.run_command(tool_build_cmd)

    # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
    number_of_genomes = len(os.listdir(pos_genome_path))
    number_of_negative_genomes = len(os.listdir(neg_genome_path))
    depth_of_tree = int(math.log(number_of_genomes, 2))

    # benchmark on test files.
    with open(result_csv, "w+") as result_file:
        # create header for the output CSV
        result_file.write("replicate" + ",")
        result_file.write("tree depth" + ",")
        result_file.write("genome count" + ",")
        result_file.write("error" + ",")
        result_file.write("time" + ",")
        result_file.write("filter recall" + ",")
        result_file.write("filter precision" + "\n")

        for replicates in range(0,3):
            for error_rate in [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]:
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

                for depth in range(0, depth_of_tree+1):
                    # get true maps
                    truth_map = utils.get_true_maps(pos_reads_path)

                    # run tool.
                    output_path = f"{test_name}"
                    run_cmd = phagefilter.run(combined_test_path, output_path, filter_reads=True, depth=depth)
                    run_result: BenchmarkResult = utils.run_command(run_cmd)

                    # benchmark filtering
                    result_map = phagefilter.parse_output(output_path, genomes_path=combined_test_path, filter_reads=True)
                    filter_recall, filter_precision = utils.get_filter_metrics(true_map=truth_map, out_map=result_map)

                    # save to file
                    result_file.write(f"{replicates}" + ",")
                    result_file.write(f"{depth}" + ",")
                    result_file.write(f"{number_of_genomes}" + ",")
                    result_file.write(f"{error_rate}" + ",")
                    result_file.write(f"{run_result.elapsed_time}" + ",")
                    result_file.write(f"{filter_recall}" + ",")
                    result_file.write(f"{filter_precision}" + "\n")

                # remove output
                utils.delete_files_with_string(test_name)