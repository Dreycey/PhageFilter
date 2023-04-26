"""
Description:
    This file contains utility functions and classes for benchmarking genomics tools.

Functions:
    run_command: a function that runs a subcommand from the command line, measuring its runtime and memory usage.
    get_true_maps: a function that returns a dictionary of true genomes and read counts given a fasta path.
    get_classification_metrics: a function that uses true and predicted maps to obtain metrics of the classification accuracy.
    get_readcount_metrics: a function that determines the absolute difference in true read counts and a tool's predicted read counts.

Classes:
    Experiment: a class that selects a random set of genomes from a given directory, stores them in a temporary directory, and provides a method for accessing the genome directory.
    BenchmarkResult: a data class that stores the elapsed time and max memory usage of a command.
"""
# standard libraries
import os
from pathlib import Path
from typing import List, Tuple, Dict
import subprocess
import resource
import sys
import time
from multiprocessing import Pool
import random
import tempfile
import shutil
from collections import defaultdict
import logging
# third party libraries
from dataclasses import dataclass




class Experiment:
    def __init__(self, num_genomes: int, source_genomes_dir: Path, tmp_name = "genome/"):
        # save tmp name
        self.tmp_name = tmp_name

        # Keep track of the seed for reproducability
        self.seed = random.randint(0, sys.maxsize)
        random.seed(self.seed)

        # create a tmp directory if it doesn't exist
        tmp_directory_path = "./benchmarking/tmp/"
        if not os.path.exists(tmp_directory_path):
            os.makedirs(tmp_directory_path)
            
        # Make a temporary directory that can hold the randomly sampled genome files
        self.tmp_dir = tempfile.mkdtemp(dir=tmp_directory_path)
        print(self.tmp_dir)
        genome_files = list(source_genomes_dir.iterdir())
        selected_genomes = random.sample(genome_files, num_genomes)

        # make genome subdirectory before copying into it
        os.mkdir(self.genome_dir())

        # copy genome files
        [shutil.copy2(file, self.genome_dir()) for file in selected_genomes]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        for file in os.listdir(self.genome_dir()):
            os.unlink(os.path.join(self.genome_dir(), file))
        os.rmdir(self.genome_dir())
        os.rmdir(self.tmp_dir)

    def genome_dir(self) -> Path:
        return Path(os.path.join(self.tmp_dir, self.tmp_name))
    
@dataclass
class BenchmarkResult:
    # Time in nano-seconds
    elapsed_time: int
    # Memory in bytes?
    max_memory: int


def run_command(arguments: List) -> BenchmarkResult:
    """_summary_
    run a subcommand from the command line.
    returns the running time and memory, along with the output
    path.
    Each subcommand is wrapped in a new process to ensure memory is properly measured.

    Args:
        arguments (List): A nested list of commands to run. Must be nested since some input arguments
                          may have multiple commands.

    Returns:
        BenchmarkResult: a dataclass containing timing and memory information.
    """
    # Run each experiment in a different process so we only measure the max memory usage for that process
    with Pool(1) as pool:
        return pool.apply(_run_command, (arguments,))


def _run_command(arguments: List) -> BenchmarkResult:
    """_summary_
    Actually runs a subcommand from the command line.
    returns the running time and memory, along with the output
    path.

    Args:
        arguments (List): A nested list of commands to run. Must be nested since some input arguments
                          may have multiple commands.

    Returns:
        BenchmarkResult: a dataclass containing timing and memory information.
    """
    start_time = time.monotonic_ns()
    for command in arguments:
        result = subprocess.run(command)
    end_time = time.monotonic_ns()
    end_usage = resource.getrusage(resource.RUSAGE_CHILDREN)

    # Checks return code, raises an error if there's a non-zero return code
    result.check_returncode()

    print(f"Got final resource usage for subprocess: {end_usage}")

    elapsed_time = end_time - start_time
    # Max resident set size used in kilobytes (of the largest child, not the maximum RSS of the process tree).
    # Alternatively could use valgrind massif for a potentially more comprehensive measurement.
    used_memory = end_usage.ru_maxrss

    return BenchmarkResult(elapsed_time, used_memory)


def get_true_maps(fasta_read_path: Path) -> Dict[str,  int]:
    """_summary_
    Given a fasta path, returns a dictionary of true genomes and read counts

    Args:
        fasta_read_path (Path): Path to fasta file with simulated reads 

    Returns:
        Dict[str,  int]: A map from NCBI ID to read count
    """
    name2counts = defaultdict(int)

    with open(fasta_read_path) as read_file:
        for line in read_file:
            if line[0] == "@":
                name = "_".join(line.strip("@").strip("\n").split("_")[0:-1])
                name2counts[name] += 1

    return name2counts

def parse_fasta(file_name):
    """
    fasta to genome string.
    """
    print(f"Reading fasta file: {file_name}")
    genome = ""
    with open(file_name) as f:
        line = f.readline()
        count = 0
        while line:
            if count > 0:
                genome += line.strip("\n")
            else:
                name = line.strip(">").strip("\n").split(" ")[0]
            line = f.readline()
            count += 1
    return genome, name

def compute_metrics(TP: int, FP: int, FN: int) -> Dict[str, float]:
    """
    Computes precision and recall metrics given a True Positive count (TP), False Positive count (FP),
    and False Negative count (FN).
    """
    assert TP >= 0, "True Positive count cannot be negative"
    assert FP >= 0, "False Positive count cannot be negative"
    assert FN >= 0, "False Negative count cannot be negative"
    recall = TP / (TP + FN) if TP + FN != 0 else 0
    precision = TP / (TP + FP) if TP + FP != 0 else 0
    return {
        'recall': recall,
        'precision': precision,
    }

def get_filter_metric_counts(true_map: Dict[str, int], out_map: Dict[str, int]) -> Dict[str, int]:
    """
    Given a true map, mapping genome names to number of read counts, and an out_map doing the same,
    the True Positives, False Positives, and False Negatives for filtered reads may be calculated.
    They are returned in a dictionary with the keys 'TP', 'FP', and 'FN' respectively.
    """
    # get true positive count
    # We don't count hits above the true count to count as True Positives
    TP = 0
    for hit in true_map.keys():
        TP += min(out_map.get(hit, 0), true_map[hit])
        if out_map.get(hit, 0) > true_map[hit]:
            logging.warning(f"Count of {hit} in the out_map ({out_map.get(hit,0)}) exceeded the count in the true_map"
                            f"({true_map[hit]}). This should not happen for read filtering.")
    # get false positive count
    # We make sure to exclude "negative" FP counts from subtracting from the total FP count
    FP = sum([max(0, out_map[hit] - true_map.get(hit, 0)) for hit in out_map.keys()])
    # get false negative count
    # We also make sure to exclude "negative" FN counts from subtracting from the total FN count
    FN = sum([max(0, true_map[hit] - out_map.get(hit, 0)) for hit in true_map.keys()])  # diff
    return {
        'TP': TP,
        'FP': FP,
        'FN': FN,
    }

def get_filter_metrics(true_map: Dict[str, int], out_map: Dict[str, int]) -> Tuple[float, float]:
    """
    Given a true map, mapping genome names to number of read counts, and an out_map doing the same,
    the recall and precision for filtered reads may be calculated.
    """
    counts = get_filter_metric_counts(true_map, out_map)
    metrics = compute_metrics(counts['TP'], counts['FP'], counts['FN'])
    return metrics['recall'], metrics['precision']

def get_classification_metric_counts(true_map: Dict[str, int], out_map: Dict[str, int]) -> Dict[str, int]:
    """
    Given a true map, mapping genome names to number of read counts, and an out_map doing the same,
    the True Positives, False Positives, and False Negatives for classified reads may be calculated.
    They are returned in a dictionary with the keys 'TP', 'FP', and 'FN' respectively. A genome is considered
    classified if it is detected in at least one read.
    """
    TP = len(true_map.keys() & out_map.keys())
    FP = len(out_map.keys() - true_map.keys())
    FN = len(true_map.keys() - out_map.keys())
    return {
        'TP': TP,
        'FP': FP,
        'FN': FN,
    }

def get_classification_metrics(true_map: Dict[str, int], out_map: Dict[str, int]) -> Tuple[float, float]:
    """_summary_
    uses true map and PhageFilter map to obtain metrics
    of the classification accuracy.

    Args:
        true_map (Dict[str, int]): a map of the true genomes and read counts.
        out_map (Dict[str, int]): a map of a tools predicted genomes and read counts.

    Returns:
        Tuple[float, float]: An output of recall and precision
    """
    TP = len(true_map.keys() & out_map.keys())
    FP = len(out_map.keys() - true_map.keys())
    FN = len(true_map.keys() - out_map.keys())

    recall = TP / (TP + FN) if TP + FN != 0 else 0
    precision = TP / (TP + FP) if TP + FP != 0 else 0

    return recall, precision


def get_readcount_metrics(true_map: Dict[str, int], out_map: Dict[str, int]) -> List[int]:
    """_summary_
    Method for determining the absolute diffrence in true read
    counts and a tools predicted read counts. This is important
    for understanding incorrectly mapped reads in addition to 
    being able to accurately predict abundances.

    Args:
        true_map (Dict[str, int]): a map of the true genomes and read counts.
        out_map (Dict[str, int]): a map of a tools predicted genomes and read counts.

    Returns:
        int (List[int]): A list of absolute differences between true and predicted
                         read counts for genomes correctly predicted.
    """
    absolute_difference = []
    for genome, count in out_map.items():
        if genome in true_map.keys():
            absolute_difference.append(abs(count - true_map[genome]))
    return absolute_difference

def delete_files_with_string(target_string, search_path="."):
    """
    This method deletes any files or directories containing 
    the target string. This function is useful for cleaning
    up the output from tools.
    """
    for root, dirs, files in os.walk(search_path):
        for name in files:
            if target_string in name:
                file_path = os.path.join(root, name)
                try:
                    os.remove(file_path)
                    print(f"Deleted file: {file_path}")
                except Exception as e:
                    print(f"Error deleting file: {file_path} - {e}")
                    
        for name in dirs:
            if target_string in name:
                dir_path = os.path.join(root, name)
                try:
                    shutil.rmtree(dir_path)
                    print(f"Deleted directory: {dir_path}")
                except Exception as e:
                    print(f"Error deleting directory: {dir_path} - {e}")