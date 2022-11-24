"""
Benchmarking module for PhageFilter.

Install:
    1. Yaml - https://pyyaml.org/wiki/PyYAMLDocumentation
    2. Kraken2 - https://github.com/DerrickWood/kraken2/wiki

Examples:
* running genome count benchmarking
```
python benchmarking/bench.py genomecount -g examples/genomes/viral_genome_dir/ -r res_genomes.csv
```

* running parameterization benchmarking
```
python benchmarking/bench.py parameterization -g examples/genomes/viral_genome_dir/ -f examples/test_reads/simulated_reads.fa -r res_parameterization.csv
```

* relative performance benchmarking
```
python benchmarking/bench.py relative_performance -g examples/genomes/viral_genome_dir/ -r res_relative_performance.csv -t examples/test_reads/ -c benchmarking/config.yaml
```
"""
from abc import ABC, abstractmethod
from enum import Enum, auto
import os
import subprocess
import random
import resource
import shutil
import sys
import tempfile
import time
from typing import List, Tuple, Dict
from pathlib import Path
from dataclasses import dataclass
import argparse
import yaml

GENOME_SUBDIR = "genome/"


class ToolName(Enum):
    PhageFilter = auto(),
    Kraken2 = auto()


@dataclass
class BenchmarkResult:
    # Time in nano-seconds
    elapsed_time: int
    # Memory in bytes?
    max_memory: int

class Experiment:
    def __init__(self, num_genomes: int, source_genomes_dir: Path):
        # Keep track of the seed for reproducability
        self.seed = random.randint(0, sys.maxsize)
        random.seed(self.seed)

        # Make a temporary directory that can hold the randomly sampled genome files
        self.tmp_dir = tempfile.mkdtemp()

        genome_files = list(source_genomes_dir.iterdir())
        selected_genomes = random.sample(genome_files, num_genomes)
        print(self.genome_dir())

        # make genome subdirectory before copying into it
        os.mkdir(self.genome_dir())
        for file in selected_genomes:
            shutil.copy2(file, self.genome_dir())

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        for file in os.listdir(self.genome_dir()):
            os.unlink(os.path.join(self.genome_dir(), file))
        os.rmdir(self.genome_dir())
        os.rmdir(self.tmp_dir)

    def genome_dir(self) -> str:
        return os.path.join(self.tmp_dir, GENOME_SUBDIR)


def run_command(arguments: List) -> BenchmarkResult:
    """
    run a subcommand from the command line.
    returns the running time and memory, along with the output
    path.
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
    name2counts = {}
    with open(fasta_read_path) as read_file:
        line = read_file.readline()
        count = 0
        while line:
            if line[0] == ">":
                name = "_".join(line.strip(">").strip("\n").split("_")[0:-1])
                if name in name2counts:
                    name2counts[name] += 1
                else:
                    name2counts[name] = 1
            line = read_file.readline()
    return name2counts


def get_classification_metrics(true_map, out_map):
    """
    uses true map and PhageFilter map to obtain metrics
    of the classification accuracy.
    """
    TP, FP, FN = (0, 0, 0)
    for genome in out_map.keys():
        if genome in true_map.keys():
            TP += 1
        else:
            FP += 1
    for genome in true_map.keys():
        if genome not in out_map.keys():
            FN += 1
    # get recall and precision
    recall = TP / (TP + FN)
    precision = TP / (TP + FP)
    return recall, precision

class ToolOp(ABC):
    """
    Defines a tool and sub-operations.
    """

    @abstractmethod
    def parse_output():
        """
        parses an output file/directory (depends on tool)
        """
        pass

    @abstractmethod
    def build():
        """
        run tool, based on input arguments, it outputs a CMD-line array.
        """
        pass

    @abstractmethod
    def run():
        """
        run tool, based on input arguments, it outputs a CMD-line array.
        """
        pass


class PhageFilter(ToolOp):

    def __init__(self, kmer_size: int, filter_thresh: float, threads=4):
        """_summary_

        Args:
            kmer_size (int): _description_
            filter_thresh (float): _description_
            threads (int, optional): _description_. Defaults to 4.
        """
        self.k = kmer_size
        self.theta = filter_thresh
        self.threads = threads
        self.db_path = None

    def parse_output(self, output_path: Path, genomes_path: Path=None) -> Dict[str, int]:
        """_summary_
        parses an output file/directory (depends on tool)
        returns a dictionary of the output of PhageFilter.

        Args:
            output_path (Path): Path where the output of PhageFilter
                                will be stored.

        Returns:
            Dict[str, int]: A map from NCBI ID to read count
        """
        name2counts = {}
        with open(output_path) as out_file:
            line = out_file.readline()
            count = 0
            while line:
                name, count = line.strip("\n").split(",")
                name2counts[name] = int(count)
                line = out_file.readline()
        return name2counts

    def build(self, db_path: Path, genomes_path: Path):
        """
        run tool, based on input arguments, it outputs a CMD-line array.
        """
        build_cmd = ["cargo", "run", "--", "build"]
        build_cmd += ["--genomes", f"{genomes_path}"]
        build_cmd += ["--db-path", f"{db_path}"]
        build_cmd += ["--kmer-size", f"{self.k}"]
        build_cmd += ["-t", f"{self.threads}"]
        self.db_path = db_path
        return [build_cmd]

    def run(self, fasta_file: Path, output_path: Path):
        """_summary_
        run tool, based on input arguments, it outputs a CMD-line array.

        Args:
            fasta_file (Path): Path to simualted reads.
            output_path (Path): Desired path for the output file.

        Returns:
            N/A.
        """
        if not self.db_path:
            print("Must first build (PhageFilter)")
            exit()
        run_cmd = ["cargo", "run", "--", "query"]
        run_cmd += ["--reads", f"{fasta_file}"]
        run_cmd += ["--db-path", f"{self.db_path}"]
        run_cmd += ["--cuttoff-threshold", f"{self.theta}"]
        run_cmd += ["--out", f"{output_path}"]
        run_cmd += ["-t", f"{self.threads}"]
        return [run_cmd]


class Kraken2(ToolOp):

    def __init__(self, kmer_size: int = 35, threads=4):
        """_summary_

        Args:
            kmer_size (int): _description_
            filter_thresh (float): _description_
            threads (int, optional): _description_. Defaults to 4.
        """
        self.k = kmer_size
        self.threads = threads
        self.db_path = None

    def parse_output(self, output_path: Path, genomes_path: Path) -> Dict[str, int]:
        """_summary_
        parses an output file/directory (depends on tool)
        returns a dictionary of the output of PhageFilter.

        Args:
            output_path (Path): Path where the output of PhageFilter
                                will be stored.

        Returns:
            Dict[str, int]: A map from NCBI ID to read count
        """
        taxid2ncbi = self.get_taxid2ncbi(genomes_path)
        name2counts = {}
        with open(output_path) as out_file:
            line = out_file.readline()
            count = 0
            while line:
                count, tax_level, taxid = line.strip("\n").split("\t")[2:5]
                if tax_level == "S" or tax_level == "S1" and int(count) > 0:
                    if taxid in taxid2ncbi:
                        for ncbi_id in taxid2ncbi[taxid]:
                            name2counts[ncbi_id] = int(count)
                line = out_file.readline()
        return name2counts

    def build(self, db_path: Path, genomes_path: Path, minimizer_len: int = 31, minimizer_spacing: int = 7):
        """
        run tool, based on input arguments, it outputs a CMD-line array.
        """
        build_cmds = []
        # get map from NCBI ID to taxid
        self.taxid2ncbi = self.get_taxid2ncbi(genomes_path)

        # make DB
        build_cmds.append(["mkdir", "-p", f"{db_path}/taxonomy/"])

        # download NCBI taxonomy
        taxdump_ftp='https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip'
        if not os.path.exists("new_taxdump.zip"):
            build_cmds.append(["wget", f"{taxdump_ftp}"])
        build_cmds.append(["unzip", "new_taxdump.zip", "-d", f"{db_path}taxonomy/"])

        # download NCBI accession2taxid
        nucl_wgs_ftp='https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        if not os.path.exists("nucl_gb.accession2taxid.gz"):
            build_cmds.append(["wget", f"{nucl_wgs_ftp}"])
        build_cmds.append(["gzip", "-d", "-k", "nucl_gb.accession2taxid.gz"])
        build_cmds.append(["mv", "nucl_gb.accession2taxid", f"{db_path}taxonomy/nucl_gb.accession2taxid"])

        # add genomic files
        for genome in os.listdir(genomes_path):
            build_cmd = ["kraken2-build", "--add-to-library", f"{genomes_path / Path(genome)}"]
            build_cmd += ["--db", f"{db_path}"]
            build_cmds.append(build_cmd)

        # build command
        build_cmd = ["kraken2-build", "--build"]
        build_cmd += ["--db", f"{db_path}"]
        build_cmd += ["--kmer-len", f"{self.k}"]
        build_cmd += ["--minimizer-len", f"{minimizer_len}"]
        build_cmd += ["--minimizer-spaces", f"{minimizer_spacing}"]
        build_cmd += ["--fast-build"]
        build_cmds.append(build_cmd)

        # update DB path
        self.db_path = db_path

        return build_cmds

    def run(self, fasta_file: Path, output_path: Path):
        """_summary_
        run tool, based on input arguments, it outputs a CMD-line array.

        Args:
            fasta_file (Path): Path to simualted reads.
            output_path (Path): Desired path for the output file.

        Returns:
            N/A.
        """
        if not self.db_path:
            print("Must first build (Kraken2)")
            exit()
        run_cmd = ["kraken2", "--db", f"{self.db_path}", f"{fasta_file}", "--report", f"{output_path}"] #"--output", f"{output_path}"
        return [run_cmd]

    @staticmethod
    def get_taxid2ncbi(genomes_path: Path) -> Dict[str, int]:
        """_summary_
        This function is used for mapping taxonomy IDs to NCBI accessions. This is useful when
        comparing the output of Kraken2 to PhageFilter.
        
        Args:
            genomes_path (Path): Path to the directory of genomes used to build the Kraken2 DB

        Returns:
            Dict[str, int]: An output dictionary mapping NCBI taxomony IDs to NCBI IDs
        """
        ncbi2tax = {}
        for genome in os.listdir(genomes_path):
            with open(os.path.join(genomes_path, genome), 'r') as f:
                line = f.readline()
                while line:
                    if line.startswith(">"):
                        ncbi = line.strip(">").strip("\n").split("|kraken:taxid|")[1].strip()
                        taxid = line.strip(">").strip("\n").split(" ")[0].strip()
                        if ncbi in ncbi2tax:
                            ncbi2tax[ncbi].append(taxid)
                        else:
                            ncbi2tax[ncbi] = [taxid]
                    line = f.readline()
        return ncbi2tax

class BenchmarkingTests:
    """_summary_
    Different benchmarking tests are contained within this class (as a module, essentially). 
    """

    @staticmethod
    def benchtest_genomecount(phagefilter: PhageFilter, genome_path: Path, phagefilter_db: Path, result_csv: Path, variation_count=10):
        """_summary_
        Performs benchmarking of PhageFilter to obtain an estimate on the impact on time and memory usage
        of having N genomes in the database.

        Args:
            phagefilter (PhageFilter): instance of PhageFilter
            genome_path (Path): Path to the genome directory
            phagefilter_db (Path): Path to the DB for benchmarking (will rewrite for each combination)
            result_csv (Path): Path to desired output file.
        """
        # segment the number of genomes checked so that 10 variations (or 'variation_count') are tested.
        number_of_genomes = len(os.listdir(genome_path))
        step_size = int(number_of_genomes/variation_count)+1

        # benchmark impact of number of genomes on build.
        genomecount2Result: Dict[int, BenchmarkResult] = {}
        for genome_count in range(step_size, number_of_genomes, step_size):
            with Experiment(genome_count, genome_path) as exp:
                # run tool on tmp build directory
                pf_build_cmd = phagefilter.build(
                    phagefilter_db, exp.genome_dir())
                genomecount2Result[genome_count] = run_command(pf_build_cmd)

        # save to output file
        with open(result_csv, "w+") as output_file:
            output_file.write(f"genome count, time (ns), memory (bytes)\n")
            for genome_count, result in genomecount2Result.items():
                output_file.write(
                    f"{genome_count}, {result.elapsed_time}, {result.max_memory}\n")

    @staticmethod
    def benchtest_parameter_sweep(phagefilter: PhageFilter, input_fasta: Path, phagefilter_db: Path, genome_path: Path, result_csv: Path):
        """_summary_
        Performs a parameter sweep for different combinations of
        kmer size and theta. For each combination it saves the output
        to a given result_csv.

        Args:
            phagefilter (PhageFilter): instance of PhageFilter
            input_fasta (Path): Path to the simulated reads file (Fasta or Fastq)
            phagefilter_db (Path): Path to the DB for benchmarking (will rewrite for each combination)
            genome_path (Path): Path to the genome directory
            result_csv (Path): Path to desired output file.

        Returns:
            N/A
        """
        # perform a parameterization for kmer_size and theta.
        truth_map = get_true_maps(input_fasta)
        result_file = open(result_csv, "w+")
        result_file.write("kmer size, theta, time, memory, recall, precision\n")
        for kmer_size in [10, 15, 20, 25]:
            # build a tree for each kmer_size
            phagefilter.k = kmer_size  # update kmer_size
            pf_build_cmd = phagefilter.build(phagefilter_db, genome_path)
            run_command(pf_build_cmd)
            # test tree on kmer size for different thresholds.
            for theta in [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]:
                output_file = f"phagefilter_{kmer_size}_{theta}.csv"
                # update theta.
                phagefilter.theta = theta
                # query
                pf_run_cmd = phagefilter.run(input_fasta, output_file)
                run_result: BenchmarkResult = run_command(pf_run_cmd)
                # benchmark
                result_map = phagefilter.parse_output(output_file)
                recall, precision = get_classification_metrics(true_map=truth_map, out_map=result_map)
                # save to file
                result_file.write(
                    f"{kmer_size}, {theta}, {run_result.elapsed_time}, {run_result.max_memory}, {recall}, {precision}\n")
        # close result file
        result_file.close()

    @staticmethod
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
        configuration = yaml.load(open(f"{config}", "r"), Loader=yaml.SafeLoader)

        # create tool instances
        kraken2 = Kraken2(kmer_size=configuration["Kraken2"]["kmer_size"])
        phagefilter = PhageFilter(kmer_size=configuration["PhageFilter"]["kmer_size"], filter_thresh=configuration["PhageFilter"]["theta"])
        tools = {"Kraken2": kraken2, "PhageFilter": phagefilter}
        
        # build DBs, if not exists
        for toolname, tool in tools.items():
            tool_DB = configuration[toolname]["database_name"]
            if not os.path.exists(tool_DB):
                tool_build_cmd = tool.build(tool_DB, genome_path)
                tool_build_result: BenchmarkResult = run_command(tool_build_cmd)
            else:
                tool.db_path = tool_DB

        # benchmark on test files.
        with open(result_csv, "w+") as result_file:
            result_file.write("tool name, test name, time, memory, recall, precision\n")
            for test_file in os.listdir(test_directory):
                test_file_path = os.path.join(test_directory, test_file)
                truth_map = get_true_maps(test_file_path)
                test_name = test_file.strip('.fna')
                for tool_name, tool in tools.items():
                    output_path = f"{tool_name}_{test_name}"
                    run_cmd = tool.run(test_file_path, output_path)
                    run_result: BenchmarkResult = run_command(run_cmd)
                    # benchmark
                    result_map = tool.parse_output(output_path, genomes_path=genome_path)
                    recall, precision = get_classification_metrics(true_map=truth_map, out_map=result_map)
                    # save to file
                    result_file.write(f"{tool_name}, {test_name}, {run_result.elapsed_time}, {run_result.max_memory}, {recall}, {precision}\n")

class SubparserNames(Enum):
    parameterization = "parameterization"
    genomecount = "genomecount"
    relative_performance = "relative_performance"


def parseArgs(argv=None) -> argparse.Namespace:
    """
    This method takes in the arguments from the command and performs
    parsing.
    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    """
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(
        help='Choose type of benchmarking', dest='sub_parser')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")

    # if parameterization
    parameterization_parser = subparsers.add_parser(
        SubparserNames.parameterization.value)
    parameterization_parser.add_argument(
        "-f", "--input_fasta", type=Path, help="Path to the simulated reads.", required=True)
    parameterization_parser.add_argument(
        "-g", "--genome_dir", type=Path, help="Path to the genome directory", required=True)
    parameterization_parser.add_argument(
        "-r", "--result_csv", type=Path, help="Path to the result CSV (output)", required=True)
    parameterization_parser.add_argument(
        "-db", "--database_name", default=Path("tree2/"), nargs='?', type=Path, help="path to the DB to use [default 'tree2/']", required=False)
    parameterization_parser.add_argument(
        "-k", "--kmer_size", default=20, nargs='?', type=int, help="size of kmer to use [default 20]", required=False)
    parameterization_parser.add_argument(
        "-t", "--threads", help="number of threads to use [Default 4]", required=False)

    # if genomecount
    genomecount_parser = subparsers.add_parser(
        SubparserNames.genomecount.value)
    genomecount_parser.add_argument(
        "-g", "--genome_dir", type=Path, help="Path to the genome directory", required=True)
    genomecount_parser.add_argument(
        "-r", "--result_csv", type=Path, help="Path to the result CSV (output)", required=True)
    genomecount_parser.add_argument(
        "-db", "--database_name", default=Path("tree1/"), nargs='?', type=Path, help="path to the DB to use [default 'tree2/']", required=False)
    genomecount_parser.add_argument(
        "-k", "--kmer_size", default=20, nargs='?', type=int, help="size of kmer to use [default 20]", required=False)
    genomecount_parser.add_argument(
        "-t", "--threads", default=4, nargs='?', type=int, help="number of threads to use [Default 4]", required=False)

    # relative_performance
    relative_performance_parser = subparsers.add_parser(
        SubparserNames.relative_performance.value)
    relative_performance_parser.add_argument(
        "-g", "--genome_dir", type=Path, help="Path to the genome directory", required=True)
    relative_performance_parser.add_argument(
        "-r", "--result_csv", type=Path, help="Path to the result CSV (output)", required=True)
    relative_performance_parser.add_argument(
        "-t", "--test_directory", type=Path, help="path to the directory with simulated test reads.", required=True)
    relative_performance_parser.add_argument(
        "-c", "--config", type=Path, help="path to the configuration file", required=True)

    return parser.parse_args(argv)


def main():
    print(ToolName.PhageFilter)

    # arguments
    args = parseArgs(sys.argv[1:])

    # run benchmark type specified
    if (args.sub_parser == SubparserNames.parameterization.value):
        print(f"Performing parameterization benchmarking...")
        # create tool interfaces
        phagefilter = PhageFilter(kmer_size=args.kmer_size, filter_thresh=1.0)

        # run test
        BenchmarkingTests.benchtest_parameter_sweep(
            phagefilter, args.input_fasta, args.database_name, args.genome_dir, args.result_csv)

    elif (args.sub_parser == SubparserNames.genomecount.value):
        print(f"Performing genome count benchmarking...")
        # create tool interfaces
        phagefilter = PhageFilter(kmer_size=args.kmer_size, filter_thresh=1.0)

        # run test
        BenchmarkingTests.benchtest_genomecount(phagefilter, args.genome_dir,
                                                args.database_name, args.result_csv)

    elif (args.sub_parser == SubparserNames.relative_performance.value):
        print(f"Performing relative performance benchmarking...")
        # run test
        BenchmarkingTests.benchtest_relative_performance(
            args.genome_dir, args.config, args.result_csv, args.test_directory)

    else:
        print(__doc__)


if __name__ == '__main__':
    main()
