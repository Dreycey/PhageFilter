"""
Description:
    This script is used to generate exact substrings of
    a given genome (in fasta format).

For more information, see:
    python3 benchmarking/scripts/simulate_reads.py -h

Usage:
    python3 benchmarking/scripts/simulate_reads.py single_genome -g <genome path - fasta> --c <number of total reads> -l <read length> -n <number of genomes> -e <error rate> -o <outfile name>
    python3 benchmarking/scripts/simulate_reads.py multi_genome -g <genome directory> -c <number of total reads> -l <read length> -n <number of genomes> -e <error rate> -o <outfile prefix>
Example:
    python3 benchmarking/scripts/simulate_reads.py single_genome -g examples/genomes/viral_genome_dir/ -c 100 -o sim_reads.fq -l 100 -e 0.1
    python3 benchmarking/scripts/simulate_reads.py single_genome -g examples/genomes/viral_genome_dir/GCF_000912115.1_ViralProj219123_genomic.fna -c 100 -o sim_reads.fq -l 100 -e 0.1
"""
import os
import sys
import random
from pathlib import Path
import argparse
from enum import Enum
import tempfile
import shutil

GENOME_SUBDIR = "genomes"


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


def simulate_reads(genome, name, read_count, outfile, readlength=100, error_rate=0.0):
    """
    This does not actually simulate reads, as
    it creates a fasta of perfect substrings. If an error
    rate is specified, it doesn't actually simulate reads (per se), but
    gives imperfect reads modeling error rates.
    """
    print(f"Simulating reads to: {outfile}")
    reads_added = 1
    with open(outfile, "a+") as out:
        while reads_added < read_count+1:
            start = random.randint(0, len(genome)-readlength)
            stop = start + readlength
            read = list(genome[start:stop])
            for base2change in range(0, len(read)):
                if random.uniform(0, 1) < error_rate:
                    #base2change = random.randint(0, len(read)-1)
                    read[base2change] = random.choice(["A", "C", "T", "G"])
            read = "".join(read)
            quality = ''.join(["#"]*len(read))
            out.write(f"@{name}_{reads_added}\n{read}\n+\n{quality}\n")
            reads_added += 1


def multi_simulate(genome_directory, number_of_genomes, read_count, outfile, readlength=100, error_rate=0.0):
    """
    Simulates reads for multiple genomes.
    """
    print(genome_directory, number_of_genomes,
          read_count, outfile, readlength, error_rate)
    outfile = Path(
        str(outfile) + f"_c{read_count}_n{number_of_genomes}_e{error_rate}.fq")
    # create subset of sampled genomes, and simulate reads from each.
    read_count_per_genome = int(read_count / number_of_genomes)
    with Experiment(number_of_genomes, genome_directory) as exp:
        for fasta in os.listdir(exp.genome_dir()):
            full_path = os.path.join(exp.genome_dir(), fasta)
            genome, name = parse_fasta(full_path)
            simulate_reads(genome=genome, name=name, read_count=read_count_per_genome,
                           outfile=outfile, readlength=readlength, error_rate=error_rate)


class SubparserNames(Enum):
    single_genome = "single_genome"
    multi_genome = "multi_genome"


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
        help='Choose type of simulation', dest='sub_parser')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")

    # if single_genome
    single_genome = subparsers.add_parser(
        SubparserNames.single_genome.value)
    single_genome.add_argument(
        "-g", "--genome_path", type=Path, help="Path to the genome (fasta)", required=True)
    single_genome.add_argument(
        "-c", "--read_count", type=int, help="Number of total reads to simulate.", required=True)
    single_genome.add_argument(
        "-o", "--output_file", type=Path, help="path to output file (fasta; single ended)", required=True)
    single_genome.add_argument(
        "-l", "--read_length", default=100, nargs='?', type=int, help="size of the read length [Default: 100]", required=False)
    single_genome.add_argument(
        "-e", "--error", default=0.01, nargs='?', type=float, help="Average error rate per base. [Default: 0.01]", required=False)

    # if multi_genome
    multi_genome = subparsers.add_parser(
        SubparserNames.multi_genome.value)
    multi_genome.add_argument(
        "-g", "--genome_path", type=Path, help="Path to the genome directory (fasta)", required=True)
    multi_genome.add_argument(
        "-c", "--read_count", type=int, help="Number of total reads to simulate.", required=True)
    multi_genome.add_argument(
        "-o", "--output_file", type=Path, help="path to output file (fasta; single ended)", required=True)
    multi_genome.add_argument(
        "-n", "--number_of_genomes", type=int, help="number of genomes to use.", required=True)
    multi_genome.add_argument(
        "-l", "--read_length", default=100, nargs='?', type=int, help="size of the read length [Default: 100]", required=False)
    multi_genome.add_argument(
        "-e", "--error", default=0.01, nargs='?', type=float, help="Average error rate per base. [Default: 0.01]", required=False)

    return parser.parse_args(argv)


if __name__ == "__main__":
    # arguments
    args = parseArgs(sys.argv[1:])

    # run simulation type specified
    if (args.sub_parser == SubparserNames.single_genome.value):
        print(f"Performing single genome simulation...")
        # simulate reads - single genome.
        genome, name = parse_fasta(args.genome_path)
        simulate_reads(genome=genome, name=name, read_count=args.read_count,
                       outfile=args.output_file, readlength=args.read_length, error_rate=args.error)

    elif (args.sub_parser == SubparserNames.multi_genome.value):
        print(f"Performing multi genome simulation...")
        # simulate reads - multiple genomes.
        multi_simulate(genome_directory=args.genome_path, number_of_genomes=args.number_of_genomes,
                       read_count=args.read_count, outfile=args.output_file, readlength=args.read_length, error_rate=args.error)

    else:
        print(__doc__)
