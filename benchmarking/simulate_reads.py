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
from bench.simulate_reads import simulate_reads, multi_simulate

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
