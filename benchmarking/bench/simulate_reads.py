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
import re
# in house libraries
from bench.utils import Experiment, parse_fasta



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
                    read[base2change] = random.choice(["A", "C", "T", "G"])
            read = "".join(read)
            quality = ''.join(["#"]*len(read))
            out.write(f"@{name}_{reads_added}\n{read}\n+\n{quality}\n")
            reads_added += 1


def multi_simulate(genome_directory: Path, number_of_genomes, read_count, outfile, readlength=100, error_rate=0.0):
    """
    Simulates reads for multiple genomes.
    """
    #print(genome_directory, number_of_genomes,read_count, outfile, readlength, error_rate)
    simulate_reads_genomes = "simread_genomes/"
    outfile = Path(str(outfile) + f"_c{read_count}_n{number_of_genomes}_e{error_rate}.fq")
    # create subset of sampled genomes, and simulate reads from each.
    read_count_per_genome = int(read_count / number_of_genomes)
    with Experiment(number_of_genomes, genome_directory, simulate_reads_genomes) as exp:
        for fasta in os.listdir(exp.genome_dir()):
            full_path = os.path.join(exp.genome_dir(), fasta)
            genome, name = parse_fasta(full_path)
            simulate_reads(genome=genome, name=name, read_count=read_count_per_genome,
                           outfile=outfile, readlength=readlength, error_rate=error_rate)
    return outfile

def get_read_counts(simreads_path: Path):
    """
    This method obtains the read counts from a simulated reads
    file.
    """
    number_of_reads_regex = re.search(r'_c(\d+)_', simreads_path)
    if number_of_reads_regex:
        number_of_reads = int(number_of_reads_regex.group(1))
    else:
        print(f"ERROR: not able to parse read count from {simreads_path}")
        exit(1)

    return number_of_reads