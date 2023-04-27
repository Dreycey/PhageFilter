"""
Description:
    This module defines functions for simualting reads. The current implemented method
    does not produce reads simulated from a particular instrument, rather equally probable
    fragments/substrings from the genome with defined error rates. The errors do not simulate
    indels, rather simple substitutions as these are most common (ref below). This allows for
    controlling the error rate, length, and other parameters for baseline testing.

ref - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3305602/

"""
import os
import sys
import random
from pathlib import Path
import argparse
from enum import Enum
from typing import List, Tuple, Dict
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


def multi_simulate(genome_directory: Path, number_of_genomes, read_count, outfile, readlength=100, error_rate=0.0, out_dir="simread_genomes/"):
    """
    Simulates reads for multiple genomes.
    """
    # create name for  the file.
    outfile = Path(str(outfile) + f"_c{read_count}_n{number_of_genomes}_e{error_rate}.fq")

    # remove file if already exists
    if os.path.isfile(outfile):
        os.remove(outfile)
    
    # TODO: change the following to log statements.
    print(f"multi_simulate() - genome_directory: {os.listdir(genome_directory)}")
    print(f"multi_simulate() - length of genome_directory: {len(os.listdir(genome_directory))}")
    print(f"multi_simulate() - number of genomes to sample: {number_of_genomes}")

    # create subset of sampled genomes, and simulate reads from each.
    read_count_per_genome = int(read_count / number_of_genomes)
    with Experiment(number_of_genomes, genome_directory, out_dir) as exp:
        genome_directory = exp.genome_dir()
        for fasta in os.listdir(genome_directory):
            print(os.listdir(genome_directory))
            full_path = os.path.join(genome_directory, fasta)
            genome, name = parse_fasta(full_path)
            simulate_reads(genome=genome, name=name, read_count=read_count_per_genome,
                           outfile=outfile, readlength=min(len(genome), readlength), error_rate=error_rate)
    return outfile

def combine_files(fasta_paths_list: List[Path], output_file):
    """
    Combine a list of fasta files into one.
    """
    # remove file if already exists
    if os.path.isfile(output_file):
        os.remove(output_file)
    # combine fasta files.
    with open(output_file, 'w') as outfile:
        for file_path in fasta_paths_list:
            if os.path.isfile(file_path):
                with open(file_path, 'r') as infile:
                    outfile.write(infile.read())
    return Path(output_file)

class SimReadParser:
    """ class contains methods for parsing simulated read files """

    @staticmethod
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

    def get_genome_counts(simreads_path: Path):
        """
        This method obtains the genome counts from a simulated reads
        file.
        """
        regex_match = re.search(r'_n(\d+)_', simreads_path)
        if regex_match:
            return int(regex_match.group(1))
        else:
            print(f"ERROR: not able to parse genome count from {simreads_path}")
            exit(1)

    @staticmethod
    def get_error_rate(simreads_path: Path):
        """
        This method obtains the error rate from a simulated reads file
        file.
        """
        regex_match = re.search(r'_e([\d.]+).fq', simreads_path)
        if regex_match:
            return float(regex_match.group(1))
        else:
            print(f"ERROR: not able to parse (error rate): {simreads_path}")
            exit(1)