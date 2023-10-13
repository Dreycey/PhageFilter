"""
Script to evaluate viral genome composition. This script is used
as a basic way to get information on the composition of viruses
in a directory of genomes.

The output from this script is an ordered csv list of genomes
ordered by the count.

Usage:
    python benchmarking/scripts/viral_genome_composition.py examples/genomes/viral_genome_dir/
"""
from pathlib import Path
import sys
from typing import List
import os
from collections import Counter




def get_fasta_names(fasta_path : Path):
    """
    Grab the names inside a fasta file.
    """
    viral_names = []
    with open(fasta_path, "r") as fasta_file:
        fasta_line = fasta_file.readline()
        while fasta_line:
            if fasta_line.startswith(">"):
                viral_names.append(fasta_line.strip(">").strip("\n"))
            fasta_line = fasta_file.readline()
    return viral_names

def get_virus_name(fasta_string_list : List[str]):
    """
    Parse out the viral name from a fasta string.

    NOTE:
        This logic is particular the NCBI strings. This
        can be changed as needed/desired.
    """
    for index, fasta_string in enumerate(fasta_string_list):
        fasta_string_list[index] = " ".join(fasta_string.split(",")[0].split(" ")[1:])
    return fasta_string_list

def get_viral_names_from_directory(genome_directory_path : Path):
    """
    Obtain the genomic names from all genomic files in a directory.
    """
    viral_names = []
    for genome_file in os.listdir(genome_directory_path):
        preparsed_viral_names = get_fasta_names(os.path.join(genome_directory_path, genome_file))
        viral_names += get_virus_name(preparsed_viral_names)
    return viral_names

def genus_count(viral_names : List[str]) -> Counter:
    """
    Get counts of each genus from a list of viral names. This
    method assumes the first string is the genus.
    """
    genus_count = Counter()
    for viral_name in viral_names:
        genus_name = viral_name.split(" ")[0]
        genus_count[genus_name] += 1
    return genus_count

def print_names_as_csv(genus_counts : Counter):
    """
    Print the genus names and counts as a csv to stdout.
    """
    for name, count in genus_counts.most_common():
        print(f"{name},{count}")

if __name__ == "__main__":
    genome_directory = Path(sys.argv[1])

    # parse genomes and count genus occurances
    viral_names = get_viral_names_from_directory(genome_directory)
    genus_name_count = genus_count(viral_names)
    print_names_as_csv(genus_name_count)