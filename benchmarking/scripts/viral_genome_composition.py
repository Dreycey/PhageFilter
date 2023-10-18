"""
Script to evaluate viral genome composition. This script is used
as a basic way to get information on the composition of viruses
in a directory of genomes.

The output from this script is an ordered csv list of genomes
ordered by the count.

Usage:
    python benchmarking/scripts/viral_genome_composition.py examples/genomes/viral_genome_dir/

To transfer viruses consisting of the top 15 genomes:
    python benchmarking/scripts/viral_genome_composition.py examples/genomes/viral_genome_dir/ subselected/ 15
"""
from pathlib import Path
import sys
from typing import List, Dict, Set
import os
from collections import Counter
import shutil



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

def get_viral_names_from_directory(genome_directory_path : Path) -> Dict[Path, List[str]]:
    """
    Obtain the genomic names from all genomic files in a directory.

    Returns:
        viral_names (Dict[Path, List[str]]) : mapping from file path to parsed viral names.
                                              one name if one entry (in the case of complete genomes.)
    """
    path2viral_names = {}
    for genome_file in os.listdir(genome_directory_path):
        genome_file_path = Path(os.path.join(genome_directory_path, genome_file))
        preparsed_viral_names = get_fasta_names(genome_file_path)
        path2viral_names[genome_file_path] = get_virus_name(preparsed_viral_names)
    return path2viral_names

def genus_count(path2viral_names : Dict[Path, List[str]]) -> Counter:
    """
    Get counts of each genus from a list of viral names. This
    method assumes the first string is the genus.
    """
    genus_count = Counter()
    for genome_path, viral_names in path2viral_names.items():
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

def get_viralgenus2paths(path2viral_names : Dict[Path, List[str]], 
                         genus_name_counter : Counter, 
                         top_n : int) -> Dict[Path, str]:
    """
    Get the paths for the top N genomes.
    """
    # get top genus names.
    top_genus_names = set()
    counter = 0
    for name, count in genus_name_counter.most_common():
        if counter >= top_n: break
        top_genus_names.add(name)
        counter += 1

    # get paths to top n genomes
    viralgenus2paths = {}
    for genome_path, viral_names in path2viral_names.items():
        for viral_name in viral_names:
            genus_name = viral_name.split(" ")[0]
            if genus_name in top_genus_names:
                viralgenus2paths[genome_path] = genus_name

    return viralgenus2paths

if __name__ == "__main__":
    genome_directory = Path(sys.argv[1])
    new_directory = None
    if len(sys.argv) > 2:
        new_directory = Path(sys.argv[2])
        number_of_hosts = int(sys.argv[3])

    # parse genomes and count genus occurances
    path2viral_names = get_viral_names_from_directory(genome_directory)
    genus_name_count = genus_count(path2viral_names)
    print_names_as_csv(genus_name_count)

    if new_directory != None:
        # get directory of genomes with top N hosts.
        paths2genus = get_viralgenus2paths(path2viral_names, genus_name_count, number_of_hosts)

        # create new directory if not exists.
        if not os.path.exists(new_directory):
            os.mkdir(new_directory)

        # copy files over
        with open(new_directory.name + ".csv", "w") as outfile_csv:
            for viral_path, genus in paths2genus.items():
                shutil.copyfile(viral_path, os.path.join(new_directory, viral_path.name))
                outfile_csv.write(f"{viral_path.name},{genus}\n")
