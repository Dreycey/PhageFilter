"""
Description:
    This is used to transform the output of PhageFilter
    Species-level classifications into genus level 
    classifications for Millard lab based genomes.

Note:
    This script is used on the output of the PhageFilter output using Millard 
    lab genomes. In particular, it depends on the genus mapping file 
    (we used `1Apr2023_itol_genus_annotations.txt`).
    URL: https://millardlab.org/phage-genomes-april-2023/

- General Usage
    python3 benchmarking/scripts/genus_abundances.py \
            <path to PhageFilter classification csv> \
            1Apr2023_itol_genus_annotations.txt 

- Example Usage
    python3 benchmarking/scripts/genus_abundances.py \
            benchmarking/results/RasPiData/PRJNA550482_EPO_FULL/CLASSIFICATION.csv \
            1Apr2023_itol_genus_annotations.txt \
            >> GENUS_CLASSIFICATION.csv
    
"""
from pathlib import Path
import sys


def get_species2readcount(classification_path: Path):
    """
    given the filepath to the classifications for PhageFilter
    this function returns the species abundance mappings.
    """
    species2readcount = {}
    with open(classification_path, 'r') as f:
        f.readline()  # skip first line
        for line in f:
            species, read_count = line.split(",")
            species2readcount[species] = int(read_count.strip("\n"))
    return species2readcount


def species2genome(filename: Path):
    """
    parses the Millard lab species-to-genus tsv file
    returning a dictionary map.
    """
    species2genus = {}
    with open(filename, 'r') as speciesmap_file:
        # skip lines not containing data
        for line in speciesmap_file:
            if line.strip() == "DATA":
                break
        # get species-to-genus mapping
        for line in speciesmap_file:
            if len(line.strip()) > 0:  # skip empty lines
                species, _, genus = line.split("\t")
                species2genus[species] = genus.strip("\n")
    return species2genus


def get_genus_readcount(species2readcount, species2genus):
    """
    calculates the genus level read counts
    given the species read counts.
    """
    genus2readcount = {}
    for species, readcount in species2readcount.items():
        genus = species2genus[species]
        if genus in genus2readcount:
            genus2readcount[genus] += readcount
        else:
            genus2readcount[genus] = readcount
    return genus2readcount


if __name__ == "__main__":
    # arguments
    phagefilter_classification_path = Path(sys.argv[1])
    speciesmap_path = Path(sys.argv[2])

    # get species to genus mapping
    species2genus = species2genome(speciesmap_path)

    # get species-to-abundance mapping
    species2readcount = get_species2readcount(phagefilter_classification_path)

    # get read counts for each genus
    genus2readcount = get_genus_readcount(species2readcount, species2genus)

    # print to stdout
    print("genus,readcount")
    for genus, readcount in genus2readcount.items():
        print(f"{genus},{readcount}")
