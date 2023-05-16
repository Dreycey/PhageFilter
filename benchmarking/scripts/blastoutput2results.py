"""
Description:
    This script is used to put consolidate the data
    from the output of the megablast output and 
    gather the results.

Note:
    This script is used on the output of both the megablast output and the 
    PhageFilter output using Millard lab genomes. In particular, it depends
    on the genus mapping file (we used `1Apr2023_itol_genus_annotations.txt`).
    URL: https://millardlab.org/phage-genomes-april-2023/

Usage: 
- General Usage
    python3 benchmarking/scripts/blastoutput2results.py \
        <Path to results from megablast_validate.py> \
        <Path to fasta file used for output of megablast_validate.py> \
        1Apr2023_itol_genus_annotations.txt

- EPO output:
    python3 benchmarking/scripts/blastoutput2results.py \
        benchmarking/results/RasPiData/PRJNA550482_EPO_FULL/blastresults_EPO.csv \
        benchmarking/results/RasPiData/PRJNA550482_EPO_FULL/POS_FILTERING.fa \
        1Apr2023_itol_genus_annotations.txt

       
- MONT output:
    python3 benchmarking/scripts/blastoutput2results.py \
            benchmarking/results/RasPiData/PRJNA550482_MONT_FULL/blastresults_MONT.csv \
            benchmarking/results/RasPiData/PRJNA550482_MONT_FULL/POS_FILTERING.fa \
            1Apr2023_itol_genus_annotations.txt
"""
from pathlib import Path
import sys


def parse_csv(filename: Path, delimiter=","):
    """
    Description:
        This method parses the output of the CSV into
        dictionary.

    Args:
        filename (Path): The path to the CSV file from the megablast output
    """
    read2prediction = {}
    with open(filename, 'r') as csv_file:
        csv_file.readline()
        for line in csv_file:
            if len(line) > 0 and line != None:
                read_name, evalue, top_hit = line.split(delimiter)[:3]
                if "|" in top_hit:
                    read2prediction[read_name] = top_hit.strip("\n").split("|")[
                        3].split(".")[0]
                else:
                    read2prediction[read_name] = top_hit.strip("\n")
    return read2prediction


def parse_fasta(filename: Path):
    """
    Parse fasta to get phagefilter names and sequence.
    """
    read2prediction = {}
    with open(filename, 'r') as fasta_path:
        for line in fasta_path:
            if line[0] == ">":
                read_name, genome = line.strip(">").replace(" ", "").split("|")
                sequence = fasta_path.readline()
                read2prediction[read_name] = (
                    genome.strip("\n"), sequence.strip("\n"))
    return read2prediction


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
                species2genus[species] = genus
    return species2genus


if __name__ == '__main__':
    blastcsv_path = Path(sys.argv[1])
    fullreads_path = Path(sys.argv[2])
    speciesmap_path = Path(sys.argv[3])

    # parse output
    blast_output = parse_csv(blastcsv_path)
    phagefilter_output = parse_fasta(fullreads_path)

    # get species to genus mapping
    species2genus = species2genome(speciesmap_path)

    # combine output
    combined_output = {}
    correct_species = 0
    correct_genus = 0
    total_genus = 0
    for read_name, blast_genome in blast_output.items():
        combined_output[read_name] = (
            blast_genome, phagefilter_output[read_name])

        # fix naming convention (blast used NC_029021 for Jenk)
        if blast_genome == "NC_029021":
            blast_genome = "KP719134"

        # species level metrics
        if blast_genome in phagefilter_output[read_name][0]:
            correct_species += 1

        # genus level metrics
        if blast_genome in species2genus:
            if species2genus[blast_genome] == species2genus[phagefilter_output[read_name][0].split(",")[0]]:
                correct_genus += 1

    # print statements
    print(f"total sequences: {len(blast_output.items())}")
    print(f"species agreement count: {correct_species}")
    print(f"genus agreement count: {correct_genus}")
