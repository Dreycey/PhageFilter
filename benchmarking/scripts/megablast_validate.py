"""
Description:
    This script is used to validate the mappings
    of PhageFilter using megablast. The output of
    this script is to STDOUT and piped into a csv.

Note:
    This takes a while to execute. It may be stopped 
    depending on the number of queries.

Usage:

- General Usage
    python3 PhageFilter/benchmarking/scripts/megablast_test.py <Fasta Path> 

- Example Usage
    python3 PhageFilter/benchmarking/scripts/megablast_test.py \
            PhageFilter/benchmarking/results/RasPiData/PRJNA550482_MONT_FULL/POS_FILTERING.fa \
            >> blastresults_MONT.csv
"""
import time
from Bio import SeqIO
import requests
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def parse_fasta(file):
    """
    Parse a fasta file. Again.
    """
    sequences = []
    with open(file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append((str(record.id), str(record.seq)))
    return sequences


def megablast_top_hit(sequence, database="nt", max_hits=1):
    """
    This function uses biopython to send a sequence through
    megablast for finding a match.

    Args:
        sequence (_type_): DNA sequence
        database (str, optional): Which blast db to use. Defaults to "nt".
        max_hits (int, optional): Number of hits to return. Defaults to 1.

    Returns:
        returns tuple of:
            1. top hit name (string)
            2. e-value of top hit (float)
    """
    result_handle = NCBIWWW.qblast(
        "blastn", database, sequence, megablast=True)
    blast_records = NCBIXML.parse(result_handle)

    for record in blast_records:
        if len(record.alignments) > 0:
            top_hit = record.alignments[0]
            top_hsp = top_hit.hsps[0]
            return top_hit.title, top_hsp.expect
        else:
            return None, None


if __name__ == "__main__":
    fasta_file = sys.argv[1]
    sequences = parse_fasta(fasta_file)

    print("sequence name, Evalue, hit name")
    for name, sequence in sequences:
        top_hit_name, e_value = megablast_top_hit(sequence)
        print(f"{name},{e_value},{top_hit_name}")
