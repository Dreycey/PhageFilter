"""
Description:
    This script is used to generate exact substrings of
    a given genome (in fasta format).

Usage:
    python <fasta genome path> <number of reads/substrings> <output file>

Example:
    python src/simualte_reads.py genomes/GCF_000912255.1_ViralProj226726_genomic.fna 1000 simulated_reads.fa
"""
import sys
import random 



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
            line = f.readline()
            count += 1
    return genome
            
def simulate_reads(genome, read_count, outfile, readlength=100):
    """
    This does not actually simulate reads, as
    it creates a fasta of perfect substrings.
    """
    print(f"Simulating reads to: {outfile}")
    reads_added = 0
    with open(outfile, "a+") as out:
        while reads_added < read_count:
            start = random.randint(0, len(genome)-1)
            stop = start + readlength
            read = genome[start:stop]
            out.write(f">read_{reads_added}\n{read}\n")
            reads_added += 1

if __name__ == "__main__":
    # parse arguments
    genome = parse_fasta(sys.argv[1])
    read_count = int(sys.argv[2])
    outfile = sys.argv[3]
    # simulate reads
    simulate_reads(genome, read_count, outfile)