"""
Description:
    This script is used to generate exact substrings of
    a given genome (in fasta format).

Usage:
    python benchmarking/simulate_reads.py <fasta genome path> <number of reads/substrings> <output file> <read length>

Example:
    python benchmarking/simulate_reads.py genomes/GCF_000912255.1_ViralProj226726_genomic.fna 1000 simulated_reads.fa 100

For many genomes:
```bash
    for genome in examples/genomes/viral_genome_dir/*; \
    do python benchmarking/simulate_reads.py $file 1000 \
    simulated_reads.fa 1000; done
```
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
            else:
                name = line.strip(">").strip("\n").split(" ")[0]
            line = f.readline()
            count += 1
    return genome, name
            
def simulate_reads(genome, name, read_count, outfile, readlength=100, error_rate=0.0):
    """
    This does not actually simulate reads, as
    it creates a fasta of perfect substrings.
    """
    print(f"Simulating reads to: {outfile}")
    reads_added = 0
    with open(outfile, "a+") as out:
        while reads_added < read_count:
            start = random.randint(0, len(genome)-read_length)
            stop = start + readlength
            read = list(genome[start:stop])
            if random.uniform(0, 1) < error_rate:
                base2change = random.randint(0, len(read)-1)
                read[base2change] = random.choice(["A", "C", "T", "G"])
            read = "".join(read)
            out.write(f">{name}_{reads_added}\n{read}\n")
            reads_added += 1

if __name__ == "__main__":
    # parse arguments
    genome, name = parse_fasta(sys.argv[1])
    read_count = int(sys.argv[2])
    outfile = sys.argv[3]
    read_length = int(sys.argv[4])
    error_rate = float(sys.argv[5])
    # simulate reads
    simulate_reads(genome, name, read_count, outfile, readlength=read_length, error_rate=error_rate)