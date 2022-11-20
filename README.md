![workflow status](https://github.com/Dreycey/PhageFilter/actions/workflows/rust.yml/badge.svg)
![phage filter logo](misc/PhageFilterLogo.png)

# PhageFilter

PhageFilter uses a Sequence Bloom Tree (SBT) to filter bacteriophage reads from metagenomic files.

## Usage

```
A fast, simple, and efficient method for taxonomic classification.

Usage: phage_filter [OPTIONS] --genomes <VALUE> --reads <VALUE> --out <VALUE>

Options:
  -g, --genomes <VALUE>    Path to genomes file or directory. (Fasta)
  -r, --reads <VALUE>      Path to read file or directory of reads. (Fasta or Fastq, or dirs with both)
  -o, --out <VALUE>        Path to output file. (Fasta)
  -t, --threads <VALUE>    Number of threads to use for read matching [default: 4]
  -k, --kmer_size <VALUE>  Size of the kmer to use; use with caution! [default: 20]
  -q, --threshold <VALUE>  Filtering theshold (Number of kmers needed to pass) [default: 1.0]
  -h, --help               Print help information
  -V, --version            Print version information
```

## Examples

1. Build the sequence bloom tree

```bash
cargo run -- build --genomes examples/genomes/viral_genome_dir/ --db-path tree -t 6
```

2. Query examples

```bash
cargo run -- query -r examples/test_reads/simulated_reads.fa -o genomes_in_file.csv -d tree -c 1.0 -t 6
```
