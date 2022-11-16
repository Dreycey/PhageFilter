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

- Using on simulated reads with 7 threads.

```bash
cargo run -- -g examples/genomes/viral_genome_dir/ -r examples/test_reads/simulated_reads.fa -t 6 --out genomes_in_file.txt -k 15 -q 1
```
