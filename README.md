![workflow status](https://github.com/Dreycey/PhageFilter/actions/workflows/rust.yml/badge.svg)
![Coverage](https://img.shields.io/badge/coverage-badge-placeholder-red)

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

### Verbosity level

The user can set the verbosity level. Below are different options for verbosity, which are available using the Rust [clap-verbosity-flag](https://crates.io/crates/clap-verbosity-flag) crate. If users want information about the tree being built, or other information about the particular run, use the `-vv` level of verbosity to get warnings and info.

```
-q silences output (Errors not shown.)
-v show warnings
-vv show info
-vvv show debug
-vvvv show trace
```

## Examples

1. Build the genomic sequence bloom tree (gSBT)

```bash
cargo run -- build --genomes examples/genomes/viral_genome_dir/ --db-path tree
```

2. Add genomes to an already built gSBT

```bash
cargo run -- add --genomes PATH/TO/OTHER/GENOMES/ --db-path tree/
```

3. Query examples

```bash
cargo run -- query -r examples/test_reads/ -o genomes_in_file.csv -d tree/ -f 1.0
```
