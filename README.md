![software workflow status](https://github.com/Dreycey/PhageFilter/actions/workflows/rust.yml/badge.svg)
![Coverage](https://img.shields.io/badge/coverage-30.65%25-brightgreen.svg?style=flat-square)
![benchmarking workflow status](https://github.com/Dreycey/PhageFilter/actions/workflows/benchmarking_tests.yaml/badge.svg)
![LICENSE](https://img.shields.io/badge/license-GPL--3.0-brightgreen)

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

## Installation
PhageFilter is written in Rust and uses the Cargo package manager for installing all dependencies. After [Installing Rust](https://www.rust-lang.org/tools/install), running the following command from the root directory will install all dependencies as well as compile the software.

```
cargo build --release
```

For easier calling, it is recomended the compiled executable, `target/release/phage_filter`, be stored in a `bin/` directory and added to the users environmental path.

## Examples
These commands assume the user has moved the executable to the enironmental path. If not, call the binary (i.e. `target/release/phage_filter`) directly instead of `phage_filter`.


1. **Build** the genomic sequence bloom tree (gSBT)

```bash
target/release/phage_filter build --genomes examples/genomes/viral_genome_dir/ --db-path tree
```

2. **Query** reads

```bash
target/release/phage_filter query -r examples/test_reads/ -o genomes_in_file.csv -d tree/ -f 1.0
```

3. **Add** *new* genomes to an already built gSBT

```bash
target/release/phage_filter add --genomes PATH/TO/OTHER/GENOMES/ --db-path tree/
```
