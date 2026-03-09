![software workflow status](https://github.com/Dreycey/PhageFilter/actions/workflows/rust.yml/badge.svg)
![Coverage](https://img.shields.io/badge/coverage-59.88%25-brightgreen.svg?style=flat-square)
![benchmarking workflow status](https://github.com/Dreycey/PhageFilter/actions/workflows/benchmarking_tests.yaml/badge.svg)
![LICENSE](https://img.shields.io/badge/license-GPL--3.0-brightgreen)

![phage filter logo](misc/PhageFilterLogo.png)

# PhageFilter

PhageFilter uses a Sequence Bloom Tree (SBT) to filter bacteriophage reads from metagenomic files.

## Usage
PhageFilter has three primary commands: (1) Build, (2) Query, and (3) Add. Each of these functions contains it's own set of required arguments that can be listed using the `--help` flag.

```
Usage: phage_filter [OPTIONS] <COMMAND>

Commands:
  build  Builds the BloomTree
  add    Adds genomes to an already built BloomFilter
  query  Queries a set of reads. (ran after building the bloom tree)
```

### Supported File Formats

PhageFilter accepts FASTA and FASTQ files, including **gzip-compressed** variants.

| Format | Recognized Extensions |
|--------|----------------------|
| FASTA  | `.fa`, `.fasta`, `.fna`, `.fsa`, `.fas` |
| FASTQ  | `.fq`, `.fastq` |
| Gzip   | Any of the above with `.gz` / `.gzip` suffix (e.g. `.fq.gz`, `.fasta.gz`) |

**Format auto-detection:** By default (`--format auto`), PhageFilter inspects the first bytes of each file to determine the format:
- `>` → FASTA
- `@` → FASTQ
- Gzip magic bytes (`1f 8b`) → decompresses and re-inspects

If content sniffing is inconclusive, the file extension is used as a fallback (FASTQ for `.fq`/`.fastq`, FASTA for everything else).

**Manual override:** Use `--format fasta` or `--format fastq` (`-F`) to force a specific parser, bypassing auto-detection.

When a directory is provided as input, PhageFilter scans for files with the extensions listed above (including `.gz` compound extensions) and silently skips unrecognized files.

**Filtering output format:** When using `--pos-filter` or `--neg-filter`, the output files (`POS_FILTERING` / `NEG_FILTERING`) inherit the input format. FASTQ input produces FASTQ output (with quality scores preserved); FASTA input produces FASTA output.

### Query Output

The `query` command produces the following output files in the directory specified by `-o`:

**`CLASSIFICATION.csv`** — A summary of read counts per genome. Each row is `genome_id,count` where `genome_id` is the sequence ID (from the FASTA header) of a reference genome in the gSBT, and `count` is the number of query reads that mapped to it. Only genomes with at least one mapped read are listed.

**`POS_FILTERING.fa` / `POS_FILTERING.fq`** (when `--pos-filter` is used) — Reads that matched at least one genome in the gSBT. The header of each output read is annotated with the matched genome(s):

```
>read_id |genome_A,genome_B
ACGTACGT...
```
or for FASTQ input:
```
@read_id |genome_A,genome_B
ACGTACGT...
+
IIIIIIII...
```

The `|`-delimited suffix lists all genome IDs whose bloom filter the read passed at the given threshold (`-f`). A read may match multiple genomes. The genome IDs correspond to the sequence IDs from the reference FASTA files used during `build`.

**`NEG_FILTERING.fa` / `NEG_FILTERING.fq`** (when `--neg-filter` is used) — Reads that did **not** match any genome. These retain their original header with no annotation.

**Matching semantics:** A read is considered to match a genome when at least the fraction of its k-mers specified by `--filter-threshold` (`-f`, default 1.0) are present in that genome's bloom filter. The match is binary (pass/fail at the threshold); no per-genome score or ranking is currently reported.

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
