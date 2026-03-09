# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Does

PhageFilter is a metagenomic filtering tool that uses a **genomic Sequence Bloom Tree (gSBT)** to identify and filter bacteriophage reads from metagenomic sequencing files. Each leaf node in the tree corresponds to one genome (stored as a `.bf` Bloom filter file on disk); internal nodes store the union of their children's filters to enable efficient tree pruning during queries.

## Build & Test Commands

```bash
# Build (release — required before running benchmarks)
cargo build --release

# Build (debug)
cargo build

# Run all tests
cargo test

# Run a single test by name
cargo test test_get_kmers
cargo test test_distance

# Run tests in a specific module
cargo test --test-target bloom_filter
```

The compiled binary is at `target/release/phage_filter`.

## Usage (three subcommands)

```bash
# 1. Build a gSBT from genome FASTA files
target/release/phage_filter build --genomes examples/genomes/viral_genome_dir/ --db-path tree/

# 2. Query reads against the tree; outputs CLASSIFICATION.csv and optional filtered FASTA files
target/release/phage_filter query -r examples/test_reads/ -o output/ -d tree/ -f 1.0

# 3. Add new genomes to an existing tree
target/release/phage_filter add --genomes PATH/TO/NEW/GENOMES/ --db-path tree/
```

Key query flags: `--filter-threshold` (0–1, fraction of k-mers required to match), `--search-depth` (prune tree to this depth before querying), `--pos-filter` / `--neg-filter` (write matching / non-matching reads to FASTA files).

## Architecture

```
src/
  main.rs         — CLI (clap), subcommand dispatch, rayon thread pool setup
  bloom_filter.rs — BloomFilter struct with bitvec bits, insert/contains, union/intersect,
                    Hamming distance, Drop-based auto-save to disk
  bloom_tree.rs   — BloomTree (gSBT) and BloomNode; greedy nearest-neighbor insertion;
                    save/load via bincode; prune_tree for depth-limited queries
  cache.rs        — BFLruCache (LRU cache of Arc<RwLock<BloomFilter>>); BloomFilterCache trait
                    for dependency injection; evicted BFs are reloaded from disk on demand
  file_parser.rs  — DNASequence, RecordTypes (FASTA/FASTQ), ReadQueue (block-based streaming);
                    k-mer extraction using canonical (lexicographically smallest) strand
  query.rs        — query_batch / _query_batch (recursive SBT traversal using rayon par_iter);
                    save_leaf_counts writes CLASSIFICATION.csv
  result_map.rs   — ResultMap: read-id → set-of-genome-ids mapping for filter output
```

### Key design decisions

- **Bloom filters are persisted per-node** as `.bf` files in the `--db-path` directory. The serialized tree structure (`tree.bin`) only stores node metadata and file paths, not the BF bits.
- **The LRU cache** (`BFLruCache`) keeps a bounded set of BFs in memory. On cache miss, the BF is deserialized from disk. On drop of a modified BF, it is automatically re-serialized (via `Drop` impl).
- **Canonical k-mers**: `file_parser::get_lex_less` picks the lexicographically smaller of a k-mer and its reverse complement, so both strands map to the same BF entry.
- **Rayon** is used throughout for parallel k-mer extraction, parallel read querying within a block, and parallel filter output writing.
- **Greedy SBT insertion**: a new leaf is placed adjacent to the existing node with the smallest Hamming distance between BF bit vectors.

## Benchmarking

The `benchmarking/` directory contains a Python benchmarking harness comparing PhageFilter to Kraken2, FastViromeExplorer, BioBloomTools, CLARK, and FACS. Run with:

```bash
pip install -r benchmarking/requirements.txt
python3 benchmarking/bench.py <subcommand> -g <viral_genomes_dir> -n <bacteria_genomes_dir> -c benchmarking/config.yaml -r results.csv
```

Subcommands: `performance_testing`, `genomecount`, `threads`, `readlength`, `parameterization`, `relative_performance`, `filter_performance`.

## Testing Notes

- Tests use `tmp_testing/` as a scratch directory (created/destroyed per test in `bloom_tree.rs`).
- `proptest` is used for property-based tests in `cache.rs`.
- Many tests in `query.rs` and `bloom_tree.rs` are commented out — they are integration tests that need to be updated after recent refactors.
