# Benchmarking Pipeline

## Description

This directory contains the benchmarking pipeline used for benchmarking the Memory/Time performance, filtering performance, and classification performance.

## Usage

### Benchmarking

1. performance benchmarking (when optimizing performance..)

```
python benchmarking/bench.py performance_testing -g examples/genomes/viral_genome_dir/ -r res_performance_benchmarking_O0.csv
```

#### 2. running genome count benchmarking

```
python benchmarking/bench.py genomecount -g examples/genomes/viral_genome_dir/ -r res_genomes.csv
```

#### 3. running parameterization benchmarking

```
python benchmarking/bench.py parameterization -g examples/genomes/viral_genome_dir/ -t examples/test_reads/ -r res_parameterization.csv
```

#### 4. relative performance benchmarking

```
python benchmarking/bench.py relative_performance -g examples/genomes/viral_genome_dir/ -r res_relative_performance.csv -t examples/test_reads/ -c benchmarking/config.yaml
```

### Simulating reads (really genome subfragments - not _true read_ simulation)

#### 1. Single genome simulation.

-   Usage

```
python3 benchmarking/simulate_reads.py single_genome -g <genome path - fasta> --c <number of total reads> -l <read length> -n <number of genomes> -e <error rate> -o <outfile name>
```

-   Example

```
python3 benchmarking/simulate_reads.py single_genome -g examples/genomes/viral_genome_dir/GCF_000912115.1_ViralProj219123_genomic.fna -c 100 -o sim_reads.fq -l 100 -e 0.1
```

#### 2. Multigenome simulations.

-   Usage

```
python3 benchmarking/simulate_reads.py multi_genome -g <genome directory> -c <number of total reads> -l <read length> -n <number of genomes> -e <error rate> -o <outfile prefix>
```

-   Example

```
python3 benchmarking/simulate_reads.py multi_genome -g examples/genomes/viral_genome_dir/ -c 100 -o sim_reads.fq -l 100 -e 0.1
```
