# Benchmarking Pipeline

## Description

This directory contains the benchmarking pipeline used for benchmarking the Memory/Time performance, filtering performance, and classification performance. This pipeline involves several tools, so we also include *basic* installation instructions for the tools. However, it is noted that installing the tools, except for PhageFilter, is outside of the realm of our control.

## Tool installation

### PhageFilter
Follow the steps in the README.md - run the following:
```
cargo build --release
```

### Kraken2

1. Use conda to install
```
conda install -c bioconda kraken2
```
2. Build - visit https://github.com/DerrickWood/kraken2

### FastViromeExplorer
These are the steps to build FVE. Of note, the path needs to be changed within the `config.yaml`.

1. Obtain from git and compile
```
# download and install FVE
git clone https://code.vt.edu/saima5/FastViromeExplorer.git;
cd FastViromeExplorer/;
javac -d bin src/*.java;
```

2. get dependencies

* if mac
```
cd tools-mac/;
sudo cp kallisto /usr/local/bin/;
cd ../../;
```

* if linux
```
cd tools-mac/;
sudo cp kallisto /usr/local/bin/;
cd ../../;
```

* if those fail (samtools might fail on mac)
```
# download and install samtools
git clone https://github.com/samtools/htslib.git;
cd htslib;
git submodule update --init --recursive;
make;
cd ../;
git clone https://github.com/samtools/samtools.git;
cd samtools/;
make;
make install;
sudo cp samtools /usr/local/bin/;
cd ../;
```

### FACS
The FACS version used should work (this paper resulted in a change to the FACS repo!)

1. Clone from github
```
git clone https://github.com/SciLifeLab/facs.git
```

2. Build
```
cd facs;
make;
```
* then copy `facs/facs/facs` (i.e. the executable) to your bin (`cp facs/facs /usr/local/bin/`)

### BioBloomTools
The conda installation from BioBloomTools was used for benchmarking.

```
conda install -c bioconda biobloomtools
```

### CLARK
CLARK is a tool that is able to perform metagenomic classification. It only works on mac and linux.

* Download the CLARK repository + unzip
```
wget http://clark.cs.ucr.edu/Download/CLARKV1.2.6.1.tar.gz;
tar -xzvf CLARKV1.2.6.1.tar.gz
```

* Install the toolset
```
cd CLARKSCV1.2.6.1/;
/install.sh;
```

## Usage

### Benchmarking

#### performance benchmarking (when optimizing performance..)
```
python3 benchmarking/bench.py performance_testing -g examples/genomes/viral_genome_dir/ -r res_performance_benchmarking_O0.csv
```

#### running genome count benchmarking
```
python3 benchmarking/bench.py genomecount -g examples/genomes/viral_genome_dir/ -r res_genomes.csv
```

#### perform threading tests
```
python3 benchmarking/bench.py threads -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_threading.csv
```

#### read length benchmarking
```
python3 benchmarking/bench.py readlength -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -c benchmarking/config.yaml -r res_readlength.csv
```

#### running parameterization benchmarking
```
python3 benchmarking/bench.py parameterization -g examples/genomes/viral_genome_dir/ -t examples/test_reads/ -r res_parameterization.csv
```

#### classification benchmarking
```
python3 benchmarking/bench.py relative_performance -g examples/genomes/viral_genome_dir/ -t examples/test_reads/ -c benchmarking/config.yaml -r res_relative_performance.csv
```

#### filter performance benchmarking
```
python3 benchmarking/bench.py filter_performance -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -t examples/test_reads/ -c benchmarking/config.yaml -r res_filter_performance.csv
```

#### filter memory benchmarking
```
python3 benchmarking/bench.py filter_memory -g examples/genomes/viral_genome_dir/ -n examples/genomes/bacteria_genome_dir/ -t examples/test_reads/ -c benchmarking/config.yaml -r res_filter_memory.csv
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
