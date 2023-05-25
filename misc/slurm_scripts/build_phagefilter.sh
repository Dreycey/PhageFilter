#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=05:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=SRA-download
#SBATCH --output=SRA-download-%j.out

module purge
module load intel
module load mkl

# check input args
if [[ "$#" -lt  3 ]]; then
    echo "Usage: "$#" <genomes path> <output DB directory> <genome size>";
    exit 1
fi

echo "== Starting build  =="
# Read the input file, extract the SRA accession numbers, and parallelize the download using GNU Parallel
phage_filter build --genomes ${1} \
                   --db-path ${2} \
                   --threads 30 \
                   --false-pos-rate 0.0001 \
                   --largest-genome ${3} \
                   --kmer-size 20;

echo "== Build Finished  =="
