#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=10:00:00
#SBATCH --partition=amilan
#SBATCH --job-name=PhageFilter-QUERY
#SBATCH --output=PhageFilter-query-%j.out

module purge

module load intel
module load mkl
module load openmpi
module load gnu_parallel

# check input args
if [[ "$#" -lt  2 ]]; then
    echo "Usage: "$#" <reads path> <DB directory> <output directory>";
    exit 1
fi

echo "== Starting query  =="

phage_filter query --reads ${1} \
                   --db-path ${2} \
                   --out ${3}/$(basename ${1})_OUTPUT/ \
                   --threads 50 \
                   --filter-threshold 0.7 \
                   --block-size-reads 100000 \
                   --cache-size 1000 \
                   --pos-filter

echo "== Query Finished  =="
