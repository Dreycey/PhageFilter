#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --time=00:30:00
#SBATCH --partition=atesting
#SBATCH --job-name=SRA-download
#SBATCH --output=SRA-download-%j.out

module purge

module load intel
module load mkl
module load gnu_parallel

# Check if the input file is provided
if [ -z "$1" ]; then
  echo "Usage: $0 input.csv prefetch/fastq <download path>"
  exit 1
fi

# Check if the second argument is valid
if [[ "$2" != "prefetch" && "$2" != "fastq" ]]; then
  echo "Invalid second argument: must be either 'prefetch' or 'fastq'"
  exit 1
fi

# Check if the second argument is valid
if [[ -z "$1" ]]; then
  echo "Must have storage location specified: "
  echo "Usage: $0 input.csv prefetch/fastq <download path>"
  exit 1
fi

echo "==> Starting downloads of SRA data  <=="
# Read the input file and extract the SRA accession numbers
TMP_DIR="/home/dral3008/project-directory/tmp/";
FILENAME=$(basename $1)
DOWNLOAD_PATH=$3
if [ "$2" == "prefetch" ]; then
  tail -n +1 "$1" | cut -d',' -f1 | parallel --tmpdir $TMP_DIR "echo 'Prefetching {}'; prefetch -O ${DOWNLOAD_PATH} {}"
else
  tail -n +1 "$1" | cut -d',' -f1 | parallel --tmpdir $TMP_DIR "echo 'Converting {} to fastq'; fastq-dump --clip -O ${DOWNLOAD_PATH}${FILENAME%%.*}_READS {}"
fi
echo "==> End of Downloads  <=="
