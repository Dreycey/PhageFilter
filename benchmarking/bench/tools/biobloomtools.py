"""
This contains the python wrapper for BioBloomTools
"""
from bench.tools.tool_template import ToolOp
from pathlib import Path
from typing import List, Tuple, Dict
import subprocess
import sys
import os
from collections import Counter




class BioBloomTools(ToolOp):

    def __init__(self, kmer_size: int, threads=1):
        """_summary_

        Args:
            kmer_size (int): size of the kmer to use for creating the bloom filters
            threads (int, optional): number of threads to use while running.
        """
        self.k = kmer_size
        self.threads = threads
        self.db_path = None

    def parse_output(self, output_path: Path, genomes_path: Path = None, filter_reads=False) -> Dict[str, int]:
        """_summary_
        parses an output file/directory (depends on tool)
        returns a dictionary of the output of PhageFilter.

        Args:
            output_path (Path): Path where the output of PhageFilter
                                will be stored.

        Returns:
            Dict[str, int]: A map from NCBI ID to read count
        """
        if filter_reads:
            read_counter = Counter()
            with open(f"{output_path}.fa", "r") as opened_file:
                line = opened_file.readline()
                while line:
                    if line[0] == ">":
                        genome_name = "_".join(line.strip(">").split("\t")[0].split("_")[:-1])
                        read_counter[genome_name] += 1
                    line = opened_file.readline()
            return read_counter
        else:
            name2counts = {}
            with open(output_path + "_summary.tsv") as out_file:
                out_file.readline() # skip first
                line = out_file.readline()
                count = 0
                while line:
                    name, count = line.strip("\n").split("\t")[:2]
                    name2counts[name] = int(count)
                    line = out_file.readline()

            # remove non genome columns
            del name2counts['repeat']
            del name2counts['noMatch']
            del name2counts['multiMatch']

            return name2counts

    def build(self, db_path: Path, genomes_path: Path) -> List[List[str]]:
        """_summary_
        Build the tools database.

        Args:
            db_path (Path): path to the database the tool needs.
            genomes_path (Path): path to the genomes the db will use for building.

        Returns:
            List[List[str]]: a nested list of command line arguments.
        """
        build_cmd = ["biobloommimaker"]
        build_cmd += ["--file_prefix", f"{db_path}"]
        build_cmd += ["--hash_num", f"{100}"] 
        build_cmd += ["--kmer_size", f"{self.k}"]
        build_cmd += ["--threads", f"{self.threads}"]
        for filename in os.listdir(genomes_path):
            filepath = os.path.join(genomes_path, filename)
            build_cmd += [f"{filepath}"]

        self.db_path = db_path+".bf"

        return [build_cmd]

    def run(self, fasta_file: Path, output_path: Path, filter_reads=False):
        """_summary_
        run tool, based on input arguments, it outputs a CMD-line array.

        Args:
            fasta_file (Path): Path to simualted reads.
            output_path (Path): Desired path for the output file.

        Returns:
            N/A.
        """
        if not self.db_path:
            print("Must first build (PhageFilter)")
            exit()
        run_cmd = ["biobloommicategorizer"]
        run_cmd += ["--filter", f"{self.db_path}"]
        run_cmd += ["--multi", f"{2}"] # how many reads per genome
        run_cmd += ["--prefix", f"{output_path}"] # prefix
        run_cmd += ["--min_FPR", f"{2}"] # Minimum -10*log(FPR) threshold for a match
        run_cmd += ["--threads", f"{self.threads}"]
        run_cmd += [f"{fasta_file}"]
        if filter_reads:
            run_cmd += ["--hitOnly", "--fa"]
            output_file = f"{output_path}.fa"
            with open(output_file, "w") as out:
                subprocess.run(run_cmd, stdout=out, check=True)
    
        return [run_cmd]

if __name__ == "__main__":
    BBT = BioBloomTools(kmer_size=25)
    # build BBT
    bbt_build_cmd = BBT.build(db_path=sys.argv[1], genomes_path=sys.argv[2])
    print(" ".join(bbt_build_cmd[0]))
    # run command
    bbt_run_cmd = BBT.run(fasta_file="INPUT", output_path="OUT", filter_reads=True)
    print(" ".join(bbt_run_cmd[0]))
    # parsing
    parsed_out = BBT.parse_output(output_path="OUT")
    print(parsed_out)