"""
This contains the python wrapper for FACS

GitHub: https://github.com/SciLifeLab/facs

Design Pattern: Template pattern.
"""
from bench.tools.tool_template import ToolOp
from pathlib import Path
from collections import Counter
from typing import List, Tuple, Dict
import sys
import os




class FACS(ToolOp):

    def __init__(self, kmer_size: int, filter_thresh: float, threads=1):
        """_summary_

        Args:
            kmer_size (int): size of the kmer to use for creating the bloom filters
            threads (int, optional): number of threads to use while running.
        """
        self.k = kmer_size
        self.theta = filter_thresh
        self.threads = threads
        self.db_path = None
        self.suffix = ".bloom"

    def parse_output(self, output_path: Path, genomes_path: Path = None, filter_reads=False) -> Dict[str, int]:
        """_summary_
        parses an output file/directory (depends on tool)
        returns a dictionary of the output of FACS.

        Args:
            output_path (Path): Path where the output of FACS
                                will be stored.

        Returns:
            List[]
        """
        read_set = []
        # parse output file name
        test_file_name = output_path.strip("FACS_")
        full_out_path = f"{output_path}{test_file_name}_{self.db_path[:-len(self.suffix)]}_contam.fastq"

        # create counter.
        read_counter = Counter()

        # parse output file.
        with open(full_out_path, "r") as opened_file:
            line = opened_file.readline()
            while line:
                if line[0] == "@":
                    genome_name = "_".join(line.strip("@").strip("\n").split("_")[0:2])
                    read_counter[genome_name] += 1
                line = opened_file.readline()

        return read_counter

    @staticmethod
    def combine_files_into_one(genomes_path: Path, output_file):
        """
        Combine multiple fasta files into one file for
        using with FACS build command.
        """
        # remove file if already exists
        if os.path.isfile(output_file):
            os.remove(output_file)
        # combine
        with open(output_file, 'w') as outfile:
            for filename in os.listdir(genomes_path):
                file_path = os.path.join(genomes_path, filename)
                if os.path.isfile(file_path):
                    with open(file_path, 'r') as infile:
                        outfile.write(infile.read())
                        outfile.write('\n')

    def build(self, db_path: Path, genomes_path: Path) -> List[List[str]]:
        """_summary_
        Build the tools database.

        Args:
            db_path (Path): path to the database the tool needs.
            genomes_path (Path): path to the genomes the db will use for building.

        Returns:
            List[List[str]]: a nested list of command line arguments.
        """
        # make bloom filter
        self.db_path = db_path+self.suffix
        tmp_file_name = "TMP_DELETE.fa"

        # combine genomes into one file
        FACS.combine_files_into_one(genomes_path, tmp_file_name)

        # make command
        build_cmd = ["facs", "build"]
        build_cmd += ["-r", f"{tmp_file_name}"]
        build_cmd += ["-o", f"{self.db_path}"]
        build_cmd += ["-k", f"{self.k}"]
        build_cmd += ["-e", f"{0.001}"] # false positive rate

        # remove tmp file
        rm_genome_build = [f"rm", f"{tmp_file_name}"]

        return [build_cmd, rm_genome_build]

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
            print("Must first build (FACS)")
            exit()
        run_cmd = ["facs", "remove"]
        run_cmd += ["-r", f"{self.db_path}"]
        run_cmd += ["-q", f"{fasta_file}"]
        run_cmd += ["-t", f"{self.theta}"]
        run_cmd += ["-o", f"{output_path}"]
    
        return [run_cmd]

if __name__ == "__main__":
    facs = FACS(kmer_size=25, filter_thresh=1.0)

    # build BBT
    facs_build_cmd = facs.build(db_path=sys.argv[1], genomes_path=sys.argv[2])
    print(" ".join(facs_build_cmd[0]))
    print(" ".join(facs_build_cmd[1]))
    print(" ".join(facs_build_cmd[2]))

    # run command
    facs_run_cmd = facs.run(fasta_file="INPUT", output_path="OUT", filter_reads=True)
    print(" ".join(facs_run_cmd[0]))

    # # parsing
    parsed_out = facs.parse_output(output_path="OUT")
    print(parsed_out)