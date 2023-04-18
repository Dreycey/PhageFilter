"""
This contains the python wrapper for PhageFilter
"""
from bench.tools.tool_template import ToolOp
from pathlib import Path
from typing import List, Tuple, Dict




class PhageFilter(ToolOp):

    def __init__(self, kmer_size: int, filter_thresh: float, threads=1):
        """_summary_

        Args:
            kmer_size (int): _description_
            filter_thresh (float): _description_
            threads (int, optional): _description_. Defaults to 4.
        """
        self.k = kmer_size
        self.theta = filter_thresh
        self.threads = threads
        self.db_path = None

    def parse_output(self, output_path: Path, genomes_path: Path = None, cuttoff=0.025) -> Dict[str, int]:
        """_summary_
        parses an output file/directory (depends on tool)
        returns a dictionary of the output of PhageFilter.

        Args:
            output_path (Path): Path where the output of PhageFilter
                                will be stored.

        Returns:
            Dict[str, int]: A map from NCBI ID to read count
        """
        name2counts = {}
        with open(output_path + "/CLASSIFICATION.csv") as out_file:
            line = out_file.readline()
            count = 0
            while line:
                name, count = line.strip("\n").split(",")
                name2counts[name] = int(count)
                line = out_file.readline()

        # filter based on read count threshold
        # total_reads_classified = sum(name2counts.values())
        # name2counts = {
        #     k: v for k, v in name2counts.items() if v > cuttoff*total_reads_classified}

        return name2counts

    def build(self, db_path: Path, genomes_path: Path, cache_size=100) -> List[List[str]]:
        """_summary_
        Build the tools database.

        Args:
            db_path (Path): path to the database the tool needs.
            genomes_path (Path): path to the genomes the db will use for building.

        Returns:
            List[List[str]]: a nested list of command line arguments.
        """
        build_cmd = ["./target/release/phage_filter", "build"]
        build_cmd += ["--genomes", f"{genomes_path}"]
        build_cmd += ["--db-path", f"{db_path}"]
        build_cmd += ["--kmer-size", f"{self.k}"]
        build_cmd += ["--cache-size", f"{cache_size}"]
        build_cmd += ["--false-pos-rate", f"{0.0001}"]
        build_cmd += ["--largest-genome", f"{100000}"]
        build_cmd += ["--threads", f"{self.threads}"]
        self.db_path = db_path

        return [build_cmd]

    def run(self, fasta_file: Path, output_path: Path, cache_size=10, filter_reads=False):
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
        run_cmd = ["./target/release/phage_filter", "query"]
        run_cmd += ["--reads", f"{fasta_file}"]
        run_cmd += ["--db-path", f"{self.db_path}"]
        run_cmd += ["--filter-threshold", f"{self.theta}"]
        run_cmd += ["--cache-size", f"{cache_size}"]
        run_cmd += ["--block-size-reads", f"{1000}"]
        run_cmd += ["--out", f"{output_path}"]
        run_cmd += ["--threads", f"{self.threads}"]
        if filter_reads:
            run_cmd += ["--pos-filter"]

        return [run_cmd]

