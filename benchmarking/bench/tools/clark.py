"""
Module containing a wrapper around CLARK.

ARE YOU  WATCHING ME?

time how long iit takes me to code this.
"""
from collections import Counter
from bench.tools.tool_template import ToolOp
from pathlib import Path
from typing import List, Tuple, Dict
import os



class Clark(ToolOp):

    def __init__(self, tool_path, kmer_size: int=31, threads=1):
        """
        Initalize the wrapper for Clark.
        """
        self.kmer_size = kmer_size
        self.threads = threads
        self.tool_path = tool_path
        self.db_path = None
        self.targets_file = None

    def parse_output(self, output_path: Path, genomes_path: Path, filter_reads=False) -> Dict[str, int]:
        """
        Abstract method to parse an output file/directory (depends on the tool).
        """
        taxid2ncbi = self.get_taxid2ncbi(genomes_path)
        output_path = output_path + ".csv"

        # get read counts
        read_counters = Counter()
        total_reads = 0
        with open(output_path, "r") as filter_file:
            filter_file.readline() # skip headers
            line = filter_file.readline()
            while line:
                _, _, taxid = line.split(",")
                taxid = taxid.strip("\n")
                if taxid in taxid2ncbi:
                    ncbi_id = taxid2ncbi[taxid][0]
                    read_counters[ncbi_id] += 1
                line = filter_file.readline()
                total_reads += 1

        # update to abundances
        if not filter_reads:
            read_counters = {name: count/total_reads for name, count in read_counters.items()}

        return read_counters
    
    def build(self, db_path: Path, genomes_path: Path) -> List[List[str]]:
        """
        method to run the tool, based on input arguments, it outputs a command-line array.

        Description:
            This works by running the tool once on a fake input file. If the database 
            doesn't already exist, it will create it. If it does, it will skip the
            creation process and run the tool.

        Args:
            db_path (Path): path to the database the tool needs.
            genomes_path (Path): path to the genomes the db will use for building.
        """
        # set targets file path && update DB path
        self.targets_file = os.path.join(db_path, "targets.txt")
        self.db_path = db_path

        # create targets file
        paths2taxids = self.get_paths2taxid(genomes_path)
        with open(self.targets_file, "w") as target_file:
            for path2taxid in paths2taxids:
                target_file.write(f"{path2taxid[0]}\t{path2taxid[1]}\n")

        # Create fake reads file
        with open("FAKE_READS.fa", "w") as f:
            f.write(">FAKE_READS\n" + "A"*100 + "\n")

        return self.run(fasta_file="FAKE_READS.fa", output_path="FAKE_OUT")

    def get_paths2taxid(self, genome_dir: Path) -> Tuple[str, str]:
        """ Obtain the paths for the given genome with taxids.

        Args:
            genome_dir (Path): Path to the directory containing the genomes.

        Returns:
            Tuple[str, str]: A tuple containing the path to the genome and the taxid.
        """
        genome_paths = []
        for genome in os.listdir(genome_dir):
            genome_path = os.path.join(genome_dir, genome)
            with open(genome_path, 'r') as genome_file:
                line = genome_file.readline()
                taxid = line.strip(">").strip("\n").split("|kraken:taxid|")[1].strip()
            genome_paths.append((genome_path, taxid))
        return genome_paths

    def run(self, fasta_file: Path, output_path: Path, filter_reads=False, light_mode=False):
        """
        Run Clark on an input set of genomes for classification.
        """
        if self.db_path is None:
            raise Exception(f"The db_path must exist")
        elif self.targets_file is None:
            raise Exception(f"The targets file must exist")

        # run classification & read filtering
        run_cmd = [self.tool_path]
        run_cmd += ["-T", self.targets_file]
        run_cmd += ["-D", f"{self.db_path}"]
        run_cmd += ["-k", f"{self.kmer_size}"]
        run_cmd += ["-n", f"{self.threads}"]
        run_cmd += ["-O", f"{fasta_file}"]
        run_cmd += ["-R", f"{output_path}"]

        return [run_cmd]

if __name__ == "__main__":
    # instantiate a weapper for clark.
    tool_path = "/Users/dreyceyalbin/Desktop/PHAGEFILT/raw_phagefilter/PhageFilter/examples/genomes/viral_genome_dir/"
    clark_wrapper = Clark(tool_path=tool_path, kmer_size=31, threads=1)

    genome_path = "/Users/dreyceyalbin/Desktop/PHAGEFILT/raw_phagefilter/PhageFilter/examples/genomes/viral_genome_dir/"
    genomes_path = clark_wrapper.build(db_path=Path("db_path"), genomes_path=genome_path)

    test_file_path = "/Users/dreyceyalbin/Desktop/PHAGEFILT/raw_phagefilter/PhageFilter/examples/test_reads/sim_reads_c10000_n20_e0.01.fq"
    genomes_run = clark_wrapper.run(fasta_file="fasta_file_path", output_path="out_path", filter_reads=True)
