"""
This contains the python wrapper for kraken2
"""
from bench.tools.tool_template import ToolOp
from pathlib import Path
from typing import List, Tuple, Dict
import os



class Kraken2(ToolOp):

    def __init__(self, kmer_size: int = 35, threads=1):
        """_summary_

        Args:
            kmer_size (int): _description_
            filter_thresh (float): _description_
            threads (int, optional): _description_. Defaults to 4.
        """
        self.k = kmer_size
        self.threads = threads
        self.db_path = None

    def parse_output(self, output_path: Path, genomes_path: Path) -> Dict[str, int]:
        """_summary_
        parses an output file/directory (depends on tool)
        returns a dictionary of the output of PhageFilter.

        Args:
            output_path (Path): Path where the output of PhageFilter
                                will be stored.

        Returns:
            Dict[str, int]: A map from NCBI ID to read count
        """
        taxid2ncbi = self.get_taxid2ncbi(genomes_path)
        name2counts = {}
        with open(output_path) as out_file:
            line = out_file.readline()
            count = 0
            while line:
                count, tax_level, taxid = line.strip("\n").split("\t")[2:5]
                if tax_level == "S" or tax_level == "S1" and int(count) > 0:
                    if taxid in taxid2ncbi:
                        for ncbi_id in taxid2ncbi[taxid]:
                            name2counts[ncbi_id] = int(count)
                line = out_file.readline()
        return name2counts

    def build(self, db_path: Path, genomes_path: Path, minimizer_len: int = 31, minimizer_spacing: int = 7):
        """_summary_
        Build the tools database.

        Args:
            db_path (Path): path to the database the tool needs.
            genomes_path (Path): path to the genomes the db will use for building.
            minimizer_len (int, optional): length for the Kraken minimizer. Defaults to 31.
            minimizer_spacing (int, optional): allows for errors in the kmers being mapped. Defaults to 7.

        Returns:
            List[List[str]]: a nested list of command line arguments.
        """
        build_cmds = []
        # get map from NCBI ID to taxid
        self.taxid2ncbi = self.get_taxid2ncbi(genomes_path)

        # make DB
        build_cmds.append(["mkdir", "-p", f"{db_path}/taxonomy/"])

        # download NCBI taxonomy
        taxdump_ftp = 'https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip'
        if not os.path.exists("new_taxdump.zip"):
            build_cmds.append(["wget", f"{taxdump_ftp}"])
        build_cmds.append(["unzip", "new_taxdump.zip", "-d", f"{db_path}taxonomy/"])

        # download NCBI accession2taxid
        nucl_wgs_ftp = 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        if not os.path.exists("nucl_gb.accession2taxid.gz"):
            build_cmds.append(["wget", f"{nucl_wgs_ftp}"])
        build_cmds.append(["gzip", "-d", "nucl_gb.accession2taxid.gz"])
        build_cmds.append(["mv", "nucl_gb.accession2taxid", f"{db_path}taxonomy/nucl_gb.accession2taxid"])

        # add genomic files
        for genome in os.listdir(genomes_path):
            build_cmd = ["kraken2-build", "--add-to-library",
                         f"{genomes_path / Path(genome)}"]
            build_cmd += ["--db", f"{db_path}"]
            build_cmds.append(build_cmd)

        # build command
        build_cmd = ["kraken2-build", "--build"]
        build_cmd += ["--db", f"{db_path}"]
        build_cmd += ["--kmer-len", f"{self.k}"]
        build_cmd += ["--minimizer-len", f"{minimizer_len}"]
        build_cmd += ["--minimizer-spaces", f"{minimizer_spacing}"]
        build_cmds.append(build_cmd)

        # update DB path
        self.db_path = db_path

        return build_cmds

    def run(self, fasta_file: Path, output_path: Path) -> List[List[str]]:
        """_summary_
        run tool, based on input arguments, it outputs a CMD-line array.

        Args:
            fasta_file (Path): Path to simualted reads.
            output_path (Path): Desired path for the output file.

        Returns:
            List[List[str]]: a nested list of command line arguments.
        """
        if not self.db_path:
            print("Must first build (Kraken2)")
            exit()
        run_cmd = ["kraken2", "--db", f"{self.db_path}",
                   f"{fasta_file}", "--report", f"{output_path}"]
        return [run_cmd]

    @staticmethod
    def get_taxid2ncbi(genomes_path: Path) -> Dict[str, int]:
        """_summary_
        This function is used for mapping taxonomy IDs to NCBI accessions. This is useful when
        comparing the output of Kraken2 to PhageFilter.

        Args:
            genomes_path (Path): Path to the directory of genomes used to build the Kraken2 DB

        Returns:
            Dict[str, int]: An output dictionary mapping NCBI taxomony IDs to NCBI IDs
        """
        ncbi2tax = {}
        for genome in os.listdir(genomes_path):
            with open(os.path.join(genomes_path, genome), 'r') as f:
                line = f.readline()
                while line:
                    if line.startswith(">"):
                        ncbi = line.strip(">").strip("\n").split(
                            "|kraken:taxid|")[1].strip()
                        taxid = line.strip(">").strip(
                            "\n").split(" ")[0].strip()
                        if ncbi in ncbi2tax:
                            ncbi2tax[ncbi].append(taxid)
                        else:
                            ncbi2tax[ncbi] = [taxid]
                    line = f.readline()
        return ncbi2tax