"""
This contains the python wrapper for kraken2
"""
from bench.tools.tool_template import ToolOp
from pathlib import Path
from typing import List, Tuple, Dict
import os


class FastViromeExplorer(ToolOp):

    def __init__(self, tool_path, list_file_path, kmer_size: int = 31):
        """_summary_

        Args:
            kmer_size (int): _description_
            filter_thresh (float): _description_
            threads (int, optional): _description_. Defaults to 4.
        """
        self.k = kmer_size
        self.tool_path = tool_path
        self.list_file_path = list_file_path
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
        name2counts = {}
        with open(os.path.join(output_path, "abundance.tsv")) as out_file:
            out_file.readline() # skip first line
            line = out_file.readline()
            count = 0
            while line:
                ncbi_id = line.strip("\n").split("\t")[0]
                count = float(line.strip("\n").split("\t")[3])
                if count > 0:
                    name2counts[ncbi_id] = count
                line = out_file.readline()
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
        ref_db_path = "tmp_delete.fa"

        # create fasta for making DB
        self.create_ref_db(genomes_path, ref_db_path=ref_db_path, list_file_path=self.list_file_path)
        build_cmd = ["kallisto", "index", "-k", f"{self.k}", "-i", f"{db_path}", f"{ref_db_path}"]

        # update attributes.
        self.db_path = db_path

        return [build_cmd]

    def create_ref_db(self, genomes_path: Path, ref_db_path: Path, list_file_path: Path, fasta_line_len=80):
        """_summary_

        Args:
            genomes_path (Path): _description_
        """
        ref_db = open(ref_db_path, "w+")
        list_file = open(list_file_path, "w+")
        for genome in os.listdir(genomes_path):
            genome_path = os.path.join(genomes_path, genome)
            # parse fasta.
            genome_seq, length, ncbi_id = self.parse_fasta(genome_path)
            # save information to list file.
            list_file.write(f"{ncbi_id}\tN/A\tN/A\t{length}\n")
            # save genome to file.
            ref_db.write(f">{ncbi_id} {length}\n")
            for genome_index in range(0, length, fasta_line_len):
                ref_db.write(
                    genome_seq[genome_index:genome_index+fasta_line_len] + "\n")

    def parse_fasta(self, fasta_file: Path):
        """_summary_

        Args:
            fasta_file (Path): _description_
        """
        genome_seq = ""
        with open(fasta_file, "r") as fasta_file:
            line = fasta_file.readline()
            while line:
                if line.startswith(">"):
                    ncbi_id = line.strip(">").strip("\n").split(" ")[0]
                else:
                    genome_seq += line.strip("\n")
                line = fasta_file.readline()
        return genome_seq, len(genome_seq), ncbi_id

    def run(self, fasta_file: Path, output_path: Path, filter_reads=False) -> List[List[str]]:
        """_summary_
        run tool, based on input arguments, it outputs a CMD-line array.

        Args:
            fasta_file (Path): Path to simualted reads.
            output_path (Path): Desired path for the output file.

        Returns:
            List[List[str]]: a nested list of command line arguments.

        Notes:
            Usage:
            java -cp bin FastViromeExplorer -1 $read1File -2 $read2File -i $indexFile -o $outputDirectory
            -1: input .fastq file for read sequences (paired-end 1), mandatory field.
            -2: input .fastq file for read sequences (paired-end 2).
            -i: kallisto index file, mandatory field.
            -db: reference database file in fasta/fa format.
            -o: output directory. Default option is the project directory.
            -l: virus list containing all viruses present in the reference database along with their length.
            -cr: the value of ratio criteria, default: 0.3.
            -co: the value of coverage criteria, default: 0.1.
            -cn: the value of number of reads criteria, default: 10.
            -salmon: use salmon instead of kallisto, default: false. To use salmon pass '-salmon true' as parameter.
            -reportRatio: default: false. To get ratio pass '-reportRatio true' as parameter.
        """
        if not self.db_path or not self.list_file_path:
            print("Must first build (FastViromeExplorer)")
            exit()
        run_cmd = ["java", "-cp", f"{self.tool_path}bin", f"FastViromeExplorer"]
        run_cmd += ["-i", f"{self.db_path}"]
        run_cmd += ["-l", f"{self.list_file_path}"]
        run_cmd += ["-1", f"{fasta_file}"]
        run_cmd += ["-o", f"{output_path}"]
        run_cmd += ["-cn", "1"]
        run_cmd += ["-co", "0.001"]
        run_cmd += ["-reportRatio", "true"]
        return [run_cmd] 