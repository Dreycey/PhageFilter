"""
Description:
    This script provides classes for working with various bioinformatics tools, 
    such as PhageFilter, Kraken2, and FastViromeExplorer. Each class inherits 
    from the abstract base class ToolOp, which defines a common interface for
    parsing output, building the tool's database, and running the tool. The use 
    of an abstract base class enables the implementation of the Template Pattern, 
    which promotes code reusability and modularity.

Classes:
    ToolOp: Abstract base class for bioinformatics tools.
        parse_output: Parse the output of PhageFilter and return a dictionary of NCBI ID to read count.
        build: Build the tool's database.
        run: Run the tool and output command-line arguments.
    PhageFilter: Class for working with the PhageFilter tool.
    Kraken2: Class for working with the Kraken2 tool.
    FastViromeExplorer: Class for working with the FastViromeExplorer tool.
"""

# standard libraries
from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Tuple, Dict
import os




class ToolOp(ABC):
    """
    Abstract base class that defines a tool and its sub-operations.
    """

    @abstractmethod
    def parse_output():
        """
        Abstract method to parse an output file/directory (depends on the tool).
        """
        pass

    @abstractmethod
    def build():
        """
        Abstract method to run the tool, based on input arguments, it outputs a command-line array.
        """
        pass

    @abstractmethod
    def run():
        """
        Abstract method to run the tool, based on input arguments, it outputs a command-line array.
        """
        pass


class PhageFilter(ToolOp):

    def __init__(self, kmer_size: int, filter_thresh: float, threads=4):
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

    def parse_output(self, output_path: Path, genomes_path: Path = None, cuttoff=0.01) -> Dict[str, int]:
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
        with open(output_path) as out_file:
            line = out_file.readline()
            count = 0
            while line:
                name, count = line.strip("\n").split(",")
                name2counts[name] = int(count)
                line = out_file.readline()

        # filter based on read count threshold
        total_reads_classified = sum(name2counts.values())
        # TODO: filtering based on different read count thresholds for testing.
        filtered_name2counts = {
            k: v for k, v in name2counts.items() if v > cuttoff*total_reads_classified}

        return filtered_name2counts

    def build(self, db_path: Path, genomes_path: Path) -> List[List[str]]:
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
        build_cmd += ["-t", f"{self.threads}"]
        self.db_path = db_path
        return [build_cmd]

    def run(self, fasta_file: Path, output_path: Path):
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
        run_cmd += ["--cuttoff-threshold", f"{self.theta}"]
        run_cmd += ["--out", f"{output_path}"]
        run_cmd += ["-t", f"{self.threads}"]
        return [run_cmd]


class Kraken2(ToolOp):

    def __init__(self, kmer_size: int = 35, threads=4):
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
        build_cmds.append(["unzip", "new_taxdump.zip",
                          "-d", f"{db_path}taxonomy/"])

        # download NCBI accession2taxid
        nucl_wgs_ftp = 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        if not os.path.exists("nucl_gb.accession2taxid.gz"):
            build_cmds.append(["wget", f"{nucl_wgs_ftp}"])
        build_cmds.append(["gzip", "-d", "-k", "nucl_gb.accession2taxid.gz"])
        build_cmds.append(["mv", "nucl_gb.accession2taxid",
                          f"{db_path}taxonomy/nucl_gb.accession2taxid"])

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
        build_cmd += ["--fast-build"]
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
        with open(os.path.join(output_path, "FastViromeExplorer-final-sorted-abundance.tsv")) as out_file:
            line = out_file.readline()
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
        self.create_ref_db(genomes_path, ref_db_path=ref_db_path,
                           list_file_path=self.list_file_path)
        build_cmd = ["kallisto", "index", "-k",
                     f"{self.k}", "-i", f"{db_path}", f"{ref_db_path}"]

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

    def run(self, fasta_file: Path, output_path: Path) -> List[List[str]]:
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
        run_cmd = ["java", "-cp",
                   f"{self.tool_path}bin", f"FastViromeExplorer"]
        run_cmd += ["-i", f"{self.db_path}"]
        run_cmd += ["-l", f"{self.list_file_path}"]
        run_cmd += ["-1", f"{fasta_file}"]
        run_cmd += ["-o", f"{output_path}"]
        run_cmd += ["-cn", "1"]
        run_cmd += ["-co", "0.00001"]
        run_cmd += ["-reportRatio", "true"]
        return [run_cmd] 