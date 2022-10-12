mod file_parser;
use std::env;

fn main() {
    println!("Hello, world!");
    // parse the command line arguments
    let args: Vec<String> = env::args().collect();
    let seq_file_path = args[1].clone().parse::<String>().unwrap();
    // obtain genomes from fasta/fastq files
    let parsed_genomes: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&seq_file_path);
}
