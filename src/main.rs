mod bloom_filter;
mod bloom_tree;
mod file_parser;
mod query;
use clap::{arg, ArgMatches, Command};
use std::fs::File;
use std::path::Path;

fn main() {
    // get values from command line.
    let matches: ArgMatches = parse_cmdline();

    let seq_file_path: &String = matches
        .get_one::<String>("genomes")
        .expect("Genomes required");
    let read_file_path: &String = matches.get_one::<String>("reads").expect("Reads required");
    let out_file_path: &String = matches.get_one::<String>("out").expect("Output required");
    let thread_count: usize = matches
        .get_one::<String>("thread_count")
        .expect("Threads required")
        .parse::<usize>()
        .unwrap();
    let kmer_size: usize = matches
        .get_one::<String>("kmer_size")
        .expect("kmer size required")
        .parse::<usize>()
        .unwrap();
    let threshold: f32 = matches
        .get_one::<String>("threshold")
        .expect("threshold required")
        .parse::<f32>()
        .unwrap();

    // number of threads to run
    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .unwrap();

    // obtain genomes from fasta/fastq files
    let parsed_genomes: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&seq_file_path);

    // build: bloom tree
    let mut bloom_node = bloom_tree::create_bloom_tree(parsed_genomes, &kmer_size);
    //let save_dir = Path::new("./tree/");
    //let mut bloom_node: bloom_tree::BloomTree = bloom_tree::BloomTree::load(save_dir);
    // open output file to write to
    let mut out_file = File::create(out_file_path).unwrap();

    // parse reads and check for presence in the bloom tree.
    print!("Querying reads...");
    let parsed_reads: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&read_file_path);
    bloom_node = query::query_batch(bloom_node, parsed_reads, threshold);

    // save the number of reads mapped to leaf nodes (i.e. genomes in the file)
    query::save_leaf_counts(&bloom_node.root.unwrap(), &mut out_file);
}

fn parse_cmdline() -> ArgMatches {
    // parse the command line arguments
    let matches = Command::new("PhageFilter")
        .version("2.0")
        .author("Dreycey Albin <albindreycey@gmail.com>")
        .about("A fast, simple, and efficient method for taxonomic classification.")
        .arg(
            arg!(--genomes <VALUE>)
                .required(true)
                .short('g')
                .long("genomes")
                .help("Path to genomes file or directory. (Fasta)"),
        )
        .arg(
            arg!(--reads <VALUE>)
                .required(true)
                .short('r')
                .long("reads")
                .help(
                    "Path to read file or directory of reads. (Fasta or Fastq, or dirs with both)",
                ),
        )
        .arg(
            arg!(--out <VALUE>)
                .required(true)
                .short('o')
                .long("out")
                .help("Path to output file. (Fasta)"),
        )
        .arg(
            arg!(--thread_count <VALUE>)
                .required(false)
                .short('t')
                .long("threads")
                .help("Number of threads to use for read matching")
                .default_value("4"),
        )
        .arg(
            arg!(--kmer_size <VALUE>)
                .required(false)
                .short('k')
                .long("kmer_size")
                .help("Size of the kmer to use; use with caution!")
                .default_value("20"),
        )
        .arg(
            arg!(--threshold <VALUE>)
                .required(false)
                .short('q')
                .long("threshold")
                .help("Filtering theshold (Number of kmers needed to pass)")
                .default_value("1.0"),
        )
        .get_matches();

    return matches;
}
