mod bloom_filter;
mod bloom_tree;
mod file_parser;
mod query;
use clap::{arg, Parser, Subcommand};
use clap_verbosity_flag::Verbosity;
use log;
use std::fs::File;
use std::path::Path;

#[derive(Parser)]
#[command(name = "MyApp")]
#[command(author = "Dreycey Albin & Kirby Linvill")]
#[command(version = "2.0")]
#[command(about = "A fast, simple and efficient metagenomic classification tool.", long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
    /// Print log statements to stdout.
    #[command(flatten)]
    verbose: Verbosity,
}

#[derive(Subcommand)]
enum Commands {
    /// Builds the BloomTree
    Build {
        /// Path to genomes file or directory. (Fasta)
        #[arg(required = true, short, long)]
        genomes: String,
        /// Path to store the tree to disk.
        #[arg(required = true, short, long)]
        db_path: String,
        /// Number of threads to use to build the bloom tree
        #[arg(required = false, default_value_t = 4, short, long)]
        threads: usize,
        /// Size of the kmer to use; use with caution!
        #[arg(required = false, default_value_t = 20, short, long)]
        kmer_size: usize,
    },
    /// Queries a set of reads (ran after building the bloom tree)
    Query {
        /// Path to read file or directory of reads. (Fasta or Fastq, or dirs with both)
        #[arg(required = true, short, long)]
        reads: String,
        /// Path to output file. (CSV)
        #[arg(required = true, short, long)]
        out: String,
        /// Path to store the tree to disk.
        #[arg(required = true, short, long)]
        db_path: String,
        /// Number of threads to use for read matching
        #[arg(required = false, default_value_t = 4, short, long)]
        threads: usize,
        /// Filtering theshold (Fraction of kmers needed to pass)
        #[arg(required = false, default_value_t = 1.0, short, long)]
        cuttoff_threshold: f32,
    },
}

fn main() {
    let cli = Cli::parse();
    // logging verbosity level
    env_logger::Builder::new()
        .filter_level(cli.verbose.log_level_filter())
        .init();
    log::info!("\n verbosity level: {} \n", cli.verbose.log_level_filter());

    // parse subcommands.
    match &cli.command {
        Commands::Build {
            genomes,
            db_path,
            threads,
            kmer_size,
        } => {
            // initial message to show used parameters.
            log::info!(
                "\n Build input-  \n\tdb:{} \n\tthreads:{} \n\tkmersize:{} \n",
                db_path,
                threads,
                kmer_size
            );
            // number of threads to run
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();
            // obtain genomes from fasta/fastq files
            let parsed_genomes: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&genomes);
            print!("Building the SBT... \n");
            // build: bloom tree
            let bloom_node = bloom_tree::create_bloom_tree(parsed_genomes, kmer_size);
            // save tree to disk
            let save_dir = Path::new(db_path);
            bloom_node.save(save_dir);
            print!("Finished. \n");
        }
        Commands::Query {
            reads,
            out,
            db_path,
            threads,
            cuttoff_threshold,
        } => {
            // initial message to show used parameters.
            log::info!(
                "\n Query input- \n\treads:{} \n\tthreads:{}, \n\tout:{}, \n\tdb_path:{}, \n\tthreads:{} \n\tcuttoff_threshold:{} \n",
                reads, threads, out, db_path, threads, cuttoff_threshold
            );
            // number of threads to run
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();
            // load bloom tree from disk
            let full_db_path: &Path = Path::new(db_path);
            let mut bloom_node: bloom_tree::BloomTree = bloom_tree::BloomTree::load(full_db_path);
            // parse reads and check for presence in the bloom tree.
            print!("Querying reads... \n");
            let parsed_reads: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&reads);
            bloom_node = query::query_batch(bloom_node, parsed_reads, *cuttoff_threshold);
            // open output file to write to
            let mut out_file = File::create(out).unwrap();
            // save the number of reads mapped to leaf nodes (i.e. genomes in the file)
            query::save_leaf_counts(&bloom_node.root.unwrap(), &mut out_file);
            print!("Finished. \n");
        }
    }
}
