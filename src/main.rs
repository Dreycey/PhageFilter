mod bloom_filter;
mod bloom_tree;
mod cache;
mod file_parser;
mod query;
use clap::{arg, Parser, Subcommand};
use clap_verbosity_flag::Verbosity;
use log;
use std::fs::File;
use std::path::Path;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "MyApp")]
#[command(author = "Dreycey Albin & Kirby Linvill")]
#[command(version = "2.0")]
#[command(about = "A fast, simple and efficient metagenomic classification tool.", long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    /// Different command options. (Build or Query)
    #[command(subcommand)]
    command: Commands,
    /// Print log statements to stdout.
    #[command(flatten)]
    verbose: Verbosity,
}

#[derive(Subcommand)]
enum Commands {
    /// Builds the BloomTree.
    Build {
        /// Path to genomes file or directory. (Fasta)
        #[arg(required = true, short, long)]
        genomes: String,
        /// Path to store the tree to disk.
        #[arg(required = true, short, long)]
        db_path: String,
        /// Number of threads to use to build the bloom tree.
        #[arg(required = false, default_value_t = 4, short, long)]
        threads: usize,
        /// Size of the kmer to use; use with caution!
        #[arg(required = false, default_value_t = 20, short, long)]
        kmer_size: usize,
    },
    /// Queries a set of reads. (ran after building the bloom tree)
    Query {
        /// Path to read file or directory of reads. (Fasta or Fastq, or dirs with both)
        #[arg(required = true, short, long)]
        reads: String,
        /// Path to output file. (CSV)
        #[arg(required = true, short, long)]
        out: String,
        /// Path of the tree stored to disk.
        #[arg(required = true, short, long)]
        db_path: String,
        /// Number of threads to use for read matching.
        #[arg(required = false, default_value_t = 4, short, long)]
        threads: usize,
        /// Size of the read blocks (how many reads in mem at a time).
        #[arg(required = false, default_value_t = usize::MAX, short, long)]
        read_block_size: usize,
        /// Filtering theshold. (Fraction of kmers needed to pass)
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
            let mut genome_queue: file_parser::ReadQueue =
                file_parser::ReadQueue::new(&genomes, 1, *kmer_size);
            let mut genome_block: Vec<file_parser::DNASequence> = genome_queue.next_block();

            // build: bloom trees
            print!("Building the SBT... \n");
            let mut bloom_tree = bloom_tree::BloomTree::new(*kmer_size, &PathBuf::from(db_path));
            while !genome_block.is_empty() {
                for genome in genome_block {
                    bloom_tree.insert(&genome);
                }
                genome_block = genome_queue.next_block();
            }

            // save tree to disk
            let save_dir = Path::new(db_path);
            bloom_tree.save(save_dir);
            print!("Finished. \n");
        }
        Commands::Query {
            reads,
            out,
            db_path,
            threads,
            read_block_size,
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
            let mut bloom_tree: bloom_tree::BloomTree = bloom_tree::BloomTree::load(full_db_path);
            // parse reads
            print!("Querying reads... \n");
            let mut readqueue =
                file_parser::ReadQueue::new(&reads, *read_block_size, bloom_tree.kmer_size);
            // Check for presence in the bloom tree; block-by-block
            let mut read_block: Vec<file_parser::DNASequence> = readqueue.next_block();
            while !read_block.is_empty() {
                bloom_tree = query::query_batch(bloom_tree, &read_block, *cuttoff_threshold);
                read_block = readqueue.next_block();
            }
            // open output file to write to
            let mut out_file = File::create(out).unwrap();
            // save the number of reads mapped to leaf nodes (i.e. genomes in the file)
            query::save_leaf_counts(&bloom_tree.root.unwrap(), &mut out_file);
            print!("Finished. \n");
        }
    }
}
