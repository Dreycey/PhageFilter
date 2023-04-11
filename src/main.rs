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
        /// Size of the LRU cache. (how many BFs in memory at once.)
        #[arg(required = false, default_value_t = 10, short, long)]
        cache_size: usize,
        /// False positive rate per bloom filter. (impacts size of the bloom filter)
        #[arg(required = false, default_value_t = 0.001, short, long)]
        false_pos_rate: f32,
        /// Largest expected genome/chromosome size. (impacts size of the bloom filter)
        #[arg(required = false, default_value_t = 1000000, short, long)]
        largest_genome: u32,
    },
    /// Adds genomes to an already built BloomFilter.
    Add {
        /// Path to genomes file or directory. (Fasta)
        #[arg(required = true, short, long)]
        genomes: String,
        /// Path to store the tree to disk.
        #[arg(required = true, short, long)]
        db_path: String,
        /// Number of threads to use to build the bloom tree.
        #[arg(required = false, default_value_t = 4, short, long)]
        threads: usize,
        /// Size of the LRU cache. (how many BFs in memory at once.)
        #[arg(required = false, default_value_t = 10, short, long)]
        cache_size: usize,
    },
    /// Queries a set of reads. (ran after building the bloom tree)
    Query {
        /// Path to read file or directory of reads. (Fasta or Fastq, or dirs with both)
        #[arg(required = true, short, long)]
        reads: String,
        /// Path to output directory.
        #[arg(required = true, short, long)]
        out: String,
        /// Path of the tree stored to disk.
        #[arg(required = true, short, long)]
        db_path: String,
        /// Number of threads to use for read matching.
        #[arg(required = false, default_value_t = 4, short, long)]
        threads: usize,
        /// Size of the read blocks (how many reads in mem at a time).
        #[arg(required = false, default_value_t = 100, short, long)]
        block_size_reads: usize,
        /// Filtering theshold. (Fraction of kmers needed to pass)
        #[arg(required = false, default_value_t = 1.0, short, long)]
        filter_threshold: f32,
        /// Size of the LRU cache. (how many BFs in memory at once.)
        #[arg(required = false, default_value_t = 10, short, long)]
        cache_size: usize,
        /// Filter reads matching genomes in the gSBT.
        #[arg(required = false, default_value_t = false, long)]
        pos_filter: bool,
        /// Filter reads NOT matching genomes in the gSBT.
        #[arg(required = false, default_value_t = false, long)]
        neg_filter: bool,
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
            cache_size,
            false_pos_rate,
            largest_genome,
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
                file_parser::ReadQueue::new(&genomes, 1, *kmer_size, false);
            let mut genome_block: Vec<file_parser::DNASequence> = genome_queue.next_block();

            // create a new cache
            let bloomfilter_cache = cache::BFLruCache::new(*cache_size);

            // build: bloom trees
            print!("Building the SBT... \n");
            let mut bloom_tree = bloom_tree::BloomTree::new(
                *kmer_size,
                &PathBuf::from(db_path),
                bloomfilter_cache,
                *false_pos_rate,
                *largest_genome,
            );
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
        Commands::Add {
            genomes,
            db_path,
            threads,
            cache_size,
        } => {
            // initial message to show used parameters.
            log::info!(
                "\n Build input-  \n\tdb:{} \n\tthreads:{}\n",
                db_path,
                threads
            );
            // number of threads to run
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            // build: bloom trees
            print!("Adding new genomes to the SBT... \n");

            // create a new cache
            let bloomfilter_cache = cache::BFLruCache::new(*cache_size);

            // load bloom tree from disk
            let full_db_path: &Path = Path::new(db_path);
            let mut bloom_tree: bloom_tree::BloomTree =
                bloom_tree::BloomTree::load(full_db_path, bloomfilter_cache);

            // obtain genomes from fasta/fastq files
            let mut genome_queue: file_parser::ReadQueue =
                file_parser::ReadQueue::new(&genomes, 1, bloom_tree.kmer_size, false);
            let mut genome_block: Vec<file_parser::DNASequence> = genome_queue.next_block();

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
            block_size_reads,
            filter_threshold: cuttoff_threshold,
            cache_size,
            pos_filter,
            neg_filter,
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

            // create a new cache
            let bloomfilter_cache = cache::BFLruCache::new(*cache_size);

            // load bloom tree from disk
            let full_db_path: &Path = Path::new(db_path);
            let mut bloom_tree: bloom_tree::BloomTree =
                bloom_tree::BloomTree::load(full_db_path, bloomfilter_cache);

            // parse reads
            println!("filtering reads: {}", pos_filter);
            println!("negative filtering reads: {}", neg_filter);
            println!("Querying reads...");
            let filtering_option = (*pos_filter || *neg_filter);
            println!("filtering: {}", filtering_option);
            let mut readqueue = file_parser::ReadQueue::new(
                &reads,
                *block_size_reads,
                bloom_tree.kmer_size,
                filtering_option,
            );

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
            println!("Finished.");
        }
    }
}
