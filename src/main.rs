mod bloom_filter;
mod bloom_tree;
mod cache;
mod file_parser;
mod query;
mod result_map;
use clap::{arg, Parser, Subcommand, ValueEnum};
use clap_verbosity_flag::Verbosity;
use rayon::prelude::*;
use std::fs;
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// CLI-level format selection, mapped to `file_parser::FormatOverride`.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum FormatArg {
    /// Auto-detect format by inspecting file contents (default)
    Auto,
    /// Force FASTA parsing
    Fasta,
    /// Force FASTQ parsing
    Fastq,
}

impl From<FormatArg> for file_parser::FormatOverride {
    fn from(arg: FormatArg) -> Self {
        match arg {
            FormatArg::Auto => file_parser::FormatOverride::Auto,
            FormatArg::Fasta => file_parser::FormatOverride::Fasta,
            FormatArg::Fastq => file_parser::FormatOverride::Fastq,
        }
    }
}

#[derive(Parser)]
#[command(name = "PhageFilter")]
#[command(author = "Dreycey Albin & Kirby Linvill")]
#[command(version = "2.0")]
#[command(about = "A fast, simple and memory efficient metagenomic filtering tool.", long_about = None)]
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
        db_path: PathBuf,
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
        /// Input file format. Auto-detects by default via content sniffing.
        #[arg(required = false, default_value_t = FormatArg::Auto, short = 'F', long, value_enum)]
        format: FormatArg,
    },
    /// Adds genomes to an already built BloomFilter.
    Add {
        /// Path to genomes file or directory. (Fasta)
        #[arg(required = true, short, long)]
        genomes: String,
        /// Path to store the tree to disk.
        #[arg(required = true, short, long)]
        db_path: PathBuf,
        /// Number of threads to use to build the bloom tree.
        #[arg(required = false, default_value_t = 4, short, long)]
        threads: usize,
        /// Size of the LRU cache. (how many BFs in memory at once.)
        #[arg(required = false, default_value_t = 10, short, long)]
        cache_size: usize,
        /// Input file format. Auto-detects by default via content sniffing.
        #[arg(required = false, default_value_t = FormatArg::Auto, short = 'F', long, value_enum)]
        format: FormatArg,
    },
    /// Queries a set of reads. (ran after building the bloom tree)
    Query {
        /// Path to read file or directory of reads. (Fasta or Fastq, or dirs with both)
        #[arg(required = true, short, long)]
        reads: String,
        /// Path to output directory.
        #[arg(required = true, short, long)]
        out: PathBuf,
        /// Path of the tree stored to disk.
        #[arg(required = true, short, long)]
        db_path: PathBuf,
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
        /// Depth of search within gSBT (once set, the tree depth is limited to speed search)
        #[arg(required = false, long)]
        search_depth: Option<usize>,
        /// Filter reads matching genomes in the gSBT.
        #[arg(required = false, default_value_t = false, long)]
        pos_filter: bool,
        /// Filter reads NOT matching genomes in the gSBT.
        #[arg(required = false, default_value_t = false, long)]
        neg_filter: bool,
        /// Input file format. Auto-detects by default via content sniffing.
        #[arg(required = false, default_value_t = FormatArg::Auto, short = 'F', long, value_enum)]
        format: FormatArg,
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
            format,
        } => {
            // initial message to show used parameters.
            log::info!(
                "\n Build input-  \n\tdb path:{:?} \n\tthreads:{} \n\tkmersize:{} \n",
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
                file_parser::ReadQueue::with_format(genomes, 1, *kmer_size, false, (*format).into());
            let mut genome_block: Vec<file_parser::DNASequence> = genome_queue.next_block();

            // create a new cache
            let bloomfilter_cache = cache::BFLruCache::new(*cache_size, db_path.clone());

            // build: bloom trees
            println!("Building the SBT...");
            let mut bloom_tree = bloom_tree::BloomTree::new(
                *kmer_size,
                db_path,
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
            bloom_tree.save(db_path);
            println!("Finished.");
        }

        Commands::Add {
            genomes,
            db_path,
            threads,
            cache_size,
            format,
        } => {
            // initial message to show used parameters.
            log::info!(
                "\n Build input-  \n\tdb:{:?} \n\tthreads:{}\n",
                db_path,
                threads
            );

            // number of threads to run
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            // build: bloom trees
            println!("Adding new genomes to the SBT...");

            // create a new cache
            let bloomfilter_cache = cache::BFLruCache::new(*cache_size, db_path.clone());

            // load bloom tree from disk
            let mut bloom_tree: bloom_tree::BloomTree =
                bloom_tree::BloomTree::load(db_path, bloomfilter_cache);

            // obtain genomes from fasta/fastq files
            let mut genome_queue: file_parser::ReadQueue =
                file_parser::ReadQueue::with_format(genomes, 1, bloom_tree.kmer_size, false, (*format).into());
            let mut genome_block: Vec<file_parser::DNASequence> = genome_queue.next_block();

            while !genome_block.is_empty() {
                for genome in genome_block {
                    bloom_tree.insert(&genome);
                }
                genome_block = genome_queue.next_block();
            }

            // save tree to disk
            bloom_tree.save(db_path);
            println!("Finished.");
        }

        Commands::Query {
            reads,
            out,
            db_path,
            threads,
            block_size_reads,
            filter_threshold: cutoff_threshold,
            cache_size,
            search_depth,
            pos_filter,
            neg_filter,
            format,
        } => {
            // initial message to show used parameters.
            log::info!(
                "\n Query input- \n\treads:{} \n\tthreads:{}, \n\tout:{:?}, \n\tdb_path:{:?}, \n\tthreads:{} \n\tcutoff_threshold:{} \n",
                reads, threads, out, db_path, threads, cutoff_threshold
            );

            // number of threads to run
            rayon::ThreadPoolBuilder::new()
                .num_threads(*threads)
                .build_global()
                .unwrap();

            // create a new cache
            let bloomfilter_cache = cache::BFLruCache::new(*cache_size, db_path.clone());

            // load bloom tree from disk
            let mut bloom_tree: bloom_tree::BloomTree =
                bloom_tree::BloomTree::load(db_path, bloomfilter_cache);

            // create a result map
            let mut result_map = result_map::ResultMap::new();

            // parse reads
            println!("Querying reads...");
            println!(
                "Filtering settings: positive={}; negative={}",
                pos_filter, neg_filter
            );
            let filtering_option = *pos_filter || *neg_filter;

            // changing search depth (i.e. prunes the tree to a search depth)
            if let Some(depth) = search_depth {
                if !filtering_option {
                    println!("If using a search depth, use a filtering flag (--pos-filter or --neg-filter, or both!)");
                }
                println!("Search depth settings: {}", depth);
                bloom_tree.prune_tree(*depth);
            }

            // create a read buffer
            let mut readqueue = file_parser::ReadQueue::with_format(
                reads,
                *block_size_reads,
                bloom_tree.kmer_size,
                filtering_option,
                (*format).into(),
            );

            // create an output directory
            let _ = create_and_overwrite_directory(out);

            // detect input format to match output format
            let input_is_fastq = readqueue.peek_format() == file_parser::SequenceFormat::Fastq;
            let filter_ext = if input_is_fastq { "fq" } else { "fa" };

            // open output files
            let pos_filter_file = if *pos_filter {
                Some(Arc::new(Mutex::new(
                    File::create(out.join(format!("POS_FILTERING.{}", filter_ext))).unwrap(),
                )))
            } else {
                None
            };
            let neg_filter_file = if *neg_filter {
                Some(Arc::new(Mutex::new(
                    File::create(out.join(format!("NEG_FILTERING.{}", filter_ext))).unwrap(),
                )))
            } else {
                None
            };

            // Check for presence in the bloom tree; block-by-block
            let mut read_block: Vec<file_parser::DNASequence> = readqueue.next_block();
            while !read_block.is_empty() {
                // query the read batch.
                bloom_tree = query::query_batch(
                    bloom_tree,
                    &read_block,
                    *cutoff_threshold,
                    &mut result_map,
                );

                // add reads to outfile
                if filtering_option {
                    read_block.par_iter().for_each(|read| {
                        let seq =
                            String::from_utf8(read.sequence.as_ref().unwrap().to_ascii_uppercase())
                                .unwrap();
                        if result_map.read_mapped(&read.id) {
                            if let Some(pos_file) = pos_filter_file.as_ref() {
                                let mut pos_file = pos_file.lock().unwrap();
                                let id = result_map.get_ext_id(&read.id);
                                write_record(&mut *pos_file, &id, &seq, read.quality.as_deref());
                            }
                        } else if let Some(neg_file) = neg_filter_file.as_ref() {
                            let mut neg_file = neg_file.lock().unwrap();
                            write_record(&mut *neg_file, &read.id, &seq, read.quality.as_deref());
                        }
                    });
                }

                // empty result map 
                result_map.empty_read_map();

                // get next read block
                read_block = readqueue.next_block();
            }

            // open output file to write to
            let mut out_file = File::create(out.join("CLASSIFICATION.csv")).unwrap();

            // save the number of reads mapped to leaf nodes (i.e. genomes in the file)
            query::save_leaf_counts(&bloom_tree.root.unwrap(), &mut out_file);
            println!("Finished.");
        }
    }
}  

fn create_and_overwrite_directory(dir_path: &Path) -> io::Result<()> {
    // Check if the directory exists
    if let Ok(metadata) = fs::metadata(dir_path) {
        if metadata.is_dir() {
            // Remove the directory if it exists
            fs::remove_dir_all(dir_path)?;
        }
    }

    // Create the directory
    fs::create_dir(dir_path)
}

/// Write a sequence record in FASTQ format (if quality is present) or FASTA format.
fn write_record(writer: &mut impl Write, id: &str, seq: &str, quality: Option<&[u8]>) {
    match quality {
        Some(qual) => {
            let qual_str = String::from_utf8_lossy(qual);
            writeln!(writer, "@{}\n{}\n+\n{}", id, seq, qual_str).unwrap();
        }
        None => {
            writeln!(writer, ">{}\n{}", id, seq).unwrap();
        }
    }
}
