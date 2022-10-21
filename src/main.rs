mod bloom_filter;
use bloom_filter::ASMS;
mod file_parser;
use clap::{arg, ArgMatches, Command};
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::str;

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

    // number of threads to run
    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .unwrap();

    // obtain genomes from fasta/fastq files
    let parsed_genomes: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&seq_file_path);

    // create a bloom filter
    let mut bloom_filter = bloom_filter::get_bloom_filter(parsed_genomes.len());

    // add genomes to the bloom filter
    add_to_bloom(parsed_genomes, &kmer_size, &mut bloom_filter);

    // open output file to write to
    let out_file = File::create(out_file_path).unwrap();

    // parse reads and check for presence in the bloom filter.
    let parsed_reads: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&read_file_path);
    check_if_in_bloom_filter(parsed_reads, &kmer_size, &mut bloom_filter, &out_file);
}

fn parse_cmdline() -> ArgMatches {
    // parse the command line arguments
    let matches = Command::new("PhageFilter")
        .version("1.0")
        .author("Dreycey Albin <albindreycey@gmail.com>")
        .about("A fast, simple, and efficient way to filter reads from metagenomic data")
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
                .default_value("4"), //.about("Number of threads to process"),
        )
        .arg(
            arg!(--kmer_size <VALUE>)
                .required(false)
                .short('k')
                .long("kmer_size")
                .help("Size of the kmer to use; use with caution!")
                .default_value("20"), //.about("Number of threads to process"),
        )
        .get_matches();

    return matches;
}

fn add_to_bloom(
    parsed_genomes: Vec<file_parser::RecordTypes>,
    kmer_size: &usize,
    bloom_filter: &mut bloom_filter::BloomFilter,
) {
    for genome in parsed_genomes {
        let sequence: Vec<u8> = file_parser::get_sequence(&genome);
        let id: &str = file_parser::get_id(&genome);
        let kmers = file_parser::get_kmers(&sequence, &kmer_size); // ATGC -> AT, TG, GC
        for kmer in kmers {
            bloom_filter.insert(&kmer);
        }
        println!("NEW GENOME: {}; length: {:#?}", id, sequence.len());
    }
}

fn check_if_in_bloom_filter(
    mut parsed_reads: Vec<file_parser::RecordTypes>,
    kmer_size: &usize,
    bloom_filter: &bloom_filter::BloomFilter,
    out_file: &File,
) {
    // for read in parsed_reads {
    //     println!("READ CHECK");
    //     check(&read, bloom_filter, kmer_size);
    // }

    parsed_reads
        .par_iter_mut()
        .for_each(|read| check_sequence(&read, bloom_filter, kmer_size, out_file));
}

fn check_sequence(
    read: &file_parser::RecordTypes,
    bloom_filter: &bloom_filter::BloomFilter,
    kmer_size: &usize,
    mut out_file: &File,
) {
    // map kmers
    let sequence: Vec<u8> = file_parser::get_sequence(read);
    let id: &str = file_parser::get_id(&read);
    let kmers = file_parser::get_kmers(&sequence, &kmer_size);
    let mut bit_vec: Vec<bool> = vec![];
    for kmer in kmers {
        bit_vec.push(bloom_filter.contains(&kmer));
    }

    // calculate the number of mapped kmers
    let kmers_in_bloom_filter: usize = bit_vec.iter().filter(|&n| *n == true).count();
    let kmer_freq: f32 = (kmers_in_bloom_filter as f32) / (bit_vec.len() as f32);

    // kmer frequency
    if kmer_freq > 0.01 {
        let s = match std::string::String::from_utf8(sequence) {
            Ok(v) => v,
            Err(e) => panic!("Invalid UTF-8 sequence: {}", e),
        };
        let mut owned_string: String = ">".to_owned();
        owned_string.push_str(id);
        owned_string.push_str("\n");
        owned_string.push_str(&s as &str);
        owned_string.push_str("\n");
        out_file
            .write(owned_string.as_bytes())
            .expect("problem writing to file!");
    }
}
