mod bloom_filter;
use bloom_filter::ASMS;
mod file_parser;
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::io::Write;
use std::str;

fn main() {
    // parse the command line arguments
    let args: Vec<String> = env::args().collect();
    let seq_file_path = args[1].clone().parse::<String>().unwrap();
    let read_file_path: String = args[2].clone().parse::<String>().unwrap();
    let out_file_path: String = args[3].clone().parse::<String>().unwrap();

    // number of threads to run
    rayon::ThreadPoolBuilder::new()
        .num_threads(7)
        .build_global()
        .unwrap();
    // create a bloom filter
    let mut bloom_filter = bloom_filter::get_bloom_filter();

    // obtain genomes from fasta/fastq files
    let parsed_genomes: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&seq_file_path);
    add_to_bloom(parsed_genomes, 20, &mut bloom_filter);

    // open output file to write to
    let out_file = File::create(out_file_path).unwrap();

    // parse reads.
    let parsed_reads: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&read_file_path);
    check_if_in_bloom_filter(parsed_reads, 20, &mut bloom_filter, &out_file);
    //parsed_reads.par_iter_mut().for_each(|p| check_if_in_bloom_filter(parsed_reads, 20, &mut bloom_filter););
}

fn add_to_bloom(
    parsed_genomes: Vec<file_parser::RecordTypes>,
    kmer_size: usize,
    bloom_filter: &mut bloom_filter::BloomFilter,
) {
    //let genome: file_parser::Record;
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
    kmer_size: usize,
    bloom_filter: &bloom_filter::BloomFilter,
    out_file: &File,
) {
    // for read in parsed_reads {
    //     println!("READ CHECK");
    //     check(&read, bloom_filter, kmer_size);
    // }

    parsed_reads
        .par_iter_mut()
        .for_each(|read| check(&read, bloom_filter, kmer_size, out_file));
}

fn check(
    read: &file_parser::RecordTypes,
    bloom_filter: &bloom_filter::BloomFilter,
    kmer_size: usize,
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
