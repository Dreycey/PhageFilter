mod bloom_filter;
use bloom_filter::ASMS;
mod file_parser;
use std::env;

fn main() {
    println!("Hello, world!");
    // parse the command line arguments
    let args: Vec<String> = env::args().collect();
    let seq_file_path = args[1].clone().parse::<String>().unwrap();
    let read_file_path: String = args[2].clone().parse::<String>().unwrap();

    // create a bloom filter
    let mut bloom_filter = bloom_filter::get_bloom_filter();

    // obtain genomes from fasta/fastq files
    let parsed_genomes: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&seq_file_path);
    add_to_bloom(parsed_genomes, 20, &mut bloom_filter);

    // parse reads.
    let parsed_reads: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&read_file_path);
    check_if_in_bloom_filter(parsed_reads, 20, &mut bloom_filter);
}

fn add_to_bloom(
    parsed_genomes: Vec<file_parser::RecordTypes>,
    kmer_size: usize,
    bloom_filter: &mut bloom_filter::BloomFilter,
) {
    //let genome: file_parser::Record;
    for genome in parsed_genomes {
        println!("NEW GENOME");
        let sequence: Vec<u8> = match genome {
            file_parser::RecordTypes::FastaRecord(record) => record.seq().to_vec(),
            file_parser::RecordTypes::FastqRecord(record) => record.seq().to_vec(),
        };
        let kmers = sequence.windows(kmer_size); // ATGC -> AT, TG, GC
        for kmer in kmers {
            bloom_filter.insert(&kmer);
        }
        print!("{:#?}", bloom_filter.num_bits());
    }
}

fn check_if_in_bloom_filter(
    parsed_reads: Vec<file_parser::RecordTypes>,
    kmer_size: usize,
    bloom_filter: &mut bloom_filter::BloomFilter,
) {
    for read in parsed_reads {
        println!("READ CHECK");
        let sequence: Vec<u8> = match read {
            file_parser::RecordTypes::FastaRecord(record) => record.seq().to_vec(),
            file_parser::RecordTypes::FastqRecord(record) => record.seq().to_vec(),
        };
        let kmers = sequence.windows(kmer_size);
        let mut bit_vec: Vec<bool> = vec![];
        for kmer in kmers {
            bit_vec.push(bloom_filter.contains(&kmer));
        }
        let kmers_in_bloom_filter: usize = bit_vec.iter().filter(|&n| *n == true).count();
        let kmer_freq: f32 = (kmers_in_bloom_filter as f32) / (bit_vec.len() as f32);
        println!("{}", kmer_freq);
    }
}
