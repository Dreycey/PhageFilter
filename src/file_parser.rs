use bio::alphabets::dna;
use bio::io::fasta;
use bio::io::fastq;
use std::fs;
use std::fs::metadata;

#[derive(Debug)]
pub enum RecordTypes {
    FastaRecord(fasta::Record),
    FastqRecord(fastq::Record),
}

fn get_lex_less(kmer: Vec<u8>) -> Vec<u8> {
    let kmer_revc: Vec<u8> = dna::revcomp(&kmer);
    for index in 0..kmer.len() {
        if kmer[index] < kmer_revc[index] {
            return kmer;
        } else if kmer[index] > kmer_revc[index] {
            return kmer_revc;
        }
    }
    return kmer; // if equal.
}

pub fn get_kmers(sequence: &Vec<u8>, &kmer_size: &usize) -> Vec<Vec<u8>> {
    let mut max_kmers: Vec<Vec<u8>> = vec![];
    for kmer_ind in 0..sequence.len() - kmer_size {
        let kmer: Vec<u8> = sequence[kmer_ind..=kmer_ind + kmer_size].to_vec();
        max_kmers.push(get_lex_less(kmer));
    }
    max_kmers
}

pub fn get_sequence(genome: &RecordTypes) -> Vec<u8> {
    let sequence: Vec<u8> = match genome {
        RecordTypes::FastaRecord(record) => record.seq().to_vec(),
        RecordTypes::FastqRecord(record) => record.seq().to_vec(),
    };
    sequence
}

pub fn get_id(genome: &RecordTypes) -> &str {
    let sequence: &str = match genome {
        RecordTypes::FastaRecord(record) => record.id(),
        RecordTypes::FastqRecord(record) => record.id(),
    };
    sequence
}

/// get_genomes(genomes_path: &String)
/// --
/// obtains a list of genome objects given a genome directory,
/// or a genome file (i.e. fasta or fastq)
///
/// Parameters
/// ----------
///
/// Returns
/// -------
///
/// Notes
/// -----
/// 1. TODO: This stores all of the reads into memory, using iterator
///          may be prefered
/// Examples
/// --------
///
pub fn get_genomes(genomes_path: &String) -> Vec<RecordTypes> {
    let md = metadata(genomes_path).unwrap();
    let mut records_arr: Vec<RecordTypes> = vec![];
    if md.is_dir() {
        let paths = fs::read_dir(genomes_path).unwrap();
        for path in paths {
            let mut tmp_arr: Vec<RecordTypes> =
                parse_genome_file(&path.unwrap().path().to_str().unwrap().to_string());
            records_arr.append(&mut tmp_arr);
        }
        return records_arr;
    } else {
        return parse_genome_file(genomes_path);
    }
}

/// parse_genome_file(genomes_file_path: &String)
/// --
/// obtains a list of genome objects given a genome file (i.e. fasta or fastq)
///
/// Parameters
/// ----------
///
/// Returns
/// -------
///
/// Notes
/// -----
/// 1. TODO: This stores all of the reads into memory, using iterator
///          may be prefered
/// Examples
/// --------
///
fn parse_genome_file(genomes_file_path: &String) -> Vec<RecordTypes> {
    //println!("\nParsing file: {}\n", genomes_file_path);
    if genomes_file_path.ends_with(".fa")
        || genomes_file_path.ends_with(".fasta")
        || genomes_file_path.ends_with(".fna")
    {
        return fasta_parser(&genomes_file_path);
    } else if genomes_file_path.ends_with(".fq") || genomes_file_path.ends_with(".fastq") {
        return fastq_parser(&genomes_file_path);
    } else {
        panic!("Incorrect file path: {}", genomes_file_path);
    }
}

/// fasta_parser(file_path: &String)
/// --
/// parses a fasta file.
///
/// Parameters
/// ----------
///
/// Returns
/// -------
///
/// Notes
/// -----
/// 1. TODO: This stores all of the reads into memory, using iterator
///          may be prefered
/// Examples
/// --------
///
fn fasta_parser(file_path: &String) -> Vec<RecordTypes> {
    let reader = fasta::Reader::from_file(file_path).unwrap();
    let mut records_arr: Vec<RecordTypes> = vec![];
    for result in reader.records() {
        records_arr.push(RecordTypes::FastaRecord(result.unwrap()));
    }
    return records_arr;
}

/// fastq_parser(file_path: &String)
/// --
/// parses a fastq file.
///
/// Parameters
/// ----------
///
/// Returns
/// -------
///
/// Notes
/// -----
/// 1. TODO: This stores all of the reads into memory, using iterator
///          may be prefered
/// Examples
/// --------
///
fn fastq_parser(file_path: &String) -> Vec<RecordTypes> {
    let reader = fastq::Reader::from_file(file_path).unwrap();
    let mut records_arr: Vec<RecordTypes> = vec![];
    for result in reader.records() {
        records_arr.push(RecordTypes::FastqRecord(result.unwrap()));
    }
    return records_arr;
}
