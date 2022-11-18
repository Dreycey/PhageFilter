/// Methods for parsing individual or a directory
/// of Fasta or Fastq files.
///
/// # Example Usage
///
/// ```rust
/// let parsed_genomes: Vec<file_parser::RecordTypes> = file_parser::get_genomes(&seq_file_path);
/// ```
use bio::alphabets::dna;
use bio::io::fasta;
use bio::io::fastq;
use std::fs;
use std::fs::metadata;

#[derive(Clone, Debug)]
pub enum RecordTypes {
    FastaRecord(fasta::Record),
    FastqRecord(fastq::Record),
}

/// Returns the lexographically smallest kmer between a
/// given kmer and its reverse compliment.
///
/// # Parameters
/// - `kmer`: The given kmer in u8
///
/// # Returns
/// - lexographically smallest kmer
///
/// # Panics
/// - N/A
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

/// Returns the kmers for a given sequence.
///
/// # Parameters
/// - `kmer`: The given kmer in u8
///
/// # Returns
/// - A vector containing the kmers for the given sequence, in u8.
///   the output vector is of type Vec<Vec<u8>>.
///
/// # Panics
/// - if the given sequence is not as long as the given kmer size. or if less/equal to zero.
pub fn get_kmers(sequence: &Vec<u8>, &kmer_size: &usize) -> Vec<Vec<u8>> {
    if kmer_size > sequence.len() || kmer_size <= 0 {
        // Can't get kmers of a size longer than the sequence
        // Can't get kmers of size 0
        return vec![];
    }

    let mut max_kmers: Vec<Vec<u8>> = vec![];
    // Rust ranges are exclusive on the end index, so adding 1 ensures we get the last kmer
    for kmer_ind in 0..sequence.len() - kmer_size + 1 {
        let kmer: Vec<u8> = sequence[kmer_ind..kmer_ind + kmer_size].to_vec();
        max_kmers.push(get_lex_less(kmer));
    }
    max_kmers
}

/// Returns the sequence from a RecordTypes enum.
///
/// # Parameters
/// - `genome`: May be a read or a genome, of enum type RecordTypes
///
/// # Returns
/// - The sequence from a RecordTypes enum, of type Vec<u8>
///
/// # Panics
/// - N/A
pub fn get_sequence(genome: &RecordTypes) -> Vec<u8> {
    let sequence: Vec<u8> = match genome {
        RecordTypes::FastaRecord(record) => record.seq().to_vec(),
        RecordTypes::FastqRecord(record) => record.seq().to_vec(),
    };
    sequence
}

/// Returns the sequence ID from a RecordTypes enum.
///
/// # Parameters
/// - `genome`: May be a read or a genome, of enum type RecordTypes
///
/// # Returns
/// - The sequence ID from a RecordTypes enum, of type Vec<u8>
///
/// # Panics
/// - N/A
pub fn get_id(genome: &RecordTypes) -> &str {
    let sequence: &str = match genome {
        RecordTypes::FastaRecord(record) => record.id(),
        RecordTypes::FastqRecord(record) => record.id(),
    };
    sequence
}

/// Obtains a list of genome objects given a genome directory,
/// or a genome file (i.e. fasta or fastq)
///
/// # Notes
/// 1. TODO: This stores all of the reads into memory, using iterator
///          may be prefered
///
/// # Parameters
/// - `genomes_path`: global path to the genome/read files (Fasta or Fastq, or directory).
///
/// # Returns
/// - A vector of RecordTypes (Vec<RecordTypes>)
///
/// # Panics
/// - N/A
pub fn get_genomes(genomes_path: &String) -> Vec<RecordTypes> {
    let supported_extensions = ["fa", "fasta", "fna", "fq", "fastq"];
    let md = metadata(genomes_path).unwrap();
    let mut records_arr: Vec<RecordTypes> = vec![];
    if md.is_dir() {
        let paths = fs::read_dir(genomes_path).unwrap().filter(|p| match p {
            Ok(s) => match s.path().extension() {
                Some(ext) => {
                    // Only use paths that have a supported extension
                    supported_extensions.contains(&ext.to_str().unwrap())
                }
                None => false,
            },
            // Propagate errors
            _ => true,
        });
        for path in paths {
            let mut tmp_arr: Vec<RecordTypes> =
                parse_genome_file(&path.unwrap().path().to_str().unwrap().to_string());
            records_arr.append(&mut tmp_arr);
        }
        if records_arr.len() == 0 {
            panic!("Could not load any genomes. May not have any files with a supported extension in the provided directory.");
        }
        return records_arr;
    } else {
        return parse_genome_file(genomes_path);
    }
}

/// Obtains a list of genome objects given a genome directory,
/// or a genome file (i.e. fasta or fastq)
///
/// # Notes
/// 1. TODO: This stores all of the reads into memory, using iterator
///          may be prefered
///
/// # Parameters
/// - `genomes_path`: global path to the genome/read files (Fasta or Fastq, or directory).
///
/// # Returns
/// - A vector of RecordTypes (Vec<RecordTypes>)
///
/// # Panics
/// - N/A
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

/// This method parses a fasta file.
///
/// # Notes
/// 1. TODO: This stores all of the reads into memory, using iterator
///          may be prefered
///
/// # Parameters
/// - `file_path`: global path to the individual genome/read file (Fasta).
///
/// # Returns
/// - A vector of RecordTypes (Vec<RecordTypes>)
///
/// # Panics
/// - N/A
fn fasta_parser(file_path: &String) -> Vec<RecordTypes> {
    let reader = fasta::Reader::from_file(file_path).unwrap();
    let mut records_arr: Vec<RecordTypes> = vec![];
    for result in reader.records() {
        records_arr.push(RecordTypes::FastaRecord(result.unwrap()));
    }
    return records_arr;
}

/// This method parses a fastq file.
///
/// # Notes
/// 1. TODO: This stores all of the reads into memory, using iterator
///          may be prefered
///
/// # Parameters
/// - `file_path`: global path to the individual genome/read file (Fastq).
///
/// # Returns
/// - A vector of RecordTypes (Vec<RecordTypes>)
///
/// # Panics
/// - N/A
fn fastq_parser(file_path: &String) -> Vec<RecordTypes> {
    let reader = fastq::Reader::from_file(file_path).unwrap();
    let mut records_arr: Vec<RecordTypes> = vec![];
    for result in reader.records() {
        records_arr.push(RecordTypes::FastqRecord(result.unwrap()));
    }
    return records_arr;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_kmers() {
        // Cannot get a kmer sequence from an empty vec
        assert_eq!(get_kmers(&vec![], &1), Vec::<Vec<u8>>::new());

        // Cannot get kmers of length 0
        assert_eq!(get_kmers(&vec![1, 2, 3], &0), Vec::<Vec<u8>>::new());

        assert_eq!(
            get_kmers(&vec![1, 2, 3], &1),
            vec![vec![1], vec![2], vec![3]]
        );
        assert_eq!(get_kmers(&vec![1, 2, 3], &2), vec![vec![1, 2], vec![2, 3]]);
        assert_eq!(get_kmers(&vec![1, 2, 3], &3), vec![vec![1, 2, 3]]);
    }
}
