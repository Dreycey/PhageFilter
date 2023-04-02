use bio::alphabets::dna;
use bio::io::{fasta, fastq};
use rayon::prelude::*;
use std::fs;
use std::fs::{metadata, File};
use std::io::BufReader;

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
///   the output vector is of type Vec<Vec<u8>>. An empty vector is
///   returned if the kmer is of size 0 or longer than the sequence.
///
/// # Panics
/// - if the given sequence is not as long as the given kmer size. or if less/equal to zero.
pub fn get_kmers(sequence: &Vec<u8>, &kmer_size: &usize) -> Vec<Vec<u8>> {
    if kmer_size > sequence.len() || kmer_size <= 0 {
        return vec![];
    }
    let kmer_count = sequence.len() - kmer_size + 1; // zero indexed.
                                                     // let mut max_kmers: Vec<Vec<u8>> = Vec::with_capacity(kmer_count);
                                                     // for kmer_ind in 0..kmer_count {
                                                     //     let kmer_slice = &sequence[kmer_ind..kmer_ind + kmer_size];
                                                     //     max_kmers.push(get_lex_less(kmer_slice.to_vec()));
                                                     // }
    let max_kmers: Vec<Vec<u8>> = (0..kmer_count)
        .into_par_iter()
        .map(|kmer_ind| {
            let kmer_slice = &sequence[kmer_ind..kmer_ind + kmer_size];
            get_lex_less(kmer_slice.to_vec())
        })
        .collect();
    max_kmers
}

#[derive(Debug, Clone)]
pub struct ReadClass {
    // pub sequence: Vec<u8>,
    //  pub id: String,
    pub kmers: Vec<Vec<u8>>,
}

#[derive(Debug)]
pub enum RecordTypes {
    FastaRecords(fasta::Records<BufReader<BufReader<File>>>),
    FastqRecords(fastq::Records<BufReader<BufReader<File>>>),
}

impl RecordTypes {
    /// Returns the lexographically smallest kmer between a
    /// given kmer and its reverse compliment.
    ///
    /// # Parameters
    /// - `kmer_size`: The kmer size used in the tree.
    ///
    /// # Returns
    /// - returns the read/record as a 'ReadClass'
    ///
    /// # Panics
    /// - N/A
    pub fn next_read(&mut self, kmer_size: usize) -> Option<ReadClass> {
        match self {
            RecordTypes::FastaRecords(record) => {
                let record = record.next();
                if record.is_none() {
                    return None;
                }
                let seq_record = record.unwrap().unwrap();
                let sequence = seq_record.seq().to_vec();
                let id = seq_record.id().to_string();
                let kmers = get_kmers(&sequence, &kmer_size);
                Some(ReadClass {
                    //   sequence,
                    // id,
                    kmers,
                })
            }
            RecordTypes::FastqRecords(record) => {
                let record = record.next();
                if record.is_none() {
                    return None;
                }
                let seq_record = record.unwrap().unwrap();
                let sequence = seq_record.seq().to_vec();
                let id = seq_record.id().to_string();
                let kmers = get_kmers(&sequence, &kmer_size);
                Some(ReadClass {
                    // sequence,
                    // id,
                    kmers,
                })
            }
        }
    }
}

#[derive(Debug)]
pub struct ReadQueue {
    filequeue: Vec<Box<String>>,
    block_size: usize,
    kmer_size: usize,
    curr_file: Option<Box<RecordTypes>>,
}

impl ReadQueue {
    pub fn get_next_file(&mut self) {
        if let Some(filepath) = self.filequeue.pop().as_deref() {
            let reader = BufReader::with_capacity(100, File::open(filepath).unwrap());
            if filepath.ends_with(".fq") || filepath.ends_with(".fastq") {
                self.curr_file = Some(Box::new(RecordTypes::FastqRecords(
                    fastq::Reader::new(reader).records(),
                )));
            } else {
                self.curr_file = Some(Box::new(RecordTypes::FastaRecords(
                    fasta::Reader::new(reader).records(),
                )));
            }
        } else {
            self.curr_file = None;
        }
    }

    pub fn next_block(&mut self) -> Vec<ReadClass> {
        let mut read_block: Vec<ReadClass> = vec![];
        if self.curr_file.is_none() {
            self.get_next_file();
        }
        while self.block_size > read_block.len() && !self.curr_file.is_none() {
            let result = self.curr_file.as_mut().unwrap().next_read(self.kmer_size);
            if result.is_none() {
                self.curr_file = None;
            } else {
                read_block.insert(0, result.unwrap());
            }
        }
        read_block
    }

    pub fn new(file_path: &str, block_size: usize, kmer_size: usize) -> Self {
        ReadQueue {
            filequeue: get_file_names(&file_path),
            block_size,
            kmer_size,
            curr_file: None,
        }
    }
}

fn get_file_names(file_path: &str) -> Vec<Box<String>> {
    let supported_extensions = ["fa", "fasta", "fna", "fq", "fastq"];

    // Get the metadata of the file/directory
    let md = metadata(file_path).unwrap();

    // If the input path is a file, return its name in a vector
    if md.is_file() {
        return vec![Box::new(file_path.to_string())];
    }

    // If the input path is a directory, get the names of all files with supported extensions
    let file_paths = fs::read_dir(file_path)
        .unwrap()
        .filter_map(|path| {
            // Filter out any directories or files with unsupported extensions
            path.ok().and_then(|dir_entry| {
                let path = dir_entry.path();
                path.extension().and_then(|ext| {
                    if supported_extensions.contains(&ext.to_str().unwrap()) {
                        Some(Box::new(path.to_str().unwrap().to_string()))
                    } else {
                        None
                    }
                })
            })
        })
        .collect::<Vec<Box<String>>>();

    file_paths
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

    #[test]
    fn test_get_lex_less() {
        // Test case 1: kmer and its reverse complement are equal
        let kmer = vec![b'A', b'C', b'G', b'T']; // <-- Lexically lesser
        let revcomp_kmer = vec![b'A', b'C', b'G', b'T'];
        assert_eq!(get_lex_less(kmer), vec![b'A', b'C', b'G', b'T']);

        // Test case 2: kmer is lexically less than its reverse complement
        let kmer2 = vec![b'A', b'A', b'T', b'G']; // <-- Lexically lesser
        let revcomp_kmer2 = vec![b'C', b'A', b'T', b'T'];
        assert_eq!(get_lex_less(kmer2), vec![b'A', b'A', b'T', b'G']);

        // Test case 3: kmer is lexically greater than its reverse complement
        let kmer3 = vec![b'G', b'T', b'A', b'G'];
        let revcomp_kmer3 = vec![b'C', b'T', b'A', b'C']; // <-- Lexically lesser
        assert_eq!(get_lex_less(kmer3), vec![b'C', b'T', b'A', b'C']);
    }
}
