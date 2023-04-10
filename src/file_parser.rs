use bio::alphabets::dna;
use bio::io::{fasta, fastq};
use rayon::prelude::*;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

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
fn get_lex_less(kmer: &[u8]) -> Vec<u8> {
    let kmer_revc: Vec<u8> = dna::revcomp(kmer);
    match kmer.cmp(&kmer_revc) {
        std::cmp::Ordering::Less => kmer.to_vec(),
        std::cmp::Ordering::Greater => kmer_revc,
        std::cmp::Ordering::Equal => kmer.to_vec(),
    }
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
            get_lex_less(kmer_slice)
        })
        .collect();
    max_kmers
}

#[derive(Debug, Clone)]
pub struct DNASequence {
    // pub sequence: Vec<u8>,
    pub id: String,
    pub kmers: Vec<Vec<u8>>,
}

impl DNASequence {
    pub fn new(sequence: Vec<u8>, id: String, kmer_size: usize) -> DNASequence {
        DNASequence {
            id,
            kmers: get_kmers(&sequence, &kmer_size),
        }
    }
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
    pub fn next_read(&mut self, kmer_size: usize) -> Option<DNASequence> {
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
                Some(DNASequence {
                    // sequence,
                    id,
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
                Some(DNASequence {
                    // sequence,
                    id,
                    kmers,
                })
            }
        }
    }
}

#[derive(Debug)]
pub struct ReadQueue {
    filequeue: Vec<String>,
    block_size: usize,
    kmer_size: usize,
    curr_file: Option<Box<RecordTypes>>,
}

impl ReadQueue {
    pub fn get_next_file(&mut self) {
        self.curr_file = self.filequeue.pop().map(|filepath| {
            let reader = BufReader::with_capacity(100, File::open(&filepath).unwrap());
            if filepath.ends_with(".fq") || filepath.ends_with(".fastq") {
                Box::new(RecordTypes::FastqRecords(
                    fastq::Reader::new(reader).records(),
                ))
            } else {
                Box::new(RecordTypes::FastaRecords(
                    fasta::Reader::new(reader).records(),
                ))
            }
        });
    }

    pub fn next_block(&mut self) -> Vec<DNASequence> {
        let mut read_block: Vec<DNASequence> = vec![];
        if self.curr_file.is_none() {
            self.get_next_file();
        }
        while self.block_size > read_block.len() && !self.curr_file.is_none() {
            if let Some(result) = self.curr_file.as_mut().unwrap().next_read(self.kmer_size) {
                read_block.push(result);
            } else {
                self.get_next_file();
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

fn get_file_names(file_path: &str) -> Vec<String> {
    let supported_extensions = ["fa", "fasta", "fna", "fq", "fastq"];
    let path = Path::new(file_path);

    if path.metadata().unwrap().is_file() {
        return vec![file_path.to_string()];
    }

    fs::read_dir(file_path)
        .unwrap()
        .filter_map(Result::ok)
        .map(|dir_entry| dir_entry.path())
        .filter(|path| {
            path.extension()
                .map(|ext| supported_extensions.contains(&ext.to_str().unwrap()))
                .unwrap_or(false)
        })
        .map(|path| path.to_string_lossy().into_owned())
        .collect()
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
        let kmer = &vec![b'A', b'C', b'G', b'T'][..]; // <-- Lexically lesser
        let revcomp_kmer = &vec![b'A', b'C', b'G', b'T'];
        assert_eq!(get_lex_less(&kmer), vec![b'A', b'C', b'G', b'T']);

        // Test case 2: kmer is lexically less than its reverse complement
        let kmer2 = &vec![b'A', b'A', b'T', b'G'][..]; // <-- Lexically lesser
        let revcomp_kmer2 = &vec![b'C', b'A', b'T', b'T'];
        assert_eq!(get_lex_less(&kmer2), vec![b'A', b'A', b'T', b'G']);

        // Test case 3: kmer is lexically greater than its reverse complement
        let kmer3 = &vec![b'G', b'T', b'A', b'G'][..];
        let revcomp_kmer3 = &vec![b'C', b'T', b'A', b'C']; // <-- Lexically lesser
        assert_eq!(get_lex_less(&kmer3), vec![b'C', b'T', b'A', b'C']);
    }
}
