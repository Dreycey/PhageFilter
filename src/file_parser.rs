use bio::alphabets::dna;
use bio::io::{fasta, fastq};
use flate2::read::GzDecoder;
use rayon::prelude::*;
use std::fs;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// User-facing format override for the `--format` CLI flag.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FormatOverride {
    Auto,
    Fasta,
    Fastq,
}

/// Detected sequence file format.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SequenceFormat {
    Fasta,
    Fastq,
}

/// Detect the sequence format of a file by inspecting its contents.
///
/// Strategy:
/// 1. If a format override is specified, use it directly.
/// 2. Peek at the first bytes of the (possibly gzip-decompressed) file.
///    - `>` indicates FASTA
///    - `@` indicates FASTQ
/// 3. Fall back to extension-based detection if content sniffing is inconclusive.
fn detect_format(filepath: &str, format_override: FormatOverride) -> SequenceFormat {
    match format_override {
        FormatOverride::Fasta => return SequenceFormat::Fasta,
        FormatOverride::Fastq => return SequenceFormat::Fastq,
        FormatOverride::Auto => {}
    }

    if let Ok(mut file) = File::open(filepath) {
        let mut header = [0u8; 2];
        if file.read_exact(&mut header).is_ok() {
            if header[0] == 0x1f && header[1] == 0x8b {
                // Gzip-compressed — decompress and peek at first byte
                file.seek(SeekFrom::Start(0)).unwrap();
                let mut decoder = GzDecoder::new(file);
                let mut first_byte = [0u8; 1];
                if decoder.read_exact(&mut first_byte).is_ok() {
                    return match first_byte[0] {
                        b'>' => SequenceFormat::Fasta,
                        b'@' => SequenceFormat::Fastq,
                        _ => format_from_extension(filepath),
                    };
                }
            } else {
                return match header[0] {
                    b'>' => SequenceFormat::Fasta,
                    b'@' => SequenceFormat::Fastq,
                    _ => format_from_extension(filepath),
                };
            }
        }
    }

    format_from_extension(filepath)
}

/// Infer format from the file extension, handling compound extensions like `.fq.gz`.
fn format_from_extension(filepath: &str) -> SequenceFormat {
    let path = Path::new(filepath);
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");

    let effective_ext = if ext.eq_ignore_ascii_case("gz") || ext.eq_ignore_ascii_case("gzip") {
        Path::new(path.file_stem().unwrap_or_default())
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
    } else {
        ext
    };

    match effective_ext {
        "fq" | "fastq" => SequenceFormat::Fastq,
        _ => SequenceFormat::Fasta,
    }
}

/// Open a file, transparently decompressing gzip if detected.
fn open_reader(filepath: &str) -> Box<dyn Read> {
    let mut file = File::open(filepath)
        .unwrap_or_else(|e| panic!("Failed to open '{}': {}", filepath, e));
    let mut magic = [0u8; 2];
    let is_gz = file.read_exact(&mut magic).is_ok() && magic[0] == 0x1f && magic[1] == 0x8b;
    file.seek(SeekFrom::Start(0)).unwrap();

    if is_gz {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    }
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
pub fn get_kmers(sequence: &[u8], &kmer_size: &usize) -> Vec<Vec<u8>> {
    if kmer_size > sequence.len() || kmer_size == 0 {
        return vec![];
    }
    let kmer_count = sequence.len() - kmer_size + 1;
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
    pub sequence: Option<Vec<u8>>,
    pub quality: Option<Vec<u8>>,
    pub id: String,
    pub kmers: Vec<Vec<u8>>,
}

impl DNASequence {
    pub fn new(
        sequence: Option<Vec<u8>>,
        quality: Option<Vec<u8>>,
        id: String,
        kmers: Vec<Vec<u8>>,
    ) -> DNASequence {
        DNASequence {
            id,
            kmers,
            sequence,
            quality,
        }
    }
}

pub enum RecordTypes {
    FastaRecords(fasta::Records<BufReader<Box<dyn Read>>>),
    FastqRecords(fastq::Records<BufReader<Box<dyn Read>>>),
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
    pub fn next_read(&mut self, kmer_size: usize, filtering: bool) -> Option<DNASequence> {
        match self {
            RecordTypes::FastaRecords(record) => {
                let record = record.next();
                if record.is_none() {
                    return None;
                }
                let seq_record = record.unwrap().unwrap();
                let id = seq_record.id().to_string();
                let mut sequence = Some(seq_record.seq().to_vec());
                let kmers = get_kmers(&sequence.as_ref().unwrap(), &kmer_size);
                if !filtering {
                    sequence = None;
                }
                Some(DNASequence::new(sequence, None, id, kmers))
            }
            RecordTypes::FastqRecords(record) => {
                let record = record.next();
                if record.is_none() {
                    return None;
                }
                let seq_record = record.unwrap().unwrap();
                let mut sequence = Some(seq_record.seq().to_vec());
                let id = seq_record.id().to_string();
                let kmers = get_kmers(&sequence.as_ref().unwrap(), &kmer_size);
                let mut quality = Some(seq_record.qual().to_vec());
                if !filtering {
                    sequence = None;
                    quality = None;
                }
                Some(DNASequence::new(sequence, quality, id, kmers))
            }
        }
    }
}

pub struct ReadQueue {
    filequeue: Vec<String>,
    block_size: usize,
    kmer_size: usize,
    curr_file: Option<Box<RecordTypes>>,
    filtering: bool,
    format_override: FormatOverride,
}

impl ReadQueue {
    pub fn get_next_file(&mut self) {
        self.curr_file = self.filequeue.pop().map(|filepath| {
            let reader = open_reader(&filepath);
            let format = detect_format(&filepath, self.format_override);
            match format {
                SequenceFormat::Fastq => Box::new(RecordTypes::FastqRecords(
                    fastq::Reader::new(reader).records(),
                )),
                SequenceFormat::Fasta => Box::new(RecordTypes::FastaRecords(
                    fasta::Reader::new(reader).records(),
                )),
            }
        });
    }

    pub fn next_block(&mut self) -> Vec<DNASequence> {
        let mut read_block: Vec<DNASequence> = vec![];
        if self.curr_file.is_none() {
            self.get_next_file();
        }
        while self.block_size > read_block.len() && self.curr_file.is_some() {
            if let Some(result) = self
                .curr_file
                .as_mut()
                .unwrap()
                .next_read(self.kmer_size, self.filtering)
            {
                read_block.push(result);
            } else {
                self.get_next_file();
            }
        }
        read_block
    }

    pub fn new(file_path: &str, block_size: usize, kmer_size: usize, filtering: bool) -> Self {
        Self::with_format(file_path, block_size, kmer_size, filtering, FormatOverride::Auto)
    }

    pub fn with_format(
        file_path: &str,
        block_size: usize,
        kmer_size: usize,
        filtering: bool,
        format_override: FormatOverride,
    ) -> Self {
        ReadQueue {
            filequeue: get_file_names(file_path),
            block_size,
            kmer_size,
            curr_file: None,
            filtering,
            format_override,
        }
    }

    /// Detect the sequence format of the first file in the queue.
    /// Useful for deciding the output file format before processing begins.
    pub fn peek_format(&self) -> SequenceFormat {
        self.filequeue
            .last()
            .map(|f| detect_format(f, self.format_override))
            .unwrap_or(SequenceFormat::Fasta)
    }
}

const SEQ_EXTENSIONS: &[&str] = &["fa", "fasta", "fna", "fsa", "fas", "fq", "fastq"];
const COMPRESSED_EXTENSIONS: &[&str] = &["gz", "gzip"];

fn get_file_names(file_path: &str) -> Vec<String> {
    let path = Path::new(file_path);

    if path.metadata().unwrap().is_file() {
        return vec![file_path.to_string()];
    }

    fs::read_dir(file_path)
        .unwrap()
        .filter_map(Result::ok)
        .map(|dir_entry| dir_entry.path())
        .filter(|path| has_supported_extension(path))
        .map(|path| path.to_string_lossy().into_owned())
        .collect()
}

/// Check if a path has a supported sequence file extension,
/// including compound extensions like `.fasta.gz`.
fn has_supported_extension(path: &Path) -> bool {
    let ext = match path.extension().and_then(|e| e.to_str()) {
        Some(e) => e,
        None => return false,
    };

    if SEQ_EXTENSIONS.contains(&ext) {
        return true;
    }

    // For .gz/.gzip files, check the inner extension
    if COMPRESSED_EXTENSIONS.contains(&ext) {
        if let Some(stem) = path.file_stem() {
            if let Some(inner_ext) = Path::new(stem).extension().and_then(|e| e.to_str()) {
                return SEQ_EXTENSIONS.contains(&inner_ext);
            }
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    fn create_temp_dir(test_name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("phagefilter_{}_{}", test_name, std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_fasta(path: &std::path::Path, id: &str, seq: &str) {
        let mut f = File::create(path).unwrap();
        writeln!(f, ">{}\n{}", id, seq).unwrap();
    }

    fn write_fastq(path: &std::path::Path, id: &str, seq: &str) {
        let mut f = File::create(path).unwrap();
        let qual = "+".repeat(seq.len()).chars().take(seq.len()).collect::<String>();
        // Use 'I' as quality score for every base
        writeln!(f, "@{}\n{}\n+\n{}", id, seq, "I".repeat(seq.len())).unwrap();
    }

    fn write_gzipped(path: &std::path::Path, content: &[u8]) {
        let f = File::create(path).unwrap();
        let mut encoder = GzEncoder::new(f, Compression::default());
        encoder.write_all(content).unwrap();
        encoder.finish().unwrap();
    }

    #[test]
    fn test_get_kmers() {
        // Cannot get a kmer sequence from an empty vec
        assert_eq!(get_kmers(&[], &1), Vec::<Vec<u8>>::new());

        // Cannot get kmers of length 0
        assert_eq!(get_kmers(&[1, 2, 3], &0), Vec::<Vec<u8>>::new());

        assert_eq!(
            get_kmers(&[1, 2, 3], &1),
            vec![vec![1], vec![2], vec![3]]
        );
        assert_eq!(get_kmers(&[1, 2, 3], &2), vec![vec![1, 2], vec![2, 3]]);
        assert_eq!(get_kmers(&[1, 2, 3], &3), vec![vec![1, 2, 3]]);
    }

    #[test]
    fn test_get_lex_less() {
        // Test case 1: kmer and its reverse complement are equal
        let kmer = &vec![b'A', b'C', b'G', b'T'][..];
        assert_eq!(get_lex_less(kmer), vec![b'A', b'C', b'G', b'T']);

        // Test case 2: kmer is lexically less than its reverse complement
        let kmer2 = &vec![b'A', b'A', b'T', b'G'][..];
        assert_eq!(get_lex_less(kmer2), vec![b'A', b'A', b'T', b'G']);

        let kmer3 = &vec![b'G', b'T', b'A', b'G'][..];
        assert_eq!(get_lex_less(kmer3), vec![b'C', b'T', b'A', b'C']);
    }

    #[test]
    fn test_detect_format_fasta_content() {
        let dir = create_temp_dir("test_detect_format_fasta_content");
        let path = dir.join("test.dat");
        write_fasta(&path, "seq1", "ACGTACGT");
        assert_eq!(
            detect_format(path.to_str().unwrap(), FormatOverride::Auto),
            SequenceFormat::Fasta
        );
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_detect_format_fastq_content() {
        let dir = create_temp_dir("test_detect_format_fastq_content");
        let path = dir.join("test.dat");
        write_fastq(&path, "read1", "ACGTACGT");
        assert_eq!(
            detect_format(path.to_str().unwrap(), FormatOverride::Auto),
            SequenceFormat::Fastq
        );
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_detect_format_override() {
        let dir = create_temp_dir("test_detect_format_override");
        // Write a FASTA file but override to FASTQ
        let path = dir.join("test.fasta");
        write_fasta(&path, "seq1", "ACGTACGT");
        assert_eq!(
            detect_format(path.to_str().unwrap(), FormatOverride::Fastq),
            SequenceFormat::Fastq
        );
        assert_eq!(
            detect_format(path.to_str().unwrap(), FormatOverride::Fasta),
            SequenceFormat::Fasta
        );
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_detect_format_gzipped_fasta() {
        let dir = create_temp_dir("test_detect_format_gzipped_fasta");
        let path = dir.join("test.fa.gz");
        let content = b">seq1\nACGTACGT\n";
        write_gzipped(&path, content);
        assert_eq!(
            detect_format(path.to_str().unwrap(), FormatOverride::Auto),
            SequenceFormat::Fasta
        );
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_detect_format_gzipped_fastq() {
        let dir = create_temp_dir("test_detect_format_gzipped_fastq");
        let path = dir.join("test.fq.gz");
        let content = b"@read1\nACGTACGT\n+\nIIIIIIII\n";
        write_gzipped(&path, content);
        assert_eq!(
            detect_format(path.to_str().unwrap(), FormatOverride::Auto),
            SequenceFormat::Fastq
        );
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_format_from_extension() {
        assert_eq!(format_from_extension("reads.fq"), SequenceFormat::Fastq);
        assert_eq!(format_from_extension("reads.fastq"), SequenceFormat::Fastq);
        assert_eq!(format_from_extension("reads.fq.gz"), SequenceFormat::Fastq);
        assert_eq!(format_from_extension("genome.fa"), SequenceFormat::Fasta);
        assert_eq!(format_from_extension("genome.fasta"), SequenceFormat::Fasta);
        assert_eq!(format_from_extension("genome.fna"), SequenceFormat::Fasta);
        assert_eq!(format_from_extension("genome.fa.gz"), SequenceFormat::Fasta);
        assert_eq!(format_from_extension("genome.fasta.gz"), SequenceFormat::Fasta);
    }

    #[test]
    fn test_has_supported_extension() {
        assert!(has_supported_extension(Path::new("test.fa")));
        assert!(has_supported_extension(Path::new("test.fasta")));
        assert!(has_supported_extension(Path::new("test.fna")));
        assert!(has_supported_extension(Path::new("test.fsa")));
        assert!(has_supported_extension(Path::new("test.fas")));
        assert!(has_supported_extension(Path::new("test.fq")));
        assert!(has_supported_extension(Path::new("test.fastq")));
        assert!(has_supported_extension(Path::new("test.fq.gz")));
        assert!(has_supported_extension(Path::new("test.fasta.gz")));
        assert!(has_supported_extension(Path::new("test.fna.gzip")));
        assert!(!has_supported_extension(Path::new("test.txt")));
        assert!(!has_supported_extension(Path::new("test.gz")));
        assert!(!has_supported_extension(Path::new("test")));
    }

    #[test]
    fn test_read_queue_gzipped_fasta() {
        let dir = create_temp_dir("test_read_queue_gzipped_fasta");
        let path = dir.join("genome.fa.gz");
        let content = b">seq1\nACGTACGTACGTACGTACGTACGTACGT\n";
        write_gzipped(&path, content);

        let mut queue = ReadQueue::new(path.to_str().unwrap(), 10, 5, false);
        let block = queue.next_block();
        assert_eq!(block.len(), 1);
        assert_eq!(block[0].id, "seq1");
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_read_queue_gzipped_fastq() {
        let dir = create_temp_dir("test_read_queue_gzipped_fastq");
        let path = dir.join("reads.fq.gz");
        let seq = "ACGTACGTACGTACGTACGTACGTACGT";
        let content = format!("@read1\n{}\n+\n{}\n", seq, "I".repeat(seq.len()));
        write_gzipped(&path, content.as_bytes());

        let mut queue = ReadQueue::new(path.to_str().unwrap(), 10, 5, false);
        let block = queue.next_block();
        assert_eq!(block.len(), 1);
        assert_eq!(block[0].id, "read1");
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_directory_scan_includes_gz_files() {
        let dir = create_temp_dir("test_directory_scan_includes_gz_files");

        // Create various files
        write_fasta(&dir.join("genome.fa"), "s1", "ACGT");
        write_fastq(&dir.join("reads.fq"), "r1", "ACGT");
        write_gzipped(&dir.join("compressed.fasta.gz"), b">s2\nACGT\n");
        File::create(dir.join("notes.txt")).unwrap();
        File::create(dir.join("random.gz")).unwrap();

        let files = get_file_names(dir.to_str().unwrap());
        assert!(files.iter().any(|f| f.ends_with("genome.fa")));
        assert!(files.iter().any(|f| f.ends_with("reads.fq")));
        assert!(files.iter().any(|f| f.ends_with("compressed.fasta.gz")));
        assert!(!files.iter().any(|f| f.ends_with("notes.txt")));
        // random.gz has no inner sequence extension, should be excluded
        assert!(!files.iter().any(|f| f.ends_with("random.gz")));

        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_fastq_quality_preserved() {
        let dir = create_temp_dir("test_fastq_quality_preserved");
        let path = dir.join("reads.fq");
        let seq = "ACGTACGTACGTACGTACGTACGTACGT";
        let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIII";
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, "@read1\n{}\n+\n{}", seq, qual).unwrap();
        }
        // filtering=true to retain sequence and quality
        let mut queue = ReadQueue::new(path.to_str().unwrap(), 10, 5, true);
        let block = queue.next_block();
        assert_eq!(block.len(), 1);
        assert_eq!(block[0].id, "read1");
        assert!(block[0].quality.is_some());
        assert_eq!(
            std::str::from_utf8(block[0].quality.as_ref().unwrap()).unwrap(),
            qual,
        );
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_fasta_has_no_quality() {
        let dir = create_temp_dir("test_fasta_has_no_quality");
        let path = dir.join("genome.fa");
        write_fasta(&path, "seq1", "ACGTACGTACGTACGTACGTACGTACGT");
        let mut queue = ReadQueue::new(path.to_str().unwrap(), 10, 5, true);
        let block = queue.next_block();
        assert_eq!(block.len(), 1);
        assert!(block[0].quality.is_none());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn test_peek_format() {
        let dir = create_temp_dir("test_peek_format");
        write_fasta(&dir.join("genome.fa"), "s1", "ACGT");
        write_fastq(&dir.join("reads.fq"), "r1", "ACGT");

        let queue_fa = ReadQueue::new(dir.join("genome.fa").to_str().unwrap(), 10, 5, false);
        assert_eq!(queue_fa.peek_format(), SequenceFormat::Fasta);

        let queue_fq = ReadQueue::new(dir.join("reads.fq").to_str().unwrap(), 10, 5, false);
        assert_eq!(queue_fq.peek_format(), SequenceFormat::Fastq);

        let _ = fs::remove_dir_all(&dir);
    }
}
