/// Methods for querying a bloom tree given a
/// vector of reads. These methods allow for querying
/// the bloom tree as well as saving leaf nodes that
/// have reads mapped.
///
///
/// # Example Usage
///
/// ```rust
/// let mut bloom_node = bloom_tree::create_bloom_tree(parsed_genomes, &kmer_size);
/// let mut out_file = File::create(out_file_path).unwrap();
/// bloom_node = query::query_batch(bloom_node, parsed_reads, threshold);
/// query::get_leaf_counts(&bloom_node.root.unwrap(), &mut out_file);
/// ```
use crate::bloom_filter::{BloomFilter, ASMS};
use crate::bloom_tree::{BloomNode, BloomTree};
use crate::file_parser;
use crate::result_map;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;

/// Used to evaluate if a given read contains the same kmers as a given
/// node. This is performed using the bloom filters. A given threshold
/// determines the stringency of the filter.
///
/// # Parameters
/// - `bloom_node`: A BloomNode node allocated on the heap (Box<BloomNode>)
/// - `read`: A reference to a single instance of file_parser::RecordTypes.
/// - `threshold`: The threshold for counting a read as a hit. (a float, 0-1)
/// - `kmer_size`: The kmer size used in the bloom tree.
///
/// # Returns
/// - True if the read contains the same kmers, false otherwise.
///
/// # Panics
/// - N/A
fn query_passes(
    bloom_filter: &BloomFilter,
    read: &file_parser::DNASequence,
    threshold: f32,
) -> bool {
    let num_matches = read
        .kmers
        .par_iter()
        .filter(|kmer| bloom_filter.contains(kmer))
        .count();
    num_matches >= (threshold * read.kmers.len() as f32).ceil() as usize
}

/// Given a BloomTree, this method uses `_query_batch` by passing
/// the root node of the bloom tree, given it exists. Likewise, it
/// passes the kmersize of the bloom tree to the `_query_batch`
/// method.
///
/// # Parameters
/// - `bloom_node`: A BloomNode node allocated on the heap (Box<BloomNode>)
/// - `read_set`: A vector of file_parser::RecordTypes, representing all reads.
/// - `threshold`: The threshol for counting a read as a hit. (a float, 0-1)
///
/// # Returns
/// - The BloomTree after mapping all given reads to leaf nodes.
///
/// # Panics
/// - N/A
pub(crate) fn query_batch(
    mut bloom_tree: BloomTree,
    read_set: &[file_parser::DNASequence],
    threshold: f32,
    result_map: &mut result_map::ResultMap,
) -> BloomTree {
    let root_match = bloom_tree.root.take();
    let read_refs: Vec<&file_parser::DNASequence> = read_set.iter().collect();
    bloom_tree.root = root_match.map(|root| _query_batch(
            &bloom_tree,
            root,
            &read_refs,
            threshold,
            result_map,
        ));
    bloom_tree
}

/// Given a root node, this method queries all children nodes
/// recursively, until a leaf node is found. This method performs
/// this operation in batch, allowing for a parrallelized query.
///
/// # Parameters
/// - `bloom_node`: A BloomNode node allocated on the heap (Box<BloomNode>)
/// - `read_set`: A vector of file_parser::RecordTypes, representing all reads.
/// - `threshold`: The threshol for counting a read as a hit. (a float, 0-1)
/// - `kmer_size`: The kmer size used in the bloom tree.
///
/// # Returns
/// - The BloomNode after evaluating itself and all its children (used recursively).
///
/// # Panics
/// - N/A
fn _query_batch(
    bloom_tree: &BloomTree,
    mut bloom_node: Box<BloomNode>,
    read_set: &[&file_parser::DNASequence],
    threshold: f32,
    result_map: &mut result_map::ResultMap,
) -> Box<BloomNode> {
    // get bloom filter
    let bf_bloom_node_b4 = bloom_tree
        .bf_cache
        .get_filter(&bloom_node.bloom_filter_path)
        .unwrap();
    let bf_bloom_node = bf_bloom_node_b4.read().unwrap();

    let pass: Vec<&file_parser::DNASequence> = read_set
        .par_iter()
        .filter(|&&read| query_passes(&bf_bloom_node, read, threshold))
        .copied()
        .collect();

    if !bloom_node.is_leafnode() {
        let queue = &pass;

        if !queue.is_empty() {
            if let Some(left_child) = bloom_node.left_child.take() {
                bloom_node.left_child = Some(_query_batch(
                    bloom_tree,
                    left_child,
                    &queue[..],
                    threshold,
                    result_map,
                ));
            }
            if let Some(right_child) = bloom_node.right_child.take() {
                bloom_node.right_child = Some(_query_batch(
                    bloom_tree,
                    right_child,
                    &queue[..],
                    threshold,
                    result_map,
                ));
            }
        }
    } else {
        bloom_node.mapped_reads += pass.len();

        // Create maping between reads and genomes.
        let genome_id = bloom_node.tax_id.clone().unwrap();
        if !pass.is_empty() {
            let first_read: &&file_parser::DNASequence = pass.first().unwrap();
            if first_read.sequence.is_some() {
                for read in pass.iter() {
                    result_map.add_read_map(read.id.clone(), genome_id.clone());
                }
            }
        }
    }

    bloom_node
}

/// Recursively traverses the bloom tree, saving lead node mappings to a
/// specified output file. Only saves the mapping information for leaves that
/// were mapped to at least once.
///
/// # Parameters
/// - `bloom_node`: A leaf node of the type `BloomNode`
/// - `output_file`: A `File` object that may be written to.
///
/// # Returns
/// - ()
///
/// # Panics
/// - If `get_leaf_counts` panics (if a leaf node does not have a taxonomic ID)
pub(crate) fn save_leaf_counts(bloom_node: &BloomNode, output_file: &mut File) {
    let nonzero_counts = get_leaf_counts(bloom_node)
        .into_iter()
        .filter(|(_, count)| *count > 0);
    for (id, count) in nonzero_counts {
        let tax_maps: String = format!("{},{}\n", id, count);
        output_file
            .write_all(tax_maps.as_bytes())
            .expect("problem writing to output file!");
    }
}

/// Recursively traverses the bloom tree, outputting the mapped reads for each
/// genome (identified by taxonomic ID). The genome-specific data is stored at
/// the leaves of the tree.
///
/// # Parameters
/// - `bloom_node`: A leaf node of the type `BloomNode`
///
/// # Returns
/// - Vector of (taxonomic ID, mapped read count) tuples
///
/// # Panics
/// - If a leaf node does not have a taxonomic ID
pub(crate) fn get_leaf_counts(bloom_node: &BloomNode) -> Vec<(&str, usize)> {
    if bloom_node.is_leafnode() {
        // TODO: we should use iterators instead of vectors since it will be
        //  much more efficient. May require building a new iterator type for
        //  our tree.
        vec![(
            bloom_node.tax_id.as_ref().unwrap().as_str(),
            bloom_node.mapped_reads,
        )]
    } else {
        let mut left_results = match &bloom_node.left_child {
            Some(child) => get_leaf_counts(child),
            None => vec![],
        };
        let mut right_results = match &bloom_node.right_child {
            Some(child) => get_leaf_counts(child),
            None => vec![],
        };
        left_results.append(&mut right_results);
        left_results
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bloom_tree::BloomTree;
    use crate::cache;
    use crate::file_parser;
    use crate::result_map;
    use std::fs;
    use std::path::PathBuf;

    fn setup_test_dir(name: &str) -> PathBuf {
        let dir = PathBuf::from(format!("tmp_testing_query/{}/", name));
        if dir.is_dir() {
            let _ = fs::remove_dir_all(&dir);
        }
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn sorted_counts(counts: &mut Vec<(&str, usize)>) {
        counts.sort_by_key(|(id, _)| id.to_string());
    }

    fn make_seq(seq: &str, id: &str, kmer_size: usize) -> file_parser::DNASequence {
        let bytes = seq.as_bytes().to_vec();
        let kmers = file_parser::get_kmers(&bytes, &kmer_size);
        file_parser::DNASequence::new(Some(bytes), None, id.to_string(), kmers)
    }

    fn build_four_genome_tree(kmer_size: usize, dir: &PathBuf) -> BloomTree {
        let bf_cache = cache::BFLruCache::new(10, dir.clone());
        let mut tree = BloomTree::new(kmer_size, dir, bf_cache, 0.001, 1000);
        for (id, seq) in [
            ("baseline", "ATCAG"),
            ("diff", "TTTAG"),
            ("onediff_first", "CTCAG"),
            ("onediff_mid", "ATTAG"),
        ] {
            tree.insert(&make_seq(seq, id, kmer_size));
        }
        tree
    }

    fn count_for(counts: &[(&str, usize)], name: &str) -> usize {
        counts.iter().find(|(id, _)| *id == name).map(|(_, c)| *c).unwrap_or(0)
    }

    #[test]
    fn test_query_passes() {
        let dir = setup_test_dir("passes");
        let kmer_size = 3;

        let bf_cache = cache::BFLruCache::new(10, dir.clone());
        let mut tree = BloomTree::new(kmer_size, &dir, bf_cache, 0.001, 1000);
        tree.insert(&make_seq("ATCGCA", "genome", kmer_size));

        let read_same = make_seq("ATCG", "read_same", kmer_size);
        let read_different = make_seq("AAAA", "read_different", kmer_size);

        let bf_path = tree.root.as_ref().unwrap().bloom_filter_path.clone();
        let bf_arc = tree.bf_cache.get_filter(&bf_path).unwrap();
        let bf = bf_arc.read().unwrap();

        assert!(query_passes(&bf, &read_same, 1.0));
        assert!(!query_passes(&bf, &read_different, 1.0));
        assert!(query_passes(&bf, &read_same, 0.0));
        assert!(query_passes(&bf, &read_different, 0.0));

        drop(bf);
        drop(bf_arc);
    }

    #[test]
    fn test_query_and_leaf_counts() {
        let dir = setup_test_dir("leaf_counts");
        let kmer_size = 5;
        let threshold = 0.1;

        let mut tree = build_four_genome_tree(kmer_size, &dir);

        let read = make_seq("ATCAG", "read_baseline", kmer_size);
        let mut result_map = result_map::ResultMap::new();
        tree = query_batch(tree, &[read], threshold, &mut result_map);

        let mut counts = get_leaf_counts(tree.root.as_ref().unwrap());
        sorted_counts(&mut counts);

        assert!(count_for(&counts, "baseline") >= 1, "baseline should match");
        assert_eq!(count_for(&counts, "diff"), 0, "diff should not match");
    }

    #[test]
    fn test_query_and_leaf_counts_smaller_kmer_size() {
        let dir = setup_test_dir("smaller_kmer");
        let kmer_size = 4;
        let threshold = 0.1;

        let mut tree = build_four_genome_tree(kmer_size, &dir);

        let read = make_seq("TCAG", "read_tcag", kmer_size);
        let mut result_map = result_map::ResultMap::new();
        tree = query_batch(tree, &[read], threshold, &mut result_map);

        let mut counts = get_leaf_counts(tree.root.as_ref().unwrap());
        sorted_counts(&mut counts);

        // "TCAG" kmer appears in baseline ("ATCAG") and onediff_first ("CTCAG")
        assert!(count_for(&counts, "baseline") >= 1, "baseline should match TCAG");
        assert_eq!(count_for(&counts, "diff"), 0, "diff should not match TCAG");
    }

    #[test]
    fn test_query_and_leaf_counts_multiple_reads() {
        let dir = setup_test_dir("multi_reads");
        let kmer_size = 4;
        let threshold = 0.51;

        let mut tree = build_four_genome_tree(kmer_size, &dir);

        let reads = vec![
            make_seq("TCAG", "read1", kmer_size),
            make_seq("ATCA", "read2", kmer_size),
        ];
        let mut result_map = result_map::ResultMap::new();
        tree = query_batch(tree, &reads, threshold, &mut result_map);

        let mut counts = get_leaf_counts(tree.root.as_ref().unwrap());
        sorted_counts(&mut counts);

        // With higher threshold, baseline should have the most hits
        let baseline_count = count_for(&counts, "baseline");
        let diff_count = count_for(&counts, "diff");
        assert!(baseline_count >= 1, "baseline should match at higher threshold");
        assert_eq!(diff_count, 0, "diff should not match");
    }

    #[test]
    fn test_query_and_leaf_counts_multiple_queries() {
        let dir = setup_test_dir("multi_queries");
        let kmer_size = 4;
        let threshold = 0.1;

        let mut tree = build_four_genome_tree(kmer_size, &dir);
        let mut result_map = result_map::ResultMap::new();

        // First query: read "TCAG"
        let read1 = make_seq("TCAG", "read_tcag", kmer_size);
        tree = query_batch(tree, &[read1], threshold, &mut result_map);

        // Second query: read "ATCA"
        let read2 = make_seq("ATCA", "read_atca", kmer_size);
        tree = query_batch(tree, &[read2], threshold, &mut result_map);

        let mut counts = get_leaf_counts(tree.root.as_ref().unwrap());
        sorted_counts(&mut counts);

        // Counts accumulate across queries; baseline should match both reads
        let baseline_count = count_for(&counts, "baseline");
        assert!(baseline_count >= 2, "baseline should accumulate counts from both queries, got {}", baseline_count);
        assert_eq!(count_for(&counts, "diff"), 0, "diff should not match either query");
    }
}
