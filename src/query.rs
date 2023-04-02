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
use crate::bloom_filter::ASMS;
use crate::bloom_tree::{BloomNode, BloomTree};
use crate::read_parser;
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
    bloom_node: &Box<BloomNode>,
    read: &read_parser::ReadClass,
    threshold: f32,
    kmer_size: usize,
) -> bool {
    let num_matches = read
        .kmers
        .par_iter()
        .filter(|kmer| bloom_node.bloom_filter.contains(kmer))
        .count();
    return num_matches >= (threshold * read.kmers.len() as f32).ceil() as usize;
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
    read_set: &[read_parser::ReadClass],
    threshold: f32,
) -> BloomTree {
    bloom_tree.root = match bloom_tree.root {
        Some(root) => Some(_query_batch(
            root,
            &read_set.iter().collect::<Vec<&read_parser::ReadClass>>()[..],
            threshold,
            bloom_tree.kmer_size,
        )),
        _ => None,
    };
    return bloom_tree;
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
    mut bloom_node: Box<BloomNode>,
    read_set: &[&read_parser::ReadClass],
    threshold: f32,
    kmer_size: usize,
) -> Box<BloomNode> {
    let pass: Vec<&read_parser::ReadClass> = read_set
        .par_iter()
        .filter(|&&read| query_passes(&bloom_node, read, threshold, kmer_size))
        .map(|&read| read) // Use map() to get the right item type
        .collect();

    if !bloom_node.is_leafnode() {
        let queue = &pass;

        if !queue.is_empty() {
            if let Some(left_child) = bloom_node.left_child.take() {
                bloom_node.left_child =
                    Some(_query_batch(left_child, &queue[..], threshold, kmer_size));
            }
            if let Some(right_child) = bloom_node.right_child.take() {
                bloom_node.right_child =
                    Some(_query_batch(right_child, &queue[..], threshold, kmer_size));
            }
        }
    } else {
        bloom_node.mapped_reads += pass.len();
    }

    return bloom_node;
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
            .write(tax_maps.as_bytes())
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
pub(crate) fn get_leaf_counts<'a>(bloom_node: &'a BloomNode) -> Vec<(&'a str, usize)> {
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

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::bloom_tree;
//     use crate::file_parser::RecordTypes;
//     use bio::io::fasta;

//     #[test]
//     fn test_query_passes() {
//         let kmer_size = 3;
//         // All kmers must match
//         let all_threshold = 1.0;
//         // No match required
//         let no_threshold = 0.0;

//         let genome =
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test1", None, "ATCGCA".as_ref()));
//         let read_same =
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test2", None, "ATCG".as_ref()));
//         let read_different =
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test3", None, "AAAA".as_ref()));

//         let tree = bloom_tree::create_bloom_tree(vec![genome], &kmer_size);
//         let root = tree.root.unwrap();

//         assert!(query_passes(&root, &read_same, all_threshold, kmer_size));
//         assert!(!query_passes(
//             &root,
//             &read_different,
//             all_threshold,
//             kmer_size,
//         ));
//         assert!(query_passes(&root, &read_same, no_threshold, kmer_size));
//         assert!(query_passes(
//             &root,
//             &read_different,
//             no_threshold,
//             kmer_size,
//         ));
//     }

//     #[test]
//     fn test_query_and_leaf_counts() {
//         let kmer_size = 5;
//         // Only 10% of reads must hit to count as a mapped read
//         let threshold = 0.1;

//         let records = [
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "baseline",
//                 None,
//                 "ATCAG".as_ref(),
//             )),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("diff", None, "TTTAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "onediff_first",
//                 None,
//                 "CTCAG".as_ref(),
//             )),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "onediff_mid",
//                 None,
//                 "ATTAG".as_ref(),
//             )),
//         ];

//         let mut tree = bloom_tree::create_bloom_tree(records.into_iter().collect(), &kmer_size);

//         let read = [RecordTypes::FastaRecord(fasta::Record::with_attrs(
//             "baseline",
//             None,
//             "ATCAG".as_ref(),
//         ))]
//         .into_iter()
//         .collect();
//         tree = query_batch(tree, read, threshold);

//         let mut expected_counts = [
//             ("baseline", 1),
//             ("diff", 0),
//             ("onediff_first", 0),
//             ("onediff_mid", 0),
//         ];

//         assert_eq!(
//             get_leaf_counts(&tree.root.unwrap()).sort(),
//             expected_counts.sort()
//         );
//     }

//     #[test]
//     fn test_query_and_leaf_counts_smaller_kmer_size() {
//         let kmer_size = 4;
//         // Only 10% of reads must hit to count as a mapped read
//         let threshold = 0.1;

//         let records = [
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "baseline",
//                 None,
//                 "ATCAG".as_ref(),
//             )),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("diff", None, "TTTAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "onediff_first",
//                 None,
//                 "CTCAG".as_ref(),
//             )),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "onediff_mid",
//                 None,
//                 "ATTAG".as_ref(),
//             )),
//         ];

//         let mut tree = bloom_tree::create_bloom_tree(records.into_iter().collect(), &kmer_size);

//         let read = [RecordTypes::FastaRecord(fasta::Record::with_attrs(
//             "baseline",
//             None,
//             "TCAG".as_ref(),
//         ))]
//         .into_iter()
//         .collect();
//         tree = query_batch(tree, read, threshold);

//         let mut expected_counts = [
//             ("baseline", 1),
//             ("diff", 0),
//             ("onediff_first", 1),
//             ("onediff_mid", 0),
//         ];

//         assert_eq!(
//             get_leaf_counts(&tree.root.unwrap()).sort(),
//             expected_counts.sort()
//         );
//     }

//     #[test]
//     fn test_query_and_leaf_counts_multiple_reads() {
//         let kmer_size = 4;
//         // now 51% of reads must hit to count as a mapped read
//         let threshold = 0.51;

//         let records = [
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "baseline",
//                 None,
//                 "ATCAG".as_ref(),
//             )),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("diff", None, "TTTAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "onediff_first",
//                 None,
//                 "CTCAG".as_ref(),
//             )),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "onediff_mid",
//                 None,
//                 "ATTAG".as_ref(),
//             )),
//         ];

//         let mut tree = bloom_tree::create_bloom_tree(records.into_iter().collect(), &kmer_size);

//         let read_set = [
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("baseline", None, "TCAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("baseline", None, "ATCA".as_ref())),
//         ]
//         .into_iter()
//         .collect();
//         tree = query_batch(tree, read_set, threshold);

//         // Even though there are multiple reads, since we only query once, the max number of hits
//         // is 1. Since we set the threshold above 50%, the genome that only one read maps to is not
//         // counted as hit.
//         let mut expected_counts = [
//             ("baseline", 1),
//             ("diff", 0),
//             ("onediff_first", 0),
//             ("onediff_mid", 0),
//         ];

//         assert_eq!(
//             get_leaf_counts(&tree.root.unwrap()).sort(),
//             expected_counts.sort()
//         );
//     }

//     #[test]
//     fn test_query_and_leaf_counts_multiple_queries() {
//         let kmer_size = 4;
//         // Only 10% of reads must hit to count as a mapped read
//         let threshold = 0.1;

//         let records = [
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "baseline",
//                 None,
//                 "ATCAG".as_ref(),
//             )),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("diff", None, "TTTAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "onediff_first",
//                 None,
//                 "CTCAG".as_ref(),
//             )),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "onediff_mid",
//                 None,
//                 "ATTAG".as_ref(),
//             )),
//         ];

//         let mut tree = bloom_tree::create_bloom_tree(records.into_iter().collect(), &kmer_size);

//         let all_read_sets = [
//             [RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "baseline",
//                 None,
//                 "TCAG".as_ref(),
//             ))]
//             .into_iter()
//             .collect(),
//             [RecordTypes::FastaRecord(fasta::Record::with_attrs(
//                 "baseline",
//                 None,
//                 "ATCA".as_ref(),
//             ))]
//             .into_iter()
//             .collect(),
//         ];

//         for read_set in all_read_sets {
//             tree = query_batch(tree, read_set, threshold);
//         }

//         // The first read set should map to both baseline and onediff_first, the second should only map to baseline.
//         let mut expected_counts = [
//             ("baseline", 2),
//             ("diff", 0),
//             ("onediff_first", 1),
//             ("onediff_mid", 0),
//         ];

//         assert_eq!(
//             get_leaf_counts(&tree.root.unwrap()).sort(),
//             expected_counts.sort()
//         );
//     }
// }
