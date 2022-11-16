use crate::bloom_filter::ASMS;
use crate::bloom_tree::{BloomNode, BloomTree};
use crate::file_parser;
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
    read: &file_parser::RecordTypes,
    threshold: f32,
    kmer_size: usize,
) -> bool {
    let sequence: Vec<u8> = file_parser::get_sequence(read);
    let kmers: Vec<Vec<u8>> = file_parser::get_kmers(&sequence, &kmer_size);
    let num_matches = kmers
        .iter()
        .filter(|kmer| bloom_node.bloom_filter.contains(kmer))
        .count();
    return num_matches >= (threshold * kmers.len() as f32).ceil() as usize;
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
    read_set: Vec<file_parser::RecordTypes>,
    threshold: f32,
) -> BloomTree {
    bloom_tree.root = match bloom_tree.root {
        Some(root) => Some(_query_batch(
            root,
            read_set,
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
/// - The BloomNode after evulating itself and all its children (used recursively).
///
/// # Panics
/// - N/A
fn _query_batch(
    mut bloom_node: Box<BloomNode>,
    read_set: Vec<file_parser::RecordTypes>,
    threshold: f32,
    kmer_size: usize,
) -> Box<BloomNode> {
    // check which reads pass the bloom filter check; in parallel
    let pass = read_set
        .par_iter()
        .filter(|read| query_passes(&bloom_node, read, threshold, kmer_size));

    // process the node differently, if intenal or leaf
    if !bloom_node.is_leafnode() {
        // add to queue for further checking children nodes.
        let queue: Vec<file_parser::RecordTypes> = pass.cloned().collect();

        // If children nodes, check these nodes against reads pass the parent bloom filter.
        if queue.len() > 0 {
            // check the left child
            if !bloom_node.left_child.is_none() {
                bloom_node.left_child = Some(_query_batch(
                    bloom_node.left_child.unwrap(),
                    queue.clone(),
                    threshold,
                    kmer_size,
                ));
            }
            // check the right child
            if !bloom_node.right_child.is_none() {
                bloom_node.right_child = Some(_query_batch(
                    bloom_node.right_child.unwrap(),
                    queue,
                    threshold,
                    kmer_size,
                ));
            }
        }
    } else {
        // count the number of reads mapped.
        bloom_node.mapped_reads = pass.count();
    }

    return bloom_node;
}

/// Recursively traverses the bloom tree, saving
/// lead node mappings to a specified output file.
///
/// # Parameters
/// - `bloom_node`: A leaf node of the type `BloomNode`
/// - `output_file`: A `File` object that may be written to.
///
/// # Returns
/// - ()
///
/// # Panics
/// - N/A
pub(crate) fn get_leaf_counts(bloom_node: &BloomNode, output_file: &mut File) {
    match (&bloom_node.left_child, &bloom_node.right_child) {
        // leaf node
        (None, None) => save_mappings(bloom_node, output_file),
        // traverse right child
        (None, Some(l_child)) => get_leaf_counts(l_child, output_file),
        // traverse left child
        (Some(r_child), None) => get_leaf_counts(r_child, output_file),
        // traverse left and right child
        (Some(r_child), Some(l_child)) => {
            get_leaf_counts(l_child, output_file);
            get_leaf_counts(r_child, output_file)
        }
    }
}

/// Saves the leaf counts to an output file.
///
/// # Parameters
/// - `bloom_node`: A leaf node of the type `BloomNode`
/// - `output_file`: A `File` object that may be written to.
///
/// # Returns
/// - ()
///
/// # Panics
/// - Panics if there is an error wrting to the output file.
fn save_mappings(bloom_node: &BloomNode, output_file: &mut File) {
    println!(
        "{},{}\n",
        bloom_node.tax_id.as_deref().unwrap(),
        bloom_node.mapped_reads
    );
    if bloom_node.mapped_reads > 0 {
        let tax_maps: String = format!(
            "{},{}\n",
            bloom_node.tax_id.as_deref().unwrap(),
            bloom_node.mapped_reads
        );
        output_file
            .write(tax_maps.as_bytes())
            .expect("problem writing to output file!");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bloom_tree;
    use crate::file_parser::RecordTypes;
    use bio::io::fasta;

    #[test]
    fn test_query_passes() {
        let kmer_size = 3;
        // All kmers must match
        let all_threshold = 1.0;
        // No match required
        let no_threshold = 0.0;

        let genome =
            RecordTypes::FastaRecord(fasta::Record::with_attrs("test1", None, "ATCGCA".as_ref()));
        let read_same =
            RecordTypes::FastaRecord(fasta::Record::with_attrs("test2", None, "ATCG".as_ref()));
        let read_different =
            RecordTypes::FastaRecord(fasta::Record::with_attrs("test3", None, "AAAA".as_ref()));

        let tree = bloom_tree::create_bloom_tree(vec![genome], &kmer_size);
        let root = tree.root.unwrap();

        assert!(query_passes(&root, &read_same, all_threshold, kmer_size));
        assert!(!query_passes(
            &root,
            &read_different,
            all_threshold,
            kmer_size
        ));
        assert!(query_passes(&root, &read_same, no_threshold, kmer_size));
        assert!(query_passes(
            &root,
            &read_different,
            no_threshold,
            kmer_size
        ));
    }
}
