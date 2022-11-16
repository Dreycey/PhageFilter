use crate::bloom_filter::{get_bloom_filter, BloomFilter, DistanceChecker, ASMS};
use crate::bloom_tree::{BloomNode, BloomTree};
use crate::file_parser;
use rayon::prelude::*;

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
    let cutoff = (threshold * kmers.len() as f32).ceil() as usize;
    println!(
        "possible kmers: {:#?}, cutoff: {:#?}, number of matches: {:#?}",
        kmers.len(),
        cutoff,
        num_matches
    );
    return num_matches >= (threshold * kmers.len() as f32).ceil() as usize;
}

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
