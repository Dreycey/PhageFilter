use crate::bloom_filter::{get_bloom_filter, BloomFilter, DistanceChecker, ASMS};
use crate::file_parser;
use std::collections::hash_map::RandomState;
use std::rc::Weak;

struct BloomTree<R = RandomState, S = RandomState> {
    root: Option<Box<BloomNode>>,

    // Size of kmers in the bloom filters
    kmer_size: usize,

    // Random states to seed hash functions from. If we use the same random states to seed hash functions, we'll get the same exact hash functions. We need this property since we'll be unioning bloom trees.
    hash_states: (R, S),
}

struct BloomNode {
    left_child: Option<Box<BloomNode>>,
    right_child: Option<Box<BloomNode>>,
    parent: Option<Weak<BloomNode>>,
    bloom_filter: BloomFilter,

    /// NCBI Taxonomy ID for the genome
    tax_id: Option<String>,

    /// Count of reads that have mapped to this node
    mapped_reads: u32,
}

impl BloomTree<RandomState, RandomState> {
    fn new(kmer_size: usize) -> Self {
        BloomTree {
            root: None,
            kmer_size,
            // initializes random hash state to use for the whole tree
            hash_states: (RandomState::new(), RandomState::new()),
        }
    }

    fn insert(mut self, genome: &file_parser::RecordTypes) {
        // Create new leaf.
        let node: Box<BloomNode> = self.init_leaf(genome);

        // insert the node into the tree.
        match self.root {
            None => self.root = Some(node),
            Some(mut r) => self.root = BloomTree::add_to_tree(&mut r, node),
        }

        // println!("NEW GENOME: {}; length: {:#?}", id, sequence.len());
    }

    fn init_leaf(&self, genome: &file_parser::RecordTypes) -> Box<BloomNode> {
        // parse information from the input file
        let id: &str = file_parser::get_id(&genome);
        let sequence: Vec<u8> = file_parser::get_sequence(&genome);
        // ATGC -> AT, TG, GC
        let kmers = file_parser::get_kmers(&sequence, &self.kmer_size);
        // Create new node.
        let mut node: Box<BloomNode> = Box::new(BloomNode::new(
            self.hash_states.clone(),
            Some(id.to_string()),
        ));
        // map kmers into the bloom filter.
        for kmer in kmers {
            node.bloom_filter.insert(&kmer);
        }
        return node;
    }

    fn init_internal_node(
        &self,
        current_node: Box<BloomNode>,
        new_node: Box<BloomNode>,
    ) -> Box<BloomNode> {
        // Create new node.
        let mut node: Box<BloomNode> = Box::new(BloomNode::new(self.hash_states.clone(), None));

        // Combine children's bloom filters
        node.bloom_filter.union(&current_node.bloom_filter);
        node.bloom_filter.union(&new_node.bloom_filter);

        // assigns left child to current node, and new node to right child.
        node.left_child = Some(current_node);
        node.right_child = Some(new_node);

        return node;
    }

    fn add_to_tree(
        current_node: &mut Box<BloomNode>,
        node: Box<BloomNode>,
    ) -> Option<Box<BloomNode>> {
        // if left or right child of current node is empty, add node.
        match (&mut current_node.left_child, &mut current_node.right_child) {
            (None, Some(_)) => {
                // Update future ancestor to include this bloom filter
                current_node.bloom_filter.union(&node.bloom_filter);

                current_node.left_child = Some(node);
            }
            (Some(_), None) => {
                // Update future ancestor to include this bloom filter
                current_node.bloom_filter.union(&node.bloom_filter);

                current_node.right_child = Some(node);
            }
            (Some(left), Some(right)) => {
                // Update future ancestor to include this bloom filter
                current_node.bloom_filter.union(&node.bloom_filter);

                // choose the closest child node - checks hamming distance for each child.
                let right_distance = right.bloom_filter.distance(&node.bloom_filter);
                let left_distance = left.bloom_filter.distance(&node.bloom_filter);
                // Add to the subtree with the closest distance so we can minimize false-positives in the bloom filters.
                if right_distance < left_distance {
                    current_node.right_child = BloomTree::add_to_tree(right, node);
                } else {
                    current_node.left_child = BloomTree::add_to_tree(left, node);
                }
            }
            // Currently at leaf node, need to merge.
            (None, None) => {
                // Create new 'unioned' node
                // update the 'unioned' node BF with child nodes BFs
                // update left and right child of 'unioned' node
                // update parent of existing node to point to 'unioned' node.
            }
        }
        return Some(*current_node);
    }
}

impl BloomNode {
    fn new(hash_states: (RandomState, RandomState), tax_id: Option<String>) -> Self {
        BloomNode {
            // initializes random hash state to use for the whole tree
            left_child: None,
            right_child: None,
            parent: None,
            tax_id,
            bloom_filter: get_bloom_filter(hash_states),
            mapped_reads: 0,
        }
    }
}

pub fn add_to_bloom(
    parsed_genomes: Vec<file_parser::RecordTypes>,
    kmer_size: &usize,
    bloom_filter: &mut BloomFilter,
) {
    for genome in parsed_genomes {
        let sequence: Vec<u8> = file_parser::get_sequence(&genome);
        let id: &str = file_parser::get_id(&genome);
        let kmers = file_parser::get_kmers(&sequence, &kmer_size); // ATGC -> AT, TG, GC
        for kmer in kmers {
            bloom_filter.insert(&kmer);
        }
        println!("NEW GENOME: {}; length: {:#?}", id, sequence.len());
    }
}
