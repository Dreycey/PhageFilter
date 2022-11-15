use crate::bloom_filter::{get_bloom_filter, BloomFilter, DistanceChecker, ASMS};
use crate::file_parser;
use rand::Rng;
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
    // parent: Option<Weak<BloomNode>>,
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

    fn insert(mut self, genome: &file_parser::RecordTypes) -> Self {
        // Create new leaf.
        let node: Box<BloomNode> = self.init_leaf(genome);

        // insert the node into the tree.
        match self.root {
            None => self.root = Some(node),
            Some(r) => self.root = Some(BloomTree::add_to_tree(r, node, &self.hash_states)),
        }

        // println!("NEW GENOME: {}; length: {:#?}", id, sequence.len());
        return self;
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
        println!("\tAdding {}", id);
        return node;
    }

    fn init_internal_node(
        current_node: Box<BloomNode>,
        new_node: Box<BloomNode>,
        hash_states: &(RandomState, RandomState),
    ) -> Box<BloomNode> {
        // get random number for internal node
        let mut rng = rand::thread_rng();
        let n2: u16 = rng.gen();
        // create internal node ID using random number generator
        let mut node_name = "Internal Node".to_string();
        node_name = format!("{}_{}", node_name, n2.to_string());
        // Create new node.
        let mut node: Box<BloomNode> =
            Box::new(BloomNode::new(hash_states.clone(), Some(node_name)));

        // Combine children's bloom filters
        node.bloom_filter.union(&current_node.bloom_filter);
        node.bloom_filter.union(&new_node.bloom_filter);

        // assigns left child to current node, and new node to right child.
        node.left_child = Some(current_node);
        node.right_child = Some(new_node);

        return node;
    }

    fn add_to_tree(
        mut current_node: Box<BloomNode>,
        node: Box<BloomNode>,
        hash_states: &(RandomState, RandomState),
    ) -> Box<BloomNode> {
        // TODO: delete debug:
        println!("Looking at current node: {:?}", current_node.tax_id);
        // if left or right child of current node is empty, add node.
        match (&current_node.left_child, &current_node.right_child) {
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
                    current_node.right_child = Some(BloomTree::add_to_tree(
                        current_node.right_child.unwrap(),
                        node,
                        hash_states,
                    ));
                } else {
                    current_node.left_child = Some(BloomTree::add_to_tree(
                        current_node.left_child.unwrap(),
                        node,
                        hash_states,
                    ));
                }
            }
            // Currently at leaf node, need to merge.
            (None, None) => {
                // Create new 'unioned' internal node
                current_node = BloomTree::init_internal_node(current_node, node, hash_states);
            }
        }
        return current_node;
    }
}

impl BloomNode {
    fn new(hash_states: (RandomState, RandomState), tax_id: Option<String>) -> Self {
        BloomNode {
            // initializes random hash state to use for the whole tree
            left_child: None,
            right_child: None,
            tax_id,
            bloom_filter: get_bloom_filter(hash_states),
            mapped_reads: 0,
        }
    }
}

pub fn create_bloom_tree(parsed_genomes: Vec<file_parser::RecordTypes>, kmer_size: &usize) {
    // create initial bloom tree
    let mut bloom_tree = BloomTree::new(*kmer_size);
    // add genomes to bloom tree
    for genome in parsed_genomes {
        bloom_tree = bloom_tree.insert(&genome);
        println!("NEW GENOME added to bloom tree");
    }
    // print final bloom tree
    println!("\n\n In-Order traversal of bloom tree: \n");
    traverse_bloom_tree(&bloom_tree.root.unwrap())
}

fn traverse_bloom_tree(bloom_node: &BloomNode) {
    // traverse left children
    match &bloom_node.left_child {
        None => println!("{:?}-left: NULL", bloom_node.tax_id),
        Some(child) => {
            println!("{:?}-left: {:?}", bloom_node.tax_id, child.tax_id);
            traverse_bloom_tree(child)
        }
    }

    // print current value.
    // println!("{:?}", bloom_node.tax_id);

    // traverse right children
    match &bloom_node.right_child {
        None => println!("{:?}-right: NULL", bloom_node.tax_id),
        Some(child) => {
            println!("{:?}-right: {:?}", bloom_node.tax_id, child.tax_id);
            traverse_bloom_tree(child)
        }
    }
}
