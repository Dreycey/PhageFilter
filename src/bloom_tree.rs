use crate::bloom_filter::hasher::HashSeed;
/// Methods and definitions for a BloomTree and BloomNodes
/// within the BloomTree. These methods allow for creating
/// a bloom tree for a set of given genomes.
///
///
/// # Example Usage
///
/// ```rust
/// let mut bloom_node = bloom_tree::create_bloom_tree(parsed_genomes, &kmer_size);
/// ```
use crate::bloom_filter::{create_bloom_filter, BloomFilter, DistanceChecker, ASMS};
use crate::Arc;
use crate::{
    cache::{BFLruCache, BloomFilterCache},
    file_parser,
};
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::RwLock;

const TREE_FILENAME: &str = "tree.bin";

#[derive(Debug, Deserialize, Serialize)]
pub(crate) struct BloomTree<R = HashSeed, S = HashSeed> {
    pub(crate) root: Option<Box<BloomNode>>,

    #[serde(skip, default = "create_cache")]
    pub(crate) bf_cache: Box<dyn BloomFilterCache>,

    /// Path to the directory containing the tree and corresponding bloom filters
    #[serde(skip)]
    directory: Option<PathBuf>,

    /// bloom filter params.
    false_pos_rate: f32,
    largest_expected_genome: u32,

    /// Size of kmers in the bloom filters
    pub(crate) kmer_size: usize,

    /// Random states to seed hash functions from. If we use the same random states to seed hash functions, we'll get the same exact hash functions. We need this property since we'll be unioning bloom trees.
    pub(crate) hash_states: (R, S),
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub(crate) struct BloomNode {
    pub(crate) left_child: Option<Box<BloomNode>>,
    pub(crate) right_child: Option<Box<BloomNode>>,
    pub(crate) bloom_filter_path: PathBuf,

    /// NCBI Taxonomy ID for the genome
    pub(crate) tax_id: Option<String>,

    /// Count of reads that have mapped to this node
    pub(crate) mapped_reads: usize,
}

fn create_cache() -> Box<BFLruCache> {
    // TODO: add the ability to change the size of the cache.
    Box::new(BFLruCache::new(1, PathBuf::from("./")))
}

/// RandomState doesn't support equality comparisons so we ignore the hash_states when comparing BloomTrees.
impl PartialEq for BloomTree<HashSeed, HashSeed> {
    /// Check if 2 BloomTrees are equal
    ///
    /// # Returns
    /// - a bool of whether the given BloomTrees are equal.
    fn eq(&self, other: &Self) -> bool {
        self.root == other.root && self.kmer_size == other.kmer_size
    }
}

impl PartialEq for BloomNode {
    /// Check if 2 BloomNodes are equal
    ///
    /// # Returns
    /// - a bool of whether the given BloomNodes are equal.
    fn eq(&self, other: &Self) -> bool {
        self.left_child == other.left_child
            && self.right_child == other.right_child
            && self.tax_id == other.tax_id
            && self.mapped_reads == other.mapped_reads
    }
}

impl BloomTree<HashSeed, HashSeed> {
    /// Construct a new BloomTree
    ///
    /// # Arguments
    /// - `kmer_size`: The kmer size the tree will use for both building and queries.
    ///
    /// # Returns
    /// - a BloomTree instance.
    pub fn new(
        kmer_size: usize,
        directory: &PathBuf,
        bf_cache: BFLruCache,
        false_pos_rate: f32,
        largest_expected_genome: u32,
    ) -> Self {
        std::fs::create_dir_all(directory).unwrap();
        BloomTree {
            root: None,
            bf_cache: Box::new(bf_cache),
            false_pos_rate,
            kmer_size,
            largest_expected_genome,
            // initializes random hash state to use for the whole tree
            hash_states: (HashSeed::new(), HashSeed::new()),
            directory: Some(directory.to_path_buf()),
        }
    }

    /// A wrapper function to add either create a root node,
    /// or add a new genome to the tree.
    ///
    /// # Arguments
    /// - `genome`: a reference to the genome (file_parser::RecordTypes)
    ///
    /// # Returns
    /// - self, a BloomTree instance.
    pub fn insert(&mut self, genome: &file_parser::DNASequence) {
        log::info!(
            "(bloom tree; insert()) inserting {} into the tree",
            genome.id
        );

        // Create new leaf.
        let node: Box<BloomNode> = self.init_leaf_node(genome);

        // insert the node into the tree.
        let root = self.root.take();
        match root {
            None => self.root = Some(node),
            Some(r) => self.root = Some(self.add_to_tree(r, node)),
        }
    }

    /// Initializes a leaf node with the given genome. Returns
    /// a new BloomNode with the generated bloom filter for the
    /// genomes kmers.
    ///
    /// # Arguments
    /// - `genome`: a reference to the genome (file_parser::RecordTypes)
    ///
    /// # Returns
    /// - A new BloomNode, representing a to-be leaf node.
    fn init_leaf_node(&self, genome: &file_parser::DNASequence) -> Box<BloomNode> {
        let new_leaf_node: Box<BloomNode> = self.make_bloom_node(genome.id.clone());
        // map kmers into the bloom filter.
        {
            let b4_new_leaf_node = self.get_bf(&new_leaf_node);
            let mut bf_new_leaf_node = b4_new_leaf_node.write().unwrap();

            for kmer in &genome.kmers {
                bf_new_leaf_node.insert(&kmer);
            }
        }

        log::debug!("(bloom tree; init()) Adding {}", genome.id);
        new_leaf_node
    }

    /// Adds a BloomNode to a given BloomTree.
    ///
    /// # Description
    /// This method traverses a BloomTree, adding a BloomNode where best fit
    /// in a greedy manner. This works by checking if the current node has children.
    /// If it has 2 children, then it will choose which one to check next based on the
    /// minimum hamming distance. If it has 1 child, it add the new node to the current node.
    /// If it's a leaf, then it create a new internal node, with the children being the new node
    /// and the current node.
    ///
    /// # Arguments
    /// - `current_node`: The current BloomNode being evaluated (used recursively).
    /// - `node`: The BloomNode being added to the tree.
    /// - `hash_states`: randomized hash states made using bloom_filter::hasher::HashSeed
    ///
    /// # Returns
    /// - The current node after modification to self or children.
    fn add_to_tree(
        &self,
        mut current_node: Box<BloomNode>,
        node: Box<BloomNode>,
    ) -> Box<BloomNode> {
        if let (Some(left), Some(right)) = (&current_node.left_child, &current_node.right_child) {
            // Get bloom filters from the cache
            self.node_union(&current_node, &node);

            // Calculate distances between the new node's bloom filter and left/right children
            let right_distance = self.node_distance(right, &node);
            let left_distance = self.node_distance(left, &node);

            if right_distance < left_distance {
                current_node.right_child =
                    Some(self.add_to_tree(current_node.right_child.unwrap(), node));
            } else {
                current_node.left_child =
                    Some(self.add_to_tree(current_node.left_child.unwrap(), node));
            }
        } else if current_node.left_child.is_none() && current_node.right_child.is_none() {
            current_node = self.init_internal_node(current_node, node);
        } else {
            panic!("Node with only one child encountered - should not happen.");
        }

        current_node
    }

    /// Returns an internal node, given two nodes containing
    /// genomes (i.e. to-be leaf nodes).
    ///
    /// # Arguments
    /// - `current_node`: BloomNode representating a leaf node already in the tree.
    /// - `new_node`: The BloomNode being added to the tree.
    /// - `hash_states`: randomized hash states made using bloom_filter::hasher::HashSeed
    ///
    /// # Returns
    /// - The new internal node, with children being the given leaf nodes.
    fn init_internal_node(
        &self,
        current_node: Box<BloomNode>,
        new_node: Box<BloomNode>,
    ) -> Box<BloomNode> {
        // create a new internal node.
        let mut rng = rand::thread_rng();
        let n2: u16 = rng.gen();
        let node_name = "Internal_Node".to_string() + "_" + &n2.to_string();
        let mut new_internal_node = self.make_bloom_node(node_name.clone());

        // union children nodes into new parent internal node
        self.node_union(&new_internal_node, &new_node);
        self.node_union(&new_internal_node, &current_node);

        // assign children nodes to the new internaal node.
        new_internal_node.left_child = Some(current_node);
        new_internal_node.right_child = Some(new_node);

        new_internal_node
    }

    fn node_distance(&self, node_a: &BloomNode, node_b: &BloomNode) -> usize {
        // get BF for node a
        let b4_node_a = self.get_bf(node_a);
        let bf_node_a = b4_node_a.read().unwrap();

        // get BF for node b
        let b4_node_b = self.get_bf(node_b);
        let bf_node_b = b4_node_b.read().unwrap();

        bf_node_a.distance(&bf_node_b)
    }

    fn node_union(&self, node2modify: &BloomNode, node2add: &BloomNode) {
        // get BF for node to modify.
        let b4_node2modify = self.get_bf(node2modify);
        let mut bf_node2modify = b4_node2modify.write().unwrap();

        // get BF for node beiing added.
        let b4_node2add = self.get_bf(node2add);
        let bf_node2add = b4_node2add.read().unwrap();

        // union BFs.
        bf_node2modify.union(&bf_node2add);
    }

    fn get_bf(&self, bloom_node: &BloomNode) -> Arc<RwLock<BloomFilter>> {
        self.bf_cache
            .get_filter(&bloom_node.bloom_filter_path)
            .expect("BF was not found!")
    }

    fn make_bloom_node(&self, nodeid: String) -> Box<BloomNode> {
        // bloom filter path.
        let bloom_filter_name = PathBuf::from(nodeid.clone() + ".bf");
        let full_bloom_filter_path = self.directory.clone().unwrap().join(&bloom_filter_name);

        // create and save the bloom filter.
        let bloom_filter = create_bloom_filter(
            self.hash_states.clone(),
            full_bloom_filter_path.clone(),
            self.false_pos_rate,
            self.largest_expected_genome,
        );

        // Cache the loaded Bloom filter
        // not needed; saves time adding to cache here! (Minimizes needless I/O)
        self.bf_cache.add_filter(&bloom_filter_name, bloom_filter);

        // make node.
        let bloomnode = BloomNode::new(Some(nodeid.clone()), bloom_filter_name);
        Box::new(bloomnode)
    }

    /// - if the directory does not exist.
    pub fn prune_tree(&mut self, search_depth: usize) {
        // print where the tree is being saved.
        log::info!(
            " (bloom tree; prune_tree()) pruning tree to depth of: {}",
            search_depth
        );
        // set up queue
        let mut queue: Vec<(&mut Box<BloomNode>, usize)> = vec![];
        let root = self.root.as_mut().unwrap();
        queue.push((root, 0));

        // BFS through tree
        while !queue.is_empty() {
            if let Some((node, node_depth)) = queue.pop() {
                if node_depth < search_depth {
                    let node_depth = node_depth + 1;
                    if node.left_child.is_some() {
                        queue.push((node.left_child.as_mut().unwrap(), node_depth));
                    }
                    if node.right_child.is_some() {
                        queue.push((node.right_child.as_mut().unwrap(), node_depth));
                    }
                } else {
                    node.left_child = None;
                    node.right_child = None;
                }
            }
        }
    }

    /// This method saves the tree to disk using a given directory name.
    ///
    /// # Arguments
    /// - `directory`: Directory name where the serialized tree will be stored.
    ///
    /// # Panics
    /// - if the directory does not exist.
    pub fn save(&self, directory: &Path) {
        // Create parent directories if they don't already exist
        std::fs::create_dir_all(directory).unwrap();
        // panics if it doesn't exist
        if !directory.is_dir() {
            panic!("Must provide a directory in which to store the tree");
        }
        // print where the tree is being saved.
        log::info!(
            " (bloom tree; save()) Saving bloom tree to directory {}",
            directory.canonicalize().unwrap().to_str().unwrap()
        );
        // serialize the tree
        let mut tree_file = BufWriter::new(File::create(directory.join(TREE_FILENAME)).unwrap());
        bincode::serialize_into(&mut tree_file, self).unwrap();
        tree_file.flush().unwrap();
    }

    /// This method loads a serialized tree from disk into memory.
    ///
    /// # Arguments
    /// - `directory`: Directory name where the serialized tree is located.
    ///
    /// # Panics
    /// - if the directory does not exist.
    pub fn load(directory: &Path, bf_cache: BFLruCache) -> BloomTree {
        // panics if it doesn't exist
        if !directory.is_dir() {
            panic!("Must provide a directory in where a tree has been stored");
        }
        // print where the tree is being saved.
        log::info!(
            "(bloom tree; load()) Reading bloom tree from directory {}",
            directory.canonicalize().unwrap().to_str().unwrap()
        );
        // serialized bloom tree path
        let tree_file: File = File::open(directory.join(TREE_FILENAME)).unwrap();
        // create a buffered reader for the file
        let tree_reader = BufReader::new(tree_file);
        // deserialize the tree from the buffered reader
        let mut deserialized_tree: BloomTree = bincode::deserialize_from(tree_reader).unwrap();
        // save path to the directory containing the serialized tree.
        deserialized_tree.directory = Some(directory.to_path_buf());
        // update cache information
        deserialized_tree.bf_cache = Box::new(bf_cache);

        deserialized_tree
    }
}

impl BloomNode {
    /// Creates a new instance of `BloomNode`.
    ///
    /// # Arguments
    ///
    /// * `hash_states` - The randomized hash states made using `bloom_filter::hasher::HashSeed`.
    /// * `tax_id` - The taxonomic identifier for the genome, or the NCBI accession.
    ///
    /// # Returns
    ///
    /// A new instance of `BloomNode`.
    fn new(tax_id: Option<String>, bloom_filter_path: PathBuf) -> Self {
        BloomNode {
            // initializes random hash state to use for the whole tree
            left_child: None,
            right_child: None,
            tax_id,
            bloom_filter_path,
            mapped_reads: 0,
        }
    }

    /// Returns a boolean indicating whether the node is a leaf node.
    ///
    /// # Returns
    ///
    /// `true` if the node is a leaf node; `false` otherwise.
    pub(crate) fn is_leafnode(&self) -> bool {
        self.left_child.is_none() && self.right_child.is_none()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bloom_filter::BloomFilter;
    use crate::cache;
    use std::fs;

    struct TestContext {
        dir: PathBuf,
    }

    impl TestContext {
        fn setup(name: &str) -> TestContext {
            let dir = PathBuf::from(format!("tmp_testing_{}/", name));
            if dir.is_dir() {
                let _ = fs::remove_dir_all(&dir);
            }
            fs::create_dir_all(&dir).unwrap();
            TestContext { dir }
        }

        fn get_tmp_dir(&self) -> PathBuf {
            self.dir.clone()
        }

        fn teardown(&self) {
            let _ = fs::remove_dir_all(&self.dir);
        }
    }

    #[allow(dead_code)]
    fn create_bloom_filter(false_positive_rate: f32, largest_expected_genome: u32) -> BloomFilter {
        let hash_states = (HashSeed::new(), HashSeed::new());
        BloomFilter::with_rate(false_positive_rate, largest_expected_genome, hash_states)
    }

    #[test]
    fn test_empty_tree_insert() {
        // setup context
        let context = TestContext::setup("bt_empty");

        // initialization states.
        let kmer_size = 5;
        let directory = context.get_tmp_dir();
        let cache_size = 5;
        let false_positive_rate = 0.001;
        let largest_expected_genome = 1000;
        assert!(directory.exists());

        // create a DNA sequence.
        let record_id = "tmp_bloom";
        let kmers = file_parser::get_kmers("ATCAG".as_bytes(), &kmer_size);
        let record = file_parser::DNASequence::new(
            Some("ATCAG".as_bytes().to_vec()),
            None,
            record_id.to_string(),
            kmers,
        );

        // create a new cache
        let bloomfilter_cache = cache::BFLruCache::new(cache_size, directory.clone());

        // (automated) creating tree.
        let mut tree = BloomTree::new(
            kmer_size,
            &directory,
            bloomfilter_cache,
            false_positive_rate,
            largest_expected_genome,
        );
        tree.insert(&record);

        //(manual) create a bloom filter with the genome.
        let mut expected_bf = BloomFilter::with_rate(
            false_positive_rate,
            largest_expected_genome,
            tree.hash_states.clone(),
        );
        for kmer in record.kmers {
            expected_bf.insert(&kmer);
        }

        //(manual) create corresponding BloomNode
        let expected_root = BloomNode::new(Some(record_id.to_string()), directory.join(record_id));

        // create a BloomTree manually.
        let expected_tree = BloomTree {
            root: Some(Box::new(expected_root)),
            directory: Some(directory.clone()),
            kmer_size,
            hash_states: tree.hash_states.clone(),
            bf_cache: Box::new(cache::BFLruCache::new(cache_size, directory.clone())),
            false_pos_rate: false_positive_rate,
            largest_expected_genome,
        };

        assert_eq!(tree, expected_tree);

        drop(tree);
        context.teardown();
    }

    #[test]
    fn test_one_elem_tree_insert() {
        let context = TestContext::setup("bt_one_elem");
        let kmer_size = 5;
        let directory = context.get_tmp_dir();
        let cache_size = 5;
        let false_positive_rate = 0.001;
        let largest_expected_genome = 1000;

        let kmers1 = file_parser::get_kmers("ATCAG".as_bytes(), &kmer_size);
        let record1 = file_parser::DNASequence::new(
            Some("ATCAG".as_bytes().to_vec()),
            "test1".to_string(),
            kmers1,
        );
        let kmers2 = file_parser::get_kmers("TTTAG".as_bytes(), &kmer_size);
        let record2 = file_parser::DNASequence::new(
            Some("TTTAG".as_bytes().to_vec()),
            "test2".to_string(),
            kmers2,
        );

        let bloomfilter_cache = cache::BFLruCache::new(cache_size, directory.clone());
        let mut tree = BloomTree::new(
            kmer_size,
            &directory,
            bloomfilter_cache,
            false_positive_rate,
            largest_expected_genome,
        );
        tree.insert(&record1);
        tree.insert(&record2);

        let root = tree.root.as_ref().expect("Root should exist");
        assert!(root.left_child.is_some(), "Root should have a left child");
        assert!(
            root.right_child.is_some(),
            "Root should have a right child"
        );

        let left = root.left_child.as_ref().unwrap();
        let right = root.right_child.as_ref().unwrap();

        // Both children should be leaves
        assert!(
            left.left_child.is_none() && left.right_child.is_none(),
            "Left child should be a leaf"
        );
        assert!(
            right.left_child.is_none() && right.right_child.is_none(),
            "Right child should be a leaf"
        );

        // Tax IDs of leaves match the inserted genomes
        assert_eq!(left.tax_id.as_deref(), Some("test1"));
        assert_eq!(right.tax_id.as_deref(), Some("test2"));

        drop(tree);
        context.teardown();
    }

    #[test]
    fn test_nested_tree_insert_left() {
        let context = TestContext::setup("bt_nested_left");
        let kmer_size = 5;
        let directory = context.get_tmp_dir();
        let cache_size = 5;
        let false_positive_rate = 0.001;
        let largest_expected_genome = 1000;

        // genome3 has same sequence as genome1 — should be closer in hamming distance
        let sequences = [("test1", "ATCAG"), ("test2", "TTTAG"), ("test3", "ATCAG")];
        let records: Vec<_> = sequences
            .iter()
            .map(|(id, seq)| {
                let kmers = file_parser::get_kmers(seq.as_bytes(), &kmer_size);
                file_parser::DNASequence::new(
                    Some(seq.as_bytes().to_vec()),
                    id.to_string(),
                    kmers,
                )
            })
            .collect();

        let bloomfilter_cache = cache::BFLruCache::new(cache_size, directory.clone());
        let mut tree = BloomTree::new(
            kmer_size,
            &directory,
            bloomfilter_cache,
            false_positive_rate,
            largest_expected_genome,
        );
        for record in &records {
            tree.insert(record);
        }

        // Expected tree:
        //                    root (internal)
        //                   /              \
        //         left (internal)     right (leaf: test2)
        //          /          \
        //   (leaf: test1)  (leaf: test3)
        let root = tree.root.as_ref().expect("Root should exist");
        assert!(root.left_child.is_some() && root.right_child.is_some());

        // Right child should be a leaf with test2's tax_id
        let right = root.right_child.as_ref().unwrap();
        assert!(
            right.left_child.is_none() && right.right_child.is_none(),
            "Right child should be a leaf"
        );
        assert_eq!(right.tax_id.as_deref(), Some("test2"));

        // Left child should be an internal node with two leaf children
        let left = root.left_child.as_ref().unwrap();
        assert!(
            left.left_child.is_some() && left.right_child.is_some(),
            "Left child should be an internal node"
        );

        let ll = left.left_child.as_ref().unwrap();
        let lr = left.right_child.as_ref().unwrap();
        assert!(ll.left_child.is_none() && ll.right_child.is_none());
        assert!(lr.left_child.is_none() && lr.right_child.is_none());

        let inner_ids: Vec<_> = [ll.tax_id.as_deref(), lr.tax_id.as_deref()]
            .into_iter()
            .collect();
        assert!(inner_ids.contains(&Some("test1")));
        assert!(inner_ids.contains(&Some("test3")));

        drop(tree);
        context.teardown();
    }

    #[test]
    fn test_nested_tree_insert_right() {
        let context = TestContext::setup("bt_nested_right");
        let kmer_size = 5;
        let directory = context.get_tmp_dir();
        let cache_size = 5;
        let false_positive_rate = 0.001;
        let largest_expected_genome = 1000;

        // genome3 has same sequence as genome2 — should be closer in hamming distance
        let sequences = [("test1", "ATCAG"), ("test2", "TTTAG"), ("test3", "TTTAG")];
        let records: Vec<_> = sequences
            .iter()
            .map(|(id, seq)| {
                let kmers = file_parser::get_kmers(seq.as_bytes(), &kmer_size);
                file_parser::DNASequence::new(
                    Some(seq.as_bytes().to_vec()),
                    id.to_string(),
                    kmers,
                )
            })
            .collect();

        let bloomfilter_cache = cache::BFLruCache::new(cache_size, directory.clone());
        let mut tree = BloomTree::new(
            kmer_size,
            &directory,
            bloomfilter_cache,
            false_positive_rate,
            largest_expected_genome,
        );
        for record in &records {
            tree.insert(record);
        }

        // Expected tree:
        //            root (internal)
        //           /              \
        //  left (leaf: test1)   right (internal)
        //                        /          \
        //                 (leaf: test2)  (leaf: test3)
        let root = tree.root.as_ref().expect("Root should exist");
        assert!(root.left_child.is_some() && root.right_child.is_some());

        // Left child should be a leaf with test1's tax_id
        let left = root.left_child.as_ref().unwrap();
        assert!(
            left.left_child.is_none() && left.right_child.is_none(),
            "Left child should be a leaf"
        );
        assert_eq!(left.tax_id.as_deref(), Some("test1"));

        // Right child should be an internal node with two leaf children
        let right = root.right_child.as_ref().unwrap();
        assert!(
            right.left_child.is_some() && right.right_child.is_some(),
            "Right child should be an internal node"
        );

        let rl = right.left_child.as_ref().unwrap();
        let rr = right.right_child.as_ref().unwrap();
        assert!(rl.left_child.is_none() && rl.right_child.is_none());
        assert!(rr.left_child.is_none() && rr.right_child.is_none());

        let inner_ids: Vec<_> = [rl.tax_id.as_deref(), rr.tax_id.as_deref()]
            .into_iter()
            .collect();
        assert!(inner_ids.contains(&Some("test2")));
        assert!(inner_ids.contains(&Some("test3")));

        drop(tree);
        context.teardown();
    }

    #[test]
    fn test_save_load() {
        let context = TestContext::setup("bt_save_load");
        let kmer_size = 4;
        let scratch_dir = context.get_tmp_dir();
        let cache_size = 5;
        let false_positive_rate = 0.001;
        let largest_expected_genome = 1000;

        let sequences = [("test1", "ATCAG"), ("test2", "TTTAG"), ("test3", "TTTAG")];
        let records: Vec<_> = sequences
            .iter()
            .map(|(id, seq)| {
                let kmers = file_parser::get_kmers(seq.as_bytes(), &kmer_size);
                file_parser::DNASequence::new(
                    Some(seq.as_bytes().to_vec()),
                    id.to_string(),
                    kmers,
                )
            })
            .collect();

        let bloomfilter_cache = cache::BFLruCache::new(cache_size, scratch_dir.clone());
        let mut tree = BloomTree::new(
            kmer_size,
            &scratch_dir,
            bloomfilter_cache,
            false_positive_rate,
            largest_expected_genome,
        );
        for record in &records {
            tree.insert(record);
        }

        tree.save(&scratch_dir);

        let load_cache = cache::BFLruCache::new(cache_size, scratch_dir.clone());
        let loaded_tree = BloomTree::load(&scratch_dir, load_cache);

        assert_eq!(tree, loaded_tree);

        drop(loaded_tree);
        drop(tree);
        context.teardown();
    }
}
