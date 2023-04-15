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
use crate::bloom_filter::{BloomFilter, create_bloom_filter, DistanceChecker, ASMS};
use crate::cache::{BFLruCache, BloomFilterCache};
use crate::file_parser;
use crate::Arc;
use std::sync::RwLock;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

const TREE_FILENAME: &'static str = "tree.bin";

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
    Box::new(BFLruCache::new(1))
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
    pub fn insert(&mut self, genome: &file_parser::DNASequence) -> () {
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
            let b4_new_leaf_node = self
                .bf_cache
                .get_filter(&new_leaf_node.bloom_filter_path)
                .unwrap();
            let mut bf_new_leaf_node = b4_new_leaf_node.write().unwrap();

            for kmer in &genome.kmers {
                bf_new_leaf_node.insert(&kmer);
            }
        }

        log::debug!("(bloom tree; init()) Adding {}", genome.id);
        return new_leaf_node;
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
            let right_distance = self.node_distance(&right, &node);
            let left_distance = self.node_distance(&left, &node);

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
        let node_name = "Internal Node".to_string() + "_" + &n2.to_string();
        let mut new_internal_node = self.make_bloom_node(node_name.clone());

        // union children nodes into new parent internal node
        self.node_union(&new_internal_node, &new_node);
        self.node_union(&new_internal_node, &current_node);

        // assign children nodes to the new internaal node.
        new_internal_node.left_child = Some(current_node);
        new_internal_node.right_child = Some(new_node);

        return new_internal_node;
    }
    
    fn node_distance(&self, node_a: &BloomNode, node_b: &BloomNode) -> usize {
            // get BF for node a
            let b4_node_a = self.get_bf(node_a);
            let bf_node_a = b4_node_a.read().unwrap();
    
            // get BF for node b
            let b4_node_b = self.get_bf(node_b);
            let bf_node_b = b4_node_b.read().unwrap();

            return bf_node_a.distance(&bf_node_b);
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
        self.bf_cache.get_filter(&bloom_node.bloom_filter_path).unwrap()
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

        return deserialized_tree;
    }

    fn make_bloom_node(&self, nodeid: String) -> Box<BloomNode> {
        // bloom filter path.
        let bloom_filter_path = PathBuf::from(nodeid.clone() + ".bf");
        let full_bloom_filter_path = self.directory.clone().unwrap().join(&bloom_filter_path);

        // create and save the bloom filter.
        let bloom_filter = create_bloom_filter(
            self.hash_states.clone(),
            full_bloom_filter_path.clone(),
            self.false_pos_rate,
            self.largest_expected_genome,
        );

        // Cache the loaded Bloom filter
        self.bf_cache
            .add_filter(&full_bloom_filter_path, bloom_filter);

        // make node.
        let bloomnode = BloomNode::new(Some(nodeid.clone()), full_bloom_filter_path);
        Box::new(bloomnode)
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
            bloom_filter_path: bloom_filter_path.clone(),
            mapped_reads: 0,
        }
    }

    /// Returns a boolean indicating whether the node is a leaf node.
    ///
    /// # Returns
    ///
    /// `true` if the node is a leaf node; `false` otherwise.
    pub(crate) fn is_leafnode(&self) -> bool {
        return self.left_child.is_none() && self.right_child.is_none();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bloom_filter::BloomFilter;
    use crate::cache;
    use crate::file_parser::RecordTypes;
    use bio::io::fasta;
    use std::fs;

    fn get_tmp_dir() -> std::path::PathBuf {
        PathBuf::from("tmp_testing/")
    }

    struct TestContext {}

    impl TestContext {
        fn setup() -> TestContext {
            println!("Setting up the test context");
            let tmp_direcory = get_tmp_dir();
            if tmp_direcory.is_dir() {
                fs::remove_dir_all(&tmp_direcory);
            }

            fs::create_dir(&tmp_direcory).unwrap();
            TestContext {}
        }

        fn get_tmp_dir(&self) -> std::path::PathBuf {
            PathBuf::from("tmp_testing/")
        }

        fn teardown(&self) -> () {
            println!("Tearing down the test context");
            fs::remove_dir_all(get_tmp_dir()).unwrap();
        }
    }

    fn create_bloom_filter(false_positive_rate: f32, largest_expected_genome: u32) -> BloomFilter {
        let hash_states = (HashSeed::new(), HashSeed::new());
        BloomFilter::with_rate(false_positive_rate, largest_expected_genome, hash_states)
    }

    #[test]
    fn test_empty_tree_insert() {
        // setup context
        let context = TestContext::setup();

        // initialization states.
        let kmer_size = 5;
        let directory = context.get_tmp_dir();
        let cache_size = 5;
        let false_positive_rate = 0.001;
        let largest_expected_genome = 1000;
        assert!(directory.exists());

        // create a DNA sequence.
        let record_id = "tmp_bloom";
        let kmers = file_parser::get_kmers(&"ATCAG".as_bytes().to_vec(), &kmer_size);
        let record = file_parser::DNASequence::new(
            Some("ATCAG".as_bytes().to_vec()),
            record_id.to_string(),
            kmers,
        );

        // create a new cache
        let bloomfilter_cache = cache::BFLruCache::new(cache_size);

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
            directory: Some(directory),
            kmer_size: kmer_size,
            hash_states: tree.hash_states.clone(),
            bf_cache: Box::new(cache::BFLruCache::new(cache_size)),
            false_pos_rate: false_positive_rate,
            largest_expected_genome,
        };

        assert_eq!(tree, expected_tree);

        //
        // let tree = expected_tree;
        // context.teardown();
    }
}

//     #[test]
//     fn test_one_elem_tree_insert() {
//         let kmer_size = 5;
//         let mut tree = BloomTree::new(kmer_size);

//         let records = [
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test1", None, "ATCAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test2", None, "TTTAG".as_ref())),
//         ];

//         for record in &records {
//             tree = tree.insert(record);
//         }

//         let expected_leaves = records.map(|r| tree.init_leaf(&r));
//         let mut expected_root = Box::new(BloomNode::new(
//             tree.hash_states.clone(),
//             tree.root.clone().unwrap().tax_id,
//         ));
//         // The leaves should under a parent node that has the union of their bloom filters
//         for leaf in &expected_leaves {
//             expected_root.bloom_filter.union(&leaf.bloom_filter);
//         }
//         [expected_root.left_child, expected_root.right_child] =
//             expected_leaves.map(|node| Some(node));

//         let expected_tree = BloomTree {
//             root: Some(expected_root),
//             kmer_size: kmer_size,
//             hash_states: tree.hash_states.clone(),
//         };

//         assert_eq!(tree, expected_tree);
//     }

//     #[test]
//     fn test_nested_tree_insert_left() {
//         let kmer_size = 5;
//         let mut tree = BloomTree::new(kmer_size);

//         let records = [
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test1", None, "ATCAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test2", None, "TTTAG".as_ref())),
//             // Has the same sequence as the first so the bloom filter should be closer to that leaf
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test3", None, "ATCAG".as_ref())),
//         ];

//         for record in &records {
//             tree = tree.insert(record);
//         }

//         let expected_leaves = records.map(|r| tree.init_leaf(&r));
//         // Expected tree should now be:
//         //                        (_, combined_test1_test2_test3)
//         //                              /                  \
//         //                            /                     \
//         //            (_, combined_test1_test3)       (test2, TTTAG)
//         //              /              \
//         //            /                 \
//         //    (test1, ATCAG)      (test3, ATCAG)
//         let mut expected_root = Box::new(BloomNode::new(
//             tree.hash_states.clone(),
//             tree.root.clone().unwrap().tax_id,
//         ));
//         // The leaves should under a parent node that has the union of their bloom filters
//         for leaf in &expected_leaves {
//             expected_root.bloom_filter.union(&leaf.bloom_filter);
//         }
//         let mut inner_left_node = Box::new(BloomNode::new(
//             tree.hash_states.clone(),
//             tree.root
//                 .clone()
//                 .unwrap()
//                 .left_child
//                 .clone()
//                 .unwrap()
//                 .tax_id,
//         ));
//         // The left inner node's bloom filter should be identical to the first and third leaves bloom filters (since they have the same sequence)
//         inner_left_node
//             .bloom_filter
//             .union(&expected_leaves[0].bloom_filter);
//         [
//             inner_left_node.left_child,
//             expected_root.right_child,
//             inner_left_node.right_child,
//         ] = expected_leaves.map(|node| Some(node));
//         expected_root.left_child = Some(inner_left_node);

//         let expected_tree = BloomTree {
//             root: Some(expected_root),
//             kmer_size: kmer_size,
//             hash_states: tree.hash_states.clone(),
//         };

//         assert_eq!(tree, expected_tree);
//     }

//     #[test]
//     fn test_nested_tree_insert_right() {
//         let kmer_size = 5;
//         let mut tree = BloomTree::new(kmer_size);

//         let records = [
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test1", None, "ATCAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test2", None, "TTTAG".as_ref())),
//             // Has the same sequence as the second so the bloom filter should be closer to that leaf
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test3", None, "TTTAG".as_ref())),
//         ];

//         for record in &records {
//             tree = tree.insert(record);
//         }

//         let expected_leaves = records.map(|r| tree.init_leaf(&r));
//         // Expected tree should now be:
//         //    (_, combined_test1_test2_test3)
//         //          /                  \
//         //        /                     \
//         // (test1, ATCAG)    (_, combined_test2_test3)
//         //                       /              \
//         //                     /                 \
//         //              (test2, TTTAG)     (test3, TTTAG)
//         let mut expected_root = Box::new(BloomNode::new(
//             tree.hash_states.clone(),
//             tree.root.clone().unwrap().tax_id,
//         ));
//         // The leaves should under a parent node that has the union of their bloom filters
//         for leaf in &expected_leaves {
//             expected_root.bloom_filter.union(&leaf.bloom_filter);
//         }
//         let mut inner_right_node = Box::new(BloomNode::new(
//             tree.hash_states.clone(),
//             tree.root
//                 .clone()
//                 .unwrap()
//                 .right_child
//                 .clone()
//                 .unwrap()
//                 .tax_id,
//         ));
//         // The right inner node's bloom filter should be identical to the second and third leaves bloom filters (since they have the same sequence)
//         inner_right_node
//             .bloom_filter
//             .union(&expected_leaves[1].bloom_filter);
//         [
//             expected_root.left_child,
//             inner_right_node.left_child,
//             inner_right_node.right_child,
//         ] = expected_leaves.map(|node| Some(node));
//         expected_root.right_child = Some(inner_right_node);

//         let expected_tree = BloomTree {
//             root: Some(expected_root),
//             kmer_size: kmer_size,
//             hash_states: tree.hash_states.clone(),
//         };

//         assert_eq!(
//             tree, expected_tree,
//             "Test failed for {:#?} and {:#?}",
//             tree, expected_tree
//         );
//     }

//     #[test]
//     fn test_save_load() {
//         let kmer_size = 4;
//         let rand_suffix: usize = rand::random();
//         // We add a random suffix to prevent parallel tests from conflicting with each other.
//         let scratch_dir_name = format!("test_tmp_{}/", rand_suffix);
//         let scratch_dir = Path::new(&scratch_dir_name);
//         let records1 = vec![
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test1", None, "ATCAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test2", None, "TTTAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test3", None, "TTTAG".as_ref())),
//         ];

//         let records2 = vec![
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test1", None, "GATCAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test2", None, "GTTTAG".as_ref())),
//             RecordTypes::FastaRecord(fasta::Record::with_attrs("test3", None, "GTTTAG".as_ref())),
//         ];

//         let tree1 = create_bloom_tree(records1, &kmer_size);

//         // Should be able to save and load the same tree
//         tree1.save(scratch_dir);
//         assert_eq!(BloomTree::load(scratch_dir), tree1);

//         // Saving overwrites existing saved tree
//         let tree2 = create_bloom_tree(records2, &kmer_size);
//         tree2.save(scratch_dir);
//         assert_eq!(BloomTree::load(scratch_dir), tree2);

//         // Cleanup scratch dir
//         std::fs::remove_dir_all(scratch_dir).unwrap();
//     }
// }
