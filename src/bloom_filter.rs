/// PhageFilter - bloom filter
///
/// Note:
///     These methods have been modified from the original
///     rust-bloom package.
///
///
use bit_vec::BitVec;
mod hash_iter;
use std::cmp::{max, min};
use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hash};

pub fn get_bloom_filter(genome_count: usize) -> BloomFilter {
    let max_genome_size: u32 = 1000000; // max size of phage genome. TODO: find this.
    let expected_num_items: u32 = (genome_count as u32) * max_genome_size;
    println!("number of items expected: {}\n", expected_num_items);
    // out of 100 items that are not inserted, expect 0.1 to return true for contain
    let false_positive_rate: f32 = 0.001;
    // instantiate a BloomFilter
    let filter = BloomFilter::with_rate(false_positive_rate, expected_num_items);
    return filter;
}

pub trait ASMS {
    fn insert<T: Hash>(&mut self, item: &T) -> bool;
    fn contains<T: Hash>(&self, item: &T) -> bool;
    fn clear(&mut self);
}

pub struct BloomFilter<R = RandomState, S = RandomState> {
    bits: BitVec,
    num_hashes: u32,
    hash_builder_one: R,
    hash_builder_two: S,
}
pub trait DistanceChecker {
    fn distance(&mut self, other: &Self) -> u8;
}

impl DistanceChecker for BloomFilter<RandomState, RandomState> {
    /// Calculates the distance between two
    /// bloom filters using the hamming_distance.
    fn distance(&mut self, other: &BloomFilter) -> u8 {
        let mut diff_1: BitVec = self.bits.clone();
        let mut diff_2: BitVec = other.bits.clone();
        diff_1.difference(&diff_2);
        diff_2.difference(&diff_1);
        diff_1.or(&diff_2);
        let hamming_distance: u8 = diff_1.iter().filter(|x| *x).count() as u8;
        return hamming_distance;
    }
}

impl BloomFilter<RandomState, RandomState> {
    /// Create a new BloomFilter with the specified number of bits,
    /// and hashes
    pub fn with_size(num_bits: usize, num_hashes: u32) -> BloomFilter<RandomState, RandomState> {
        BloomFilter {
            bits: BitVec::from_elem(num_bits, false),
            num_hashes: num_hashes,
            hash_builder_one: RandomState::new(),
            hash_builder_two: RandomState::new(),
        }
    }

    /// create a BloomFilter that expects to hold
    /// `expected_num_items`.  The filter will be sized to have a
    /// false positive rate of the value specified in `rate`.
    pub fn with_rate(rate: f32, expected_num_items: u32) -> BloomFilter<RandomState, RandomState> {
        let bits = needed_bits(rate, expected_num_items);
        BloomFilter::with_size(bits, optimal_num_hashes(bits, expected_num_items))
    }

    /// Get the number of bits this BloomFilter is using
    pub fn num_bits(&self) -> usize {
        self.bits.len()
    }

    /// Get the number of hash functions this BloomFilter is using
    pub fn num_hashes(&self) -> u32 {
        self.num_hashes
    }

    /// Calculates the intersection of two BloomFilters.  Only items inserted into both filters will still be present in `self`.
    ///
    /// Both BloomFilters must be using the same number of
    /// bits. Returns true if self changed.
    ///
    /// # Panics
    /// Panics if the BloomFilters are not using the same number of bits
    fn intersect(&mut self, other: &BloomFilter) -> bool {
        self.bits.and(&other.bits)
    }

    /// Calculates the union of two BloomFilters.  Items inserted into
    /// either filters will be present in `self`.
    ///
    /// Both BloomFilters must be using the same number of
    /// bits. Returns true if self changed.
    ///
    /// # Panics
    /// Panics if the BloomFilters are not using the same number of bits
    fn union(&mut self, other: &BloomFilter) -> bool {
        self.bits.or(&other.bits)
    }
}

impl<R, S> ASMS for BloomFilter<R, S>
where
    R: BuildHasher,
    S: BuildHasher,
{
    /// Insert item into this BloomFilter.
    ///
    /// If the BloomFilter did not have this value present, `true` is returned.
    ///
    /// If the BloomFilter did have this value present, `false` is returned.
    fn insert<T: Hash>(&mut self, item: &T) -> bool {
        let mut contained = true;
        for h in hash_iter::HashIter::from(
            item,
            self.num_hashes,
            &self.hash_builder_one,
            &self.hash_builder_two,
        ) {
            let idx = (h % self.bits.len() as u64) as usize;
            match self.bits.get(idx) {
                Some(b) => {
                    if !b {
                        contained = false;
                    } else {
                        contained = true;
                    }
                }
                None => {
                    panic!("Hash mod failed in insert");
                }
            }
            self.bits.set(idx, true)
        }
        !contained
    }

    /// Check if the item has been inserted into this bloom filter.
    /// This function can return false positives, but not false
    /// negatives.
    fn contains<T: Hash>(&self, item: &T) -> bool {
        for h in hash_iter::HashIter::from(
            item,
            self.num_hashes,
            &self.hash_builder_one,
            &self.hash_builder_two,
        ) {
            let idx = (h % self.bits.len() as u64) as usize;
            match self.bits.get(idx) {
                Some(b) => {
                    if !b {
                        return false;
                    }
                }
                None => {
                    panic!("Hash mod failed - out of bounds for given bit vec");
                }
            }
        }
        true
    }

    /// Remove all values from this BloomFilter
    fn clear(&mut self) {
        self.bits.clear();
    }
}

/// Return the optimal number of hashes to use for the given number of
/// bits and items in a filter
pub fn optimal_num_hashes(num_bits: usize, num_items: u32) -> u32 {
    min(
        max(
            (num_bits as f32 / num_items as f32 * core::f32::consts::LN_2).round() as u32,
            2,
        ),
        200,
    )
}

/// Return the number of bits needed to satisfy the specified false
/// positive rate, if the filter will hold `num_items` items.
pub fn needed_bits(false_pos_rate: f32, num_items: u32) -> usize {
    let ln22 = core::f32::consts::LN_2 * core::f32::consts::LN_2;
    (num_items as f32 * ((1.0 / false_pos_rate).ln() / ln22)).round() as usize
}
