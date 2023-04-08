// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.

/// A standard BloomFilter.  If an item is instered then `contains`
/// is guaranteed to return `true` for that item.  For items not
/// inserted `contains` will probably return false.  The probability
/// that `contains` returns `true` for an item that was not inserted
/// is called the False Positive Rate.
///
/// # False Positive Rate
/// The false positive rate is specified as a float in the range
/// (0,1).  If indicates that out of `X` probes, `X * rate` should
/// return a false positive.  Higher values will lead to smaller (but
/// more inaccurate) filters.
///
/// # Example Usage
///
/// ```rust
/// use bloom::{ASMS,BloomFilter};
///
/// let expected_num_items = 1000;
///
/// // out of 100 items that are not inserted, expect 1 to return true for contain
/// let false_positive_rate = 0.01;
///
/// let mut filter = BloomFilter::with_rate(false_positive_rate,expected_num_items);
/// filter.insert(&1);
/// filter.contains(&1); /* true */
/// filter.contains(&2); /* false */
/// ```
use bitvec::prelude::*;
mod hash_iter;
pub mod hasher;
use hasher::HashSeed;
use serde::{Deserialize, Serialize};
use std::any::TypeId;
use std::cmp::{max, min};
use std::fmt::{Debug, Formatter};
use std::fs::File;
use std::hash::{BuildHasher, Hash};
use std::io::Write;
use std::io::{BufReader, BufWriter};
use std::ops::Drop;
use std::path::{Path, PathBuf};

pub fn create_bloom_filter(
    hash_states: (HashSeed, HashSeed),
    full_bf_path: PathBuf,
    false_pos_rate: f32,
) -> BloomFilter {
    let expected_num_items: u32 = 1_000_000;
    log::debug!(
        "(Building BF) number of items expected: {}\n",
        expected_num_items
    );
    // instantiate a BloomFilter
    let mut filter = BloomFilter::with_rate(false_pos_rate, expected_num_items, hash_states);
    filter.file_path = Some(full_bf_path);

    return filter;
}

/// Approximate Set Membership Structure
pub trait ASMS {
    fn insert<T: Hash>(&mut self, item: &T) -> bool;
    fn contains<T: Hash>(&self, item: &T) -> bool;
    fn clear(&mut self);
}

pub trait DistanceChecker {
    fn distance(&self, other: &Self) -> usize;
}

#[derive(Clone, Deserialize, Serialize)]
pub struct BloomFilter<R = HashSeed, S = HashSeed> {
    bits: BitVec,
    num_hashes: u32,
    hash_builder_one: R,
    hash_builder_two: S,
    file_path: Option<PathBuf>,
}

impl<R, S> BloomFilter<R, S> {
    fn as_bloom_filter(&mut self) -> &mut BloomFilter {
        unsafe { &mut *(self as *mut Self as *mut BloomFilter) }
    }
}

impl<R, S> Drop for BloomFilter<R, S> {
    fn drop(&mut self) {
        let bloomfilter_inner = self.as_bloom_filter();
        let path = bloomfilter_inner.file_path.as_ref().unwrap().as_path();
        bloomfilter_inner.save_to_file(path);
    }
}

/// Equality for bloom filters is judged only using the bits field
impl PartialEq for BloomFilter<HashSeed, HashSeed> {
    fn eq(&self, other: &Self) -> bool {
        self.bits == other.bits
    }
}

impl Debug for BloomFilter<HashSeed, HashSeed> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BloomFilter")
            // TODO: Omitted bits representation for brevity because the number of bits in the
            //  bitvector is currently >1,000,000 which takes forever to print out and takes up
            //  unnecessary space
            .field("bits", &"[omitted for brevity...]")
            .field("num_hashes", &self.num_hashes)
            .field("hash_builder_one", &self.hash_builder_one)
            .field("hash_builder_two", &self.hash_builder_two)
            .finish()
    }
}

impl DistanceChecker for BloomFilter<HashSeed, HashSeed> {
    /// Calculates the distance between two bloom filters using the hamming_distance.
    fn distance(&self, other: &BloomFilter) -> usize {
        let diff = self.bits.clone() ^ &other.bits;
        diff.count_ones()
    }
}

impl BloomFilter<HashSeed, HashSeed> {
    pub fn load_from_file(bloom_filter_path: &Path) -> Self {
        // Open the file at the given path
        let file = File::open(bloom_filter_path).unwrap_or_else(|_| {
            panic!("Failed to open Bloom filter file: {:?}", bloom_filter_path)
        });

        // Create a buffered reader for the file
        let reader = BufReader::new(file);

        // Deserialize the Bloom filter from the reader
        let bloom_filter: BloomFilter = bincode::deserialize_from(reader).unwrap_or_else(|_| {
            panic!(
                "Failed to deserialize Bloom filter from file: {:?}",
                bloom_filter_path
            )
        });

        bloom_filter
    }

    pub fn save_to_file(&self, bloom_filter_path: &Path) {
        // Create the file at the given path
        let file = File::create(bloom_filter_path).unwrap_or_else(|_| {
            panic!(
                "Failed to create Bloom filter file: {:?}",
                bloom_filter_path
            )
        });

        // Create a buffered writer for the file
        let mut writer = BufWriter::new(file);

        // Serialize the Bloom filter into the writer
        bincode::serialize_into(&mut writer, &self).unwrap_or_else(|_| {
            panic!(
                "Failed to serialize Bloom filter to file: {:?}",
                bloom_filter_path
            )
        });

        // Flush the writer to ensure the data is written to the file
        writer.flush().unwrap_or_else(|_| {
            panic!("Failed to flush Bloom filter file: {:?}", bloom_filter_path)
        });
    }

    /// Create a new BloomFilter with the specified number of bits,
    /// and hashes
    pub fn with_size(
        num_bits: usize,
        num_hashes: u32,
        hash_states: (HashSeed, HashSeed),
    ) -> BloomFilter<HashSeed, HashSeed> {
        let (hash_builder_one, hash_builder_two) = hash_states;
        BloomFilter {
            // Initialize all bits to zero
            bits: bitvec![0; num_bits],
            num_hashes: num_hashes,
            hash_builder_one,
            hash_builder_two,
            file_path: None, // TODO: update
        }
    }

    /// create a BloomFilter that expects to hold
    /// `expected_num_items`.  The filter will be sized to have a
    /// false positive rate of the value specified in `rate`.
    pub fn with_rate(
        rate: f32,
        expected_num_items: u32,
        hash_states: (HashSeed, HashSeed),
    ) -> BloomFilter<HashSeed, HashSeed> {
        let bits = needed_bits(rate, expected_num_items);
        BloomFilter::with_size(
            bits,
            optimal_num_hashes(bits, expected_num_items),
            hash_states,
        )
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
    /// bits.
    ///
    /// # Panics
    /// Panics if the BloomFilters are not using the same number of bits
    fn intersect(&mut self, other: &BloomFilter) {
        self.bits &= &other.bits;
    }

    /// Calculates the union of two BloomFilters.  Items inserted into
    /// either filters will be present in `self`.
    ///
    /// Both BloomFilters must be using the same number of
    /// bits.
    ///
    /// # Panics
    /// Panics if the BloomFilters are not using the same number of bits
    pub(crate) fn union(&mut self, other: &BloomFilter) {
        self.bits |= &other.bits;
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
                    panic!("Hash mod failed");
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

#[cfg(test)]
mod tests {
    use super::*;

    // Build a static bloom filter for testing. Not meant to be inserted into.
    fn get_static_bloom_filter(bits: BitVec) -> BloomFilter {
        let state = HashSeed::new();

        BloomFilter {
            bits,
            num_hashes: 0,
            hash_builder_one: state.clone(),
            hash_builder_two: state.clone(),
            file_path: None,
        }
    }

    #[test]
    fn test_distance() {
        let b1 = get_static_bloom_filter(BitVec::from_bitslice(0b00101101.view_bits()));
        let b2 = get_static_bloom_filter(BitVec::from_bitslice(0b10100111.view_bits()));
        let expected_distance = 3;

        let b_none = get_static_bloom_filter(BitVec::from_bitslice(0b00000000.view_bits()));
        let b_all = get_static_bloom_filter(BitVec::from_bitslice(0b11111111.view_bits()));

        assert_eq!(b1.distance(&b2), expected_distance);
        assert_eq!(b2.distance(&b1), expected_distance);
        assert_eq!(b1.distance(&b1), 0);
        assert_eq!(b2.distance(&b2), 0);
        assert_eq!(b_none.distance(&b_all), 8);
    }
}
