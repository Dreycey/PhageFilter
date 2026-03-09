// utilities for hashing
use std::hash::{BuildHasher, Hash};
pub(in crate::bloom_filter) struct HashIter {
    h1: u64,
    h2: u64,
    i: u32,
    count: u32,
}

impl Iterator for HashIter {
    type Item = u64;

    fn next(&mut self) -> Option<u64> {
        if self.i == self.count {
            return None;
        }
        let r = match self.i {
            0 => self.h1,
            1 => self.h2,
            _ => {
                let p1 = self.h1.wrapping_add(self.i as u64);
                p1.wrapping_mul(self.h2)
            }
        };
        self.i += 1;
        Some(r)
    }
}

impl HashIter {
    pub fn from<T: Hash, R: BuildHasher, S: BuildHasher>(
        item: T,
        count: u32,
        build_hasher_one: &R,
        build_hasher_two: &S,
    ) -> HashIter {
        let h1 = build_hasher_one.hash_one(&item);
        let h2 = build_hasher_two.hash_one(&item);
        HashIter {
            h1,
            h2,
            i: 0,
            count,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bloom_filter::hasher::HashSeed;

    fn make_iter(item: &str, count: u32) -> (HashIter, u64, u64) {
        let s1 = HashSeed::new();
        let s2 = HashSeed::new();
        let probe = HashIter::from(item, count, &s1, &s2);
        let h1 = {
            s1.hash_one(item)
        };
        let h2 = {
            s2.hash_one(item)
        };
        (probe, h1, h2)
    }

    #[test]
    fn test_count_items() {
        for count in [0, 1, 2, 5] {
            let (iter, _, _) = make_iter("hello", count);
            assert_eq!(iter.collect::<Vec<_>>().len(), count as usize);
        }
    }

    #[test]
    fn test_first_is_h1_second_is_h2() {
        let (iter, h1, h2) = make_iter("hello", 2);
        let vals: Vec<u64> = iter.collect();
        assert_eq!(vals[0], h1);
        assert_eq!(vals[1], h2);
    }

    #[test]
    fn test_formula_for_i_ge_2() {
        let (iter, h1, h2) = make_iter("world", 5);
        let vals: Vec<u64> = iter.collect();
        for i in 2u32..5 {
            let expected = h1.wrapping_add(i as u64).wrapping_mul(h2);
            assert_eq!(vals[i as usize], expected, "mismatch at i={}", i);
        }
    }

    #[test]
    fn test_different_seeds_produce_different_sequences() {
        let s1a = HashSeed::new();
        let s2a = HashSeed::new();
        let s1b = HashSeed::new();
        let s2b = HashSeed::new();
        let a: Vec<u64> = HashIter::from("test", 5, &s1a, &s2a).collect();
        let b: Vec<u64> = HashIter::from("test", 5, &s1b, &s2b).collect();
        assert_ne!(a, b, "different seeds should produce different sequences");
    }
}
