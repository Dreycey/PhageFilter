use crate::bloom_filter::BloomFilter;
use crate::Path;
use crate::PathBuf;
///
/// This module contains the bloom filter cache.
///
use lru::LruCache;
use std::fmt;
use std::num::NonZeroUsize;
use std::sync::{Arc, RwLock};

//#[derive(Debug)]
pub struct BFLruCache {
    cache: RwLock<LruCache<PathBuf, Arc<RwLock<BloomFilter>>>>,
    capacity: usize,
}

/// Bloom filter cache trait defines the required
/// methods for a cache. This allows for ensuring the
/// BllomTree is loosely coupled with the cache. Additionally,
/// it allows for dependency injection and creating Mocks for testing.
pub trait BloomFilterCache {
    fn get_filter(&self, key: &PathBuf) -> Option<Arc<RwLock<BloomFilter>>>;
    fn get_capacity(&self) -> usize;
    fn add_filter(&self, key: &PathBuf, bloom_filter: BloomFilter);
}

impl fmt::Debug for dyn BloomFilterCache {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BloomFilterCache")
            .field("capacity", &self.get_capacity())
            .finish()
    }
}

impl BFLruCache {
    pub fn new(capacity: usize) -> Self {
        Self {
            cache: RwLock::new(LruCache::new(NonZeroUsize::new(capacity).unwrap())),
            capacity,
        }
    }
}

impl fmt::Debug for BFLruCache {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("BFLruCache")
            .field("capacity", &self.get_capacity())
            .finish()
    }
}

impl BloomFilterCache for BFLruCache {
    fn get_filter(&self, key: &PathBuf) -> Option<Arc<RwLock<BloomFilter>>> {
        let mut cache_lock = self.cache.write().unwrap();
        if let Some(filter) = cache_lock.get(key) {
            Some(filter.clone())
        } else {
            let filter_path = Path::new(key);
            if filter_path.exists() {
                let filter: BloomFilter = BloomFilter::load_from_file(&filter_path);
                let arc_filter = Arc::new(RwLock::new(filter));
                cache_lock.put(key.to_path_buf(), arc_filter.clone());
                Some(arc_filter)
            } else {
                println!("the following is a no-go: {:?}", key.to_path_buf());
                None
            }
        }
    }

    fn get_capacity(&self) -> usize {
        self.capacity
    }

    fn add_filter(&self, key: &PathBuf, bloom_filter: BloomFilter) {
        let mut cache_lock = self.cache.write().unwrap();
        let arc_filter = Arc::new(RwLock::new(bloom_filter));
        cache_lock.put(key.to_path_buf(), arc_filter);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bloom_filter;
    use crate::bloom_filter::hasher::HashSeed;
    use proptest::prelude::*;
    use proptest::sample::select;
    use std::fs;

    fn create_bloom_filter() -> BloomFilter {
        let hash_states = (HashSeed::new(), HashSeed::new());
        let false_pos_rate = 0.1;
        let largest_expected_genome = 100;
        BloomFilter::with_rate(false_pos_rate, largest_expected_genome, hash_states)
    }

    #[test]
    fn test_get_filter() {
        let key: PathBuf = PathBuf::from("TEST");
        let cache = BFLruCache::new(2);
        let bloomfilter = create_bloom_filter();
        cache.add_filter(&key, bloomfilter);
        assert!(cache.get_filter(&key).is_some());
    }

    #[test]
    fn test_bad_filter() {
        let key: PathBuf = PathBuf::from("TEST");
        let cache = BFLruCache::new(1);
        assert!(cache.get_filter(&key).is_none());
    }

    #[test]
    fn test_get_filter_from_disk() {
        // create a bloom filter
        let bf_disk_key: PathBuf = PathBuf::from("tmp_bloomfilter.bf");
        let hash_states = (HashSeed::new(), HashSeed::new());
        let false_pos_rate = 0.1;
        let largest_expected_genome = 100;
        let bloomfilter = bloom_filter::create_bloom_filter(
            hash_states,
            bf_disk_key.clone(),
            false_pos_rate,
            largest_expected_genome,
        );

        // create cache
        let cache = BFLruCache::new(1);
        cache.add_filter(&bf_disk_key, bloomfilter);

        // kick out first bloom filter
        let key: PathBuf = PathBuf::from("TEST");
        cache.add_filter(&key, create_bloom_filter());

        // test pulling from disk works
        assert!(cache.get_filter(&bf_disk_key).is_some());

        // delete temp file
        cache.add_filter(&key, create_bloom_filter());
        assert!(fs::remove_file(bf_disk_key).is_ok())
    }

    proptest! {
        #[test]
        fn test_get_capacity(capacity in select((0..=500_000).collect::<Vec<usize>>())) {
            let cache = BFLruCache::new(capacity);
            prop_assert_eq!(cache.get_capacity(), capacity);
        }

        #[test]
        fn test_cache_debug(capacity in select((10..=500_000).collect::<Vec<usize>>())) {
            let cache = BFLruCache::new(capacity);
            let debug = format!("{:?}", cache);
            let expected = format!("BFLruCache {{ capacity: {} }}", cache.get_capacity());
            prop_assert_eq!(debug, expected);
        }
    }
}
