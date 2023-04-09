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

#[derive(Debug)]
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
            .field("Cache capactity", &self.get_capacity())
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
