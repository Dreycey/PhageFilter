use crate::bloom_filter::BloomFilter;
use crate::Path;
use crate::PathBuf;
///
/// This module contains the bloom filter cache.
///
use lru::LruCache;
use std::num::NonZeroUsize;
use std::sync::{Arc, RwLock};

#[derive(Debug)]
pub struct BloomFilterCache {
    cache: RwLock<LruCache<PathBuf, Arc<RwLock<BloomFilter>>>>,
    capacity: usize,
}

impl BloomFilterCache {
    pub fn new(capacity: usize) -> Self {
        Self {
            cache: RwLock::new(LruCache::new(NonZeroUsize::new(capacity).unwrap())),
            capacity,
        }
    }

    pub fn get_filter(&self, key: &PathBuf) -> Option<Arc<RwLock<BloomFilter>>> {
        println!("trying to get filter: {:?}", key.to_path_buf());
        let mut cache_lock = self.cache.write().unwrap();
        if let Some(filter) = cache_lock.get(key) {
            Some(filter.clone())
        } else {
            let filter_path = Path::new(key);
            if filter_path.exists() {
                println!("the following is a saved file: {:?}", filter_path);
                let filter: BloomFilter = BloomFilter::load_from_file(&filter_path);
                println!("file loded");
                let arc_filter = Arc::new(RwLock::new(filter));
                cache_lock.put(key.to_path_buf(), arc_filter.clone());
                Some(arc_filter)
            } else {
                println!("the following is a no-go: {:?}", key.to_path_buf());
                None
            }
        }
    }

    pub fn add_filter(&self, key: &PathBuf, bloom_filter: BloomFilter) {
        let mut cache_lock = self.cache.write().unwrap();
        let arc_filter = Arc::new(RwLock::new(bloom_filter));
        cache_lock.put(key.to_path_buf(), arc_filter);
    }
}
