use rand;
use rustc_hash::FxHasher;
use serde::{Deserialize, Serialize};
use std::hash::{BuildHasher, Hasher};

// We make our own hash seed since std::collections::hash_map::RandomState can't be serialized and deserialized.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct HashSeed {
    seed: usize,
}

impl BuildHasher for HashSeed {
    type Hasher = FxHasher;

    fn build_hasher(&self) -> Self::Hasher {
        let mut hasher = FxHasher::default();
        // Note: the FxHasher has no seed (it always starts from zero) so we instead hash the seed to change the hashers internal state.
        hasher.write_usize(self.seed);
        hasher
    }
}

impl HashSeed {
    pub fn new() -> Self {
        HashSeed {
            seed: rand::random(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hasher_different_seeds_return_different_hashes() {
        let seed1 = HashSeed { seed: 5 };
        let seed2 = HashSeed { seed: 10 };

        let test_string = "Hello world!";

        let mut hasher1 = seed1.build_hasher();
        let mut hasher2 = seed2.build_hasher();
        hasher1.write(test_string.as_bytes());
        hasher2.write(test_string.as_bytes());

        assert_ne!(hasher1.finish(), hasher2.finish(),)
    }
}
