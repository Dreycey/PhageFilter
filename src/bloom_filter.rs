pub use bloom::{BloomFilter, ASMS};

pub fn get_bloom_filter() -> BloomFilter {
    let expected_num_items: u32 = 1000000;

    // // out of 100 items that are not inserted, expect 1 to return true for contain
    let false_positive_rate: f32 = 0.01;

    let filter = BloomFilter::with_rate(false_positive_rate, expected_num_items);

    return filter;
}
