///
/// This struct is used to connect reads to
/// genomes and map the output to a designated file.
///
///
use std::collections::HashMap;
use std::collections::HashSet;

pub struct ResultMap {
    read_map: HashMap<String, HashSet<String>>,
}

impl ResultMap {
    pub fn new() -> ResultMap {
        ResultMap {
            read_map: HashMap::new(),
        }
    }

    pub fn add_read_map(&mut self, read_id: String, genome_id: String) {
        self.read_map.entry(read_id).or_default().insert(genome_id);
    }

    pub fn get_ext_id(&self, read_id: &str) -> String {
        let genomes = self
            .read_map
            .get(read_id)
            .map(|genome_set| {
                genome_set
                    .iter()
                    .map(AsRef::as_ref)
                    .collect::<Vec<&str>>()
                    .join(",")
            })
            .unwrap_or_default();
        format!("{} |{}", read_id, genomes)
    }

    pub fn read_mapped(&self, read_id: &str) -> bool {
        self.read_map.contains_key(read_id)
    }

    pub fn empty_read_map(&mut self) {
        self.read_map.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_creates_empty_map() {
        let rm = ResultMap::new();
        assert!(!rm.read_mapped("any"));
    }

    #[test]
    fn test_add_read_map_single_insert() {
        let mut rm = ResultMap::new();
        rm.add_read_map("read1".to_string(), "genomeA".to_string());
        assert!(rm.read_mapped("read1"));
        assert_eq!(rm.read_map.get("read1").unwrap().len(), 1);
    }

    #[test]
    fn test_add_read_map_duplicate_read_different_genome() {
        let mut rm = ResultMap::new();
        rm.add_read_map("read1".to_string(), "genomeA".to_string());
        rm.add_read_map("read1".to_string(), "genomeB".to_string());
        let genomes = rm.read_map.get("read1").unwrap();
        assert_eq!(genomes.len(), 2);
        assert!(genomes.contains("genomeA"));
        assert!(genomes.contains("genomeB"));
    }

    #[test]
    fn test_get_ext_id_unmapped_read() {
        let rm = ResultMap::new();
        let ext = rm.get_ext_id("unknown");
        assert_eq!(ext, "unknown |");
    }

    #[test]
    fn test_get_ext_id_single_genome() {
        let mut rm = ResultMap::new();
        rm.add_read_map("read1".to_string(), "genomeA".to_string());
        let ext = rm.get_ext_id("read1");
        assert_eq!(ext, "read1 |genomeA");
    }

    #[test]
    fn test_get_ext_id_multiple_genomes() {
        let mut rm = ResultMap::new();
        rm.add_read_map("read1".to_string(), "genomeA".to_string());
        rm.add_read_map("read1".to_string(), "genomeB".to_string());
        let ext = rm.get_ext_id("read1");
        assert!(ext.starts_with("read1 |"));
        let suffix = &ext["read1 |".len()..];
        let mut parts: Vec<&str> = suffix.split(',').collect();
        parts.sort();
        assert_eq!(parts, vec!["genomeA", "genomeB"]);
    }

    #[test]
    fn test_read_mapped_true_and_false() {
        let mut rm = ResultMap::new();
        assert!(!rm.read_mapped("read1"));
        rm.add_read_map("read1".to_string(), "genomeA".to_string());
        assert!(rm.read_mapped("read1"));
        assert!(!rm.read_mapped("read2"));
    }

    #[test]
    fn test_empty_read_map_clears_all() {
        let mut rm = ResultMap::new();
        rm.add_read_map("read1".to_string(), "genomeA".to_string());
        rm.add_read_map("read2".to_string(), "genomeB".to_string());
        assert!(rm.read_mapped("read1"));
        rm.empty_read_map();
        assert!(!rm.read_mapped("read1"));
        assert!(!rm.read_mapped("read2"));
    }
}
