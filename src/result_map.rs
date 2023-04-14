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
        //println!("adding read {} and genome {}", read_id, genome_id);

        if let Some(mut genome_set) = self.read_map.get_mut(&read_id) {
            genome_set.insert(genome_id);
        } else {
            let mut new_genome_set = HashSet::new();
            new_genome_set.insert(genome_id);
            self.read_map.insert(read_id, new_genome_set);
        }
    }

    pub fn get_ext_id(&self, read_id: &String) -> String {
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
        read_id.clone() + " |" + &genomes
    }

    pub fn read_mapped(&self, read_id: &String) -> bool {
        self.read_map.contains_key(read_id)
    }

    pub fn empty_read_map(&mut self) {
        self.read_map.clear();
    }
}
