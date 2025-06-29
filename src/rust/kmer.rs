use std::collections::HashMap;

use serde::{Deserialize, Serialize};

/// Represents information about a k-mer occurrence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerInfo {
    pub ksize: usize,
    pub hashval: u64,
    pub encoded_kmer: String,
    pub original_kmer_to_position: HashMap<String, Vec<usize>>,
}

impl KmerInfo {
    /// Creates a new KmerInfo with pre-allocated strings for k-mers
    pub fn new(hashval: u64, ksize: usize) -> Self {
        Self {
            ksize,
            hashval,
            encoded_kmer: String::with_capacity(ksize),
            original_kmer_to_position: HashMap::new(),
        }
    }

    /// Adds a k-mer position with a pre-allocated string
    pub fn add_kmer_position(&mut self, kmer: &str, position: usize) {
        let mut key = String::with_capacity(self.ksize);
        key.push_str(kmer);
        self.original_kmer_to_position.entry(key).or_insert_with(Vec::new).push(position);
    }
}
