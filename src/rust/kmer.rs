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
            original_kmer_to_position: HashMap::with_capacity(ksize),
        }
    }

    /// Adds a k-mer position with a pre-allocated string
    pub fn add_kmer_position(&mut self, kmer: &str, position: usize) {
        self.original_kmer_to_position
            .entry(kmer.to_string())
            .or_insert_with(|| Vec::with_capacity(1))
            .push(position);
    }

    /// Get the number of unique original k-mers
    pub fn unique_kmer_count(&self) -> usize {
        self.original_kmer_to_position.len()
    }

    /// Get the total number of k-mer occurrences
    pub fn total_occurrences(&self) -> usize {
        self.original_kmer_to_position.values().map(|positions| positions.len()).sum()
    }

    /// Check if this k-mer appears at a specific position
    pub fn has_position(&self, position: usize) -> bool {
        self.original_kmer_to_position.values().any(|positions| positions.contains(&position))
    }
}
