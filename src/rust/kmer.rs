use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Represents information about a k-mer occurrence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerInfo {
    pub ksize: usize,
    pub hashval: u64,
    pub encoded_kmer: String,
    pub original_kmer_to_position: HashMap<String, Vec<usize>>,
}

impl KmerInfo {
    /// Creates a new KmerInfo with pre-allocated capacity
    pub fn new(hashval: u64, ksize: usize) -> Self {
        Self {
            ksize,
            hashval,
            encoded_kmer: String::with_capacity(ksize),
            original_kmer_to_position: HashMap::new(),
        }
    }

    /// Adds a k-mer position
    pub fn add_position(&mut self, kmer: &str, position: usize) {
        self.original_kmer_to_position
            .entry(kmer.to_string())
            .or_insert_with(|| Vec::with_capacity(1))
            .push(position);
    }

    /// Get the number of unique original k-mers
    pub fn unique_kmer_count(&self) -> usize {
        self.original_kmer_to_position.len()
    }

    /// Get the number of k-mer occurrences
    pub fn total_occurrences(&self) -> usize {
        self.original_kmer_to_position.values().map(|positions| positions.len()).sum()
    }

    /// Check if this k-mer appears at a specific position
    pub fn has_position(&self, position: usize) -> bool {
        self.original_kmer_to_position.values().any(|positions| positions.contains(&position))
    }

    /// Get all positions where this k-mer occurs
    pub fn get_positions(&self) -> Vec<usize> {
        self.original_kmer_to_position.values().flatten().cloned().collect()
    }
}
