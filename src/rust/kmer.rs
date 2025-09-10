use serde::{Deserialize, Serialize};

/// Represents information about a k-mer occurrence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerInfo {
    pub ksize: usize,
    pub hashval: u64,
    pub encoded_kmer: String,
    pub positions: Vec<usize>,
}

impl KmerInfo {
    /// Creates a new KmerInfo with pre-allocated capacity
    pub fn new(hashval: u64, ksize: usize) -> Self {
        Self { ksize, hashval, encoded_kmer: String::with_capacity(ksize), positions: Vec::new() }
    }

    /// Adds a k-mer position
    pub fn add_position(&mut self, position: usize) {
        self.positions.push(position);
    }

    /// Get the number of k-mer occurrences
    pub fn total_occurrences(&self) -> usize {
        self.positions.len()
    }

    /// Check if this k-mer appears at a specific position
    pub fn has_position(&self, position: usize) -> bool {
        self.positions.contains(&position)
    }

    /// Get all positions where this k-mer occurs
    pub fn get_positions(&self) -> &[usize] {
        &self.positions
    }
}
