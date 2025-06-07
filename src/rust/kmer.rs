use crate::uniprot::{self, UniProtEntry};
use serde::{Deserialize, Serialize};

/// Represents the position of a k-mer in a protein sequence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerPosition {
    pub protein: UniProtEntry,
    pub position: usize, // 0-based position in sequence
}

impl KmerPosition {
    /// Get all features that overlap with this k-mer
    pub fn overlapping_features(&self, kmer_length: usize) -> Vec<&uniprot::ProteinFeature> {
        let kmer_start = self.position + 1; // Convert to 1-based position
        let kmer_end = kmer_start + kmer_length - 1;

        self.protein
            .features
            .iter()
            .filter(|feature| {
                // Check if the feature overlaps with the k-mer
                !(feature.end < kmer_start || feature.start > kmer_end)
            })
            .collect()
    }
}

/// Represents information about a k-mer occurrence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerInfo {
    pub original_kmer: String,
    pub positions: Vec<KmerPosition>,
}

/// Statistics for k-mer frequency analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerStats {
    pub idf: f64,       // Inverse document frequency
    pub frequency: f64, // Raw frequency
}
