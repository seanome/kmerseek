use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::SigStore;
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

use crate::uniprot::{self, UniProtEntry};

pub const SEED: u64 = 42;

/// Represents the position of a k-mer in a protein sequence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerPosition {
    pub protein: UniProtEntry,
    pub position: usize, // 0-based position in sequence
}

// TODO: move to the protein level struct
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
    pub hashval: u64,
    pub encoded_kmer: String,
    pub original_kmer_to_position: HashMap<String, Vec<usize>>,
}

// Represents a single protein's k-mer information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerSignature {
    // The protein's signature for searching
    pub signature: SerializableSignature,
    // Hashval -> KmerInfo (encoded -> original k-mer -> positions)
    pub kmer_infos: HashMap<u64, KmerInfo>,
}

impl KmerSignature {
    pub fn new(signature: SmallSignature) -> Self {
        Self { kmer_infos: HashMap::new(), signature: SerializableSignature::from(signature) }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableSignature {
    pub location: String,
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
}

impl From<SmallSignature> for SerializableSignature {
    fn from(sig: SmallSignature) -> Self {
        Self { location: sig.location, name: sig.name, md5sum: sig.md5sum, minhash: sig.minhash }
    }
}

impl From<SigStore> for SerializableSignature {
    fn from(sig: SigStore) -> Self {
        Self {
            location: sig.filename().clone(),
            name: sig.name().clone(),
            md5sum: sig.md5sum().to_string(),
            minhash: sig.minhash().unwrap().clone(),
        }
    }
}
