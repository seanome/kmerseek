use crate::uniprot::UniProtEntry;
use serde::{Deserialize, Serialize};

use super::kmer_signature::KmerSignature;

// High-level mapping of a protein sequence to its Sourmash Signature, UniProt features, and k-mer information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Protein {
    pub uniprot_entry: UniProtEntry,
    pub signature_kmers: KmerSignature,
}

impl Protein {
    pub fn new(uniprot_entry: UniProtEntry, signature_kmers: KmerSignature) -> Self {
        Self { uniprot_entry, signature_kmers }
    }
}
