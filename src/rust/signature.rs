use std::collections::HashMap;

use anyhow::Result;
use serde::{Deserialize, Serialize};
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::SigStore;
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

use crate::{encoding::get_hash_function_from_moltype, kmer::KmerInfo};

pub const SEED: u64 = 42;

/// Trait for accessing signature information
pub trait SignatureAccess {
    fn get_minhash(&self) -> &KmerMinHash;
    fn get_name(&self) -> &str;
    fn get_location(&self) -> &str;
    fn get_md5sum(&self) -> &str;
}

impl SignatureAccess for SmallSignature {
    fn get_minhash(&self) -> &KmerMinHash {
        &self.minhash
    }

    fn get_name(&self) -> &str {
        &self.name
    }

    fn get_location(&self) -> &str {
        &self.location
    }

    fn get_md5sum(&self) -> &str {
        &self.md5sum
    }
}

impl SignatureAccess for SerializableSignature {
    fn get_minhash(&self) -> &KmerMinHash {
        &self.minhash
    }

    fn get_name(&self) -> &str {
        &self.name
    }

    fn get_location(&self) -> &str {
        &self.location
    }

    fn get_md5sum(&self) -> &str {
        &self.md5sum
    }
}

impl SignatureAccess for &SerializableSignature {
    fn get_minhash(&self) -> &KmerMinHash {
        &self.minhash
    }

    fn get_name(&self) -> &str {
        &self.name
    }

    fn get_location(&self) -> &str {
        &self.location
    }

    fn get_md5sum(&self) -> &str {
        &self.md5sum
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

/// A wrapper around SmallSignature that handles protein k-mer size conversions
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinSignature {
    signature: SerializableSignature,
    moltype: String,
    protein_ksize: u32,
    // Hashval -> KmerInfo (encoded -> original k-mer -> positions)
    kmer_infos: HashMap<u64, KmerInfo>,
}

// Represents a single protein's k-mer signature and hashval -> kmer info mapping
impl ProteinSignature {
    /// Create a new ProteinSignature with the given protein k-mer size
    pub fn new(
        name: &str,
        protein_ksize: u32,
        scaled: u32,
        moltype: &str,
        seed: u64,
    ) -> Result<Self> {
        let hash_function = get_hash_function_from_moltype(moltype)?;
        let minhash_ksize = protein_ksize * 3; // Convert protein ksize to minhash ksize

        let minhash = KmerMinHash::new(
            scaled,
            minhash_ksize,
            hash_function,
            seed,
            true, // track_abundance
            0,    // num (use scaled instead)
        );

        let signature = SerializableSignature {
            location: "".to_string(),
            name: name.to_string(),
            md5sum: "".to_string(),
            minhash,
        };

        Ok(Self {
            signature,
            moltype: moltype.to_string(),
            protein_ksize,
            kmer_infos: HashMap::new(),
        })
    }

    /// Add a protein sequence to the signature
    pub fn add_protein(&mut self, sequence: &[u8]) -> Result<()> {
        self.signature.minhash.add_protein(sequence)?;

        // Generate a simple hash-based identifier from the minhash data
        let mins = self.signature.minhash.mins();
        let mut hash = 0u64;
        for min in mins {
            hash = hash.wrapping_add(min);
        }
        self.signature.md5sum = format!("{:x}", hash);

        Ok(())
    }

    /// Get the protein k-mer size
    pub fn protein_ksize(&self) -> u32 {
        self.protein_ksize
    }

    /// Get the moltype
    pub fn moltype(&self) -> &str {
        &self.moltype
    }

    /// Get the minhash k-mer size
    pub fn minhash_ksize(&self) -> u32 {
        self.protein_ksize * 3
    }

    /// Get the underlying SmallSignature
    pub fn into_signature(self) -> SerializableSignature {
        self.signature
    }

    /// Get a reference to the underlying SmallSignature
    pub fn signature(&self) -> &SerializableSignature {
        &self.signature
    }

    /// Get a reference to the kmer infos HashMap (hashval -> kmer info)
    pub fn kmer_infos(&self) -> &HashMap<u64, KmerInfo> {
        &self.kmer_infos
    }

    /// Get a mutable reference to the kmer infos HashMap (hashval -> kmer info)
    pub fn kmer_infos_mut(&mut self) -> &mut HashMap<u64, KmerInfo> {
        &mut self.kmer_infos
    }
}
