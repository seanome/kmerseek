use std::collections::HashMap;

use anyhow::Result;
use serde::{Deserialize, Serialize};
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::SigStore;
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

use crate::{encoding::get_hash_function_from_moltype, kmer::KmerInfo};

pub const SEED: u64 = 42;
pub const PROTEIN_TO_MINHASH_RATIO: u32 = 3;

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
#[serde(crate = "serde")]
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
#[serde(crate = "serde")]
pub struct ProteinSignature {
    signature: SerializableSignature,
    moltype: String,
    protein_ksize: u32,
    // Hashval -> KmerInfo (encoded -> original k-mer -> positions)
    kmer_infos: HashMap<u64, KmerInfo>,
    // Efficient storage data (optional, for performance)
    #[serde(skip)]
    efficient_data: Option<ProteinSignatureData>,
}

// Represents a single protein's k-mer signature and hashval -> kmer info mapping
impl ProteinSignature {
    /// Create a new ProteinSignature with the given protein k-mer size
    pub fn new(name: &str, protein_ksize: u32, scaled: u32, moltype: &str) -> Result<Self> {
        let hash_function = get_hash_function_from_moltype(moltype)?;
        let minhash_ksize = protein_ksize * PROTEIN_TO_MINHASH_RATIO; // Convert protein ksize to minhash ksize

        let minhash = KmerMinHash::new(
            scaled,
            minhash_ksize,
            hash_function,
            SEED,
            true, // track_abundance
            0,    // num (use scaled instead)
        );

        let signature = SerializableSignature {
            location: String::new(),
            name: name.to_string(),
            md5sum: String::new(),
            minhash,
        };

        Ok(Self {
            signature,
            moltype: moltype.to_string(),
            protein_ksize,
            kmer_infos: HashMap::new(),
            efficient_data: None,
        })
    }

    /// Create a ProteinSignature from existing signature data
    /// This is useful for reconstructing signatures during load operations
    pub fn from_existing_data(
        signature: SerializableSignature,
        moltype: String,
        protein_ksize: u32,
        kmer_infos: HashMap<u64, KmerInfo>,
    ) -> Self {
        Self { signature, moltype, protein_ksize, kmer_infos, efficient_data: None }
    }

    /// Create a ProteinSignature from efficient storage data
    pub fn from_efficient_data(
        data: ProteinSignatureData,
        moltype: String,
        protein_ksize: u32,
        scaled: u32,
    ) -> Result<Self> {
        let hash_function = get_hash_function_from_moltype(&moltype)?;
        let minhash_ksize = protein_ksize * PROTEIN_TO_MINHASH_RATIO;

        // Reconstruct the minhash from raw data
        let mut minhash = KmerMinHash::new(
            scaled,
            minhash_ksize,
            hash_function,
            SEED,
            true, // track_abundance
            0,    // num (use scaled instead)
        );

        if let Some(abunds) = &data.abunds {
            minhash.add_many_with_abund(
                &data.mins.clone().into_iter().zip(abunds.iter().cloned()).collect::<Vec<_>>(),
            )?;
        } else {
            minhash.add_many(&data.mins)?;
        }

        // Generate md5sum from mins
        let md5sum = data.mins.iter().fold(0u64, |acc, &min| acc.wrapping_add(min));
        let md5sum = format!("{:x}", md5sum);

        let signature = SerializableSignature {
            location: String::new(),
            name: data.name.clone(),
            md5sum,
            minhash,
        };

        Ok(Self {
            signature,
            moltype,
            protein_ksize,
            kmer_infos: data.kmer_infos.clone(),
            efficient_data: Some(data),
        })
    }

    /// Convert to efficient storage format
    pub fn to_efficient_data(&self, include_raw_sequence: bool) -> ProteinSignatureData {
        let minhash = &self.signature.minhash;
        let mins = minhash.mins().to_vec();
        let abunds = minhash.abunds().map(|abunds| abunds.to_vec());

        let raw_sequence = if include_raw_sequence {
            // If we have efficient data with raw sequence, use it
            if let Some(ref data) = self.efficient_data {
                data.raw_sequence.clone()
            } else {
                None // We don't have the raw sequence stored
            }
        } else {
            None
        };

        ProteinSignatureData::new(
            self.signature.name.clone(),
            mins,
            abunds,
            self.kmer_infos.clone(),
            raw_sequence,
        )
    }

    /// Convert to efficient storage format with pre-allocated sequence capacity
    pub fn to_efficient_data_with_capacity(
        &self,
        sequence_capacity: usize,
    ) -> ProteinSignatureData {
        let minhash = &self.signature.minhash;
        let mins = minhash.mins().to_vec();
        let abunds = minhash.abunds().map(|abunds| abunds.to_vec());

        ProteinSignatureData::with_sequence_capacity(
            self.signature.name.clone(),
            mins,
            abunds,
            self.kmer_infos.clone(),
            sequence_capacity,
        )
    }

    /// Set the efficient data (useful for performance optimization)
    pub fn set_efficient_data(&mut self, data: ProteinSignatureData) {
        self.efficient_data = Some(data);
    }

    /// Get the efficient data if available
    pub fn get_efficient_data(&self) -> Option<&ProteinSignatureData> {
        self.efficient_data.as_ref()
    }

    /// Check if efficient data is available
    pub fn has_efficient_data(&self) -> bool {
        self.efficient_data.is_some()
    }

    /// Get the raw sequence if stored in efficient data
    pub fn get_raw_sequence(&self) -> Option<&str> {
        self.efficient_data.as_ref()?.get_raw_sequence()
    }

    /// Add a protein sequence to the signature
    pub fn add_protein(&mut self, sequence: &[u8]) -> Result<()> {
        self.signature.minhash.add_protein(sequence)?;

        // Generate a simple hash-based identifier from the minhash data
        let md5sum =
            self.signature.minhash.mins().iter().fold(0u64, |acc, &min| acc.wrapping_add(min));
        self.signature.md5sum = format!("{:x}", md5sum);

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
        self.protein_ksize * PROTEIN_TO_MINHASH_RATIO
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

/// Efficient storage structure for protein signatures
/// Stores raw values to avoid serialization overhead and optionally includes raw sequences
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(crate = "serde")]
pub struct ProteinSignatureData {
    /// Protein name/identifier
    pub name: String,
    /// MinHash minimum values (hashvals)
    pub mins: Vec<u64>,
    /// MinHash abundance values (if tracked)
    pub abunds: Option<Vec<u64>>,
    /// K-mer information mapping
    pub kmer_infos: HashMap<u64, KmerInfo>,
    /// Raw protein sequence (optional, for performance)
    pub raw_sequence: Option<String>,
}

impl ProteinSignatureData {
    /// Create new signature data with optional raw sequence storage
    pub fn new(
        name: String,
        mins: Vec<u64>,
        abunds: Option<Vec<u64>>,
        kmer_infos: HashMap<u64, KmerInfo>,
        raw_sequence: Option<String>,
    ) -> Self {
        Self { name, mins, abunds, kmer_infos, raw_sequence }
    }

    /// Create new signature data with pre-allocated capacity for raw sequence
    pub fn with_sequence_capacity(
        name: String,
        mins: Vec<u64>,
        abunds: Option<Vec<u64>>,
        kmer_infos: HashMap<u64, KmerInfo>,
        sequence_capacity: usize,
    ) -> Self {
        let raw_sequence = if sequence_capacity > 0 {
            Some(String::with_capacity(sequence_capacity))
        } else {
            None
        };

        Self { name, mins, abunds, kmer_infos, raw_sequence }
    }

    /// Set the raw sequence (only if storage is enabled)
    pub fn set_raw_sequence(&mut self, sequence: String) {
        if self.raw_sequence.is_some() {
            self.raw_sequence = Some(sequence);
        }
    }

    /// Get the raw sequence if stored
    pub fn get_raw_sequence(&self) -> Option<&str> {
        self.raw_sequence.as_deref()
    }

    /// Check if raw sequence storage is enabled
    pub fn has_raw_sequence_storage(&self) -> bool {
        self.raw_sequence.is_some()
    }

    /// Get the number of k-mers in this signature
    pub fn kmer_count(&self) -> usize {
        self.mins.len()
    }

    /// Get the total size of stored data in bytes (approximate)
    pub fn estimated_size(&self) -> usize {
        let mut size = 0;

        // Name size
        size += self.name.len();

        // Mins size (8 bytes per u64)
        size += self.mins.len() * 8;

        // Abunds size if present
        if let Some(ref abunds) = self.abunds {
            size += abunds.len() * 8;
        }

        // Kmer infos size (approximate)
        for kmer_info in self.kmer_infos.values() {
            size += 8; // hashval
            size += 4; // ksize
            size += kmer_info.encoded_kmer.len();
            size += kmer_info.original_kmer_to_position.len() * 16; // rough estimate
        }

        // Raw sequence size if present
        if let Some(ref seq) = self.raw_sequence {
            size += seq.len();
        }

        size
    }
}
