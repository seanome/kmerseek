use std::collections::HashMap;

use anyhow::Result;
use serde;
use serde::{Deserialize, Serialize};
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::SigStore;
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

use crate::{
    encoding::{
        get_hash_function_from_moltype, get_moltype_from_hash_function,
        get_moltype_from_hash_function_string,
    },
    kmer::KmerInfo,
};

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

impl SignatureAccess for StableSignature {
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

impl SignatureAccess for &StableSignature {
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

#[derive(Debug, Clone, Default)]
pub struct StableSignature {
    pub location: String,
    pub name: String,
    pub md5sum: String,
    pub minhash: KmerMinHash,
    pub moltype: String,
    pub ksize: u32,
}

// Custom serialization for StableSignature to avoid KmerMinHash serialization issues
impl Serialize for StableSignature {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("StableSignature", 7)?;
        state.serialize_field("location", &self.location)?;
        state.serialize_field("name", &self.name)?;
        state.serialize_field("md5sum", &self.md5sum)?;
        state.serialize_field("mins", &self.minhash.mins())?;
        state.serialize_field("abunds", &self.minhash.abunds())?;
        state.serialize_field("scaled", &self.minhash.scaled())?;
        state.serialize_field("ksize", &self.minhash.ksize())?;
        state.serialize_field("moltype", &self.moltype)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for StableSignature {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, MapAccess, Visitor};
        use std::fmt;

        struct StableSignatureVisitor;

        impl<'de> Visitor<'de> for StableSignatureVisitor {
            type Value = StableSignature;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct StableSignature")
            }

            fn visit_map<V>(self, mut map: V) -> Result<StableSignature, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut location = None;
                let mut name = None;
                let mut md5sum = None;
                let mut mins = None;
                let mut abunds = None;
                let mut scaled = None;
                let mut ksize = None;
                let mut moltype = None;

                while let Some(key) = map.next_key()? {
                    match key {
                        "location" => {
                            if location.is_some() {
                                return Err(de::Error::duplicate_field("location"));
                            }
                            location = Some(map.next_value()?);
                        }
                        "name" => {
                            if name.is_some() {
                                return Err(de::Error::duplicate_field("name"));
                            }
                            name = Some(map.next_value()?);
                        }
                        "md5sum" => {
                            if md5sum.is_some() {
                                return Err(de::Error::duplicate_field("md5sum"));
                            }
                            md5sum = Some(map.next_value()?);
                        }
                        "mins" => {
                            if mins.is_some() {
                                return Err(de::Error::duplicate_field("mins"));
                            }
                            mins = Some(map.next_value()?);
                        }
                        "abunds" => {
                            if abunds.is_some() {
                                return Err(de::Error::duplicate_field("abunds"));
                            }
                            abunds = Some(map.next_value()?);
                        }
                        "scaled" => {
                            if scaled.is_some() {
                                return Err(de::Error::duplicate_field("scaled"));
                            }
                            scaled = Some(map.next_value()?);
                        }
                        "ksize" => {
                            if ksize.is_some() {
                                return Err(de::Error::duplicate_field("ksize"));
                            }
                            ksize = Some(map.next_value()?);
                        }
                        "moltype" => {
                            if moltype.is_some() {
                                return Err(de::Error::duplicate_field("moltype"));
                            }
                            moltype = Some(map.next_value::<String>()?);
                        }
                        _ => {
                            let _ = map.next_value::<de::IgnoredAny>()?;
                        }
                    }
                }

                let location = location.ok_or_else(|| de::Error::missing_field("location"))?;
                let name = name.ok_or_else(|| de::Error::missing_field("name"))?;
                let md5sum = md5sum.ok_or_else(|| de::Error::missing_field("md5sum"))?;
                let mins: Vec<u64> = mins.ok_or_else(|| de::Error::missing_field("mins"))?;
                let abunds: Option<Vec<u64>> =
                    abunds.ok_or_else(|| de::Error::missing_field("abunds"))?;
                let scaled: u32 = scaled.ok_or_else(|| de::Error::missing_field("scaled"))?;
                let ksize: u32 = ksize.ok_or_else(|| de::Error::missing_field("ksize"))?;
                let moltype = moltype.ok_or_else(|| de::Error::missing_field("moltype"))?;

                let hash_function =
                    get_hash_function_from_moltype(&moltype).map_err(de::Error::custom)?;
                // Reconstruct KmerMinHash from the stored data with correct parameters
                // We'll use a default hash function since we don't store the specific one
                let mut minhash = KmerMinHash::new(
                    scaled,
                    ksize,
                    hash_function,
                    SEED, // seed
                    true, // track_abundance
                    0,    // num
                );

                // Add the stored data
                if let Some(abunds) = abunds {
                    let data: Vec<(u64, u64)> = mins.into_iter().zip(abunds).collect();
                    minhash.add_many_with_abund(&data).map_err(de::Error::custom)?;
                } else {
                    minhash.add_many(&mins).map_err(de::Error::custom)?;
                }

                Ok(StableSignature { location, name, md5sum, minhash, moltype, ksize })
            }
        }

        const FIELDS: &[&str] =
            &["location", "name", "md5sum", "mins", "abunds", "scaled", "ksize", "moltype"];
        deserializer.deserialize_struct("StableSignature", FIELDS, StableSignatureVisitor)
    }
}

impl From<SmallSignature> for StableSignature {
    fn from(sig: SmallSignature) -> Self {
        let moltype = get_moltype_from_hash_function(sig.minhash.hash_function())
            .expect("Invalid hash function");
        let ksize = sig.minhash.ksize();
        Self {
            location: sig.location,
            name: sig.name,
            md5sum: sig.md5sum,
            minhash: sig.minhash,
            moltype: moltype,
            ksize: ksize as u32,
        }
    }
}

impl From<SigStore> for StableSignature {
    fn from(sig: SigStore) -> Self {
        let moltype = get_moltype_from_hash_function_string(sig.hash_function());
        Self {
            location: sig.filename().clone(),
            name: sig.name().clone(),
            md5sum: sig.md5sum().to_string(),
            minhash: sig.minhash().unwrap().clone(),
            moltype: moltype.unwrap(),
            ksize: sig.minhash().unwrap().ksize() as u32,
        }
    }
}

/// A wrapper around `StableSignature` that handles protein k-mer size conversions
#[derive(Debug, Clone)]
pub struct ProteinSketch {
    signature: StableSignature,
    moltype: String,
    protein_ksize: u32,
    // Hashval -> KmerInfo (encoded -> original k-mer -> positions)
    kmer_infos: HashMap<u64, KmerInfo>,
    // Efficient storage data (optional, for performance)
    efficient_data: Option<ProteinSketchStore>,
}

// Custom serialization for ProteinSketch
impl Serialize for ProteinSketch {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("ProteinSketch", 4)?;
        state.serialize_field("signature", &self.signature)?;
        state.serialize_field("moltype", &self.moltype)?;
        state.serialize_field("protein_ksize", &self.protein_ksize)?;
        state.serialize_field("kmer_infos", &self.kmer_infos)?;
        // Skip efficient_data as it's marked with #[serde(skip)]
        state.end()
    }
}

impl<'de> Deserialize<'de> for ProteinSketch {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, MapAccess, Visitor};
        use std::fmt;

        struct ProteinSketchVisitor;

        impl<'de> Visitor<'de> for ProteinSketchVisitor {
            type Value = ProteinSketch;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct ProteinSketch")
            }

            fn visit_map<V>(self, mut map: V) -> Result<ProteinSketch, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut signature = None;
                let mut moltype = None;
                let mut protein_ksize = None;
                let mut kmer_infos = None;

                while let Some(key) = map.next_key()? {
                    match key {
                        "signature" => {
                            if signature.is_some() {
                                return Err(de::Error::duplicate_field("signature"));
                            }
                            signature = Some(map.next_value()?);
                        }
                        "moltype" => {
                            if moltype.is_some() {
                                return Err(de::Error::duplicate_field("moltype"));
                            }
                            moltype = Some(map.next_value()?);
                        }
                        "protein_ksize" => {
                            if protein_ksize.is_some() {
                                return Err(de::Error::duplicate_field("protein_ksize"));
                            }
                            protein_ksize = Some(map.next_value()?);
                        }
                        "kmer_infos" => {
                            if kmer_infos.is_some() {
                                return Err(de::Error::duplicate_field("kmer_infos"));
                            }
                            kmer_infos = Some(map.next_value()?);
                        }
                        _ => {
                            let _ = map.next_value::<de::IgnoredAny>()?;
                        }
                    }
                }

                let signature = signature.ok_or_else(|| de::Error::missing_field("signature"))?;
                let moltype = moltype.ok_or_else(|| de::Error::missing_field("moltype"))?;
                let protein_ksize =
                    protein_ksize.ok_or_else(|| de::Error::missing_field("protein_ksize"))?;
                let kmer_infos =
                    kmer_infos.ok_or_else(|| de::Error::missing_field("kmer_infos"))?;

                Ok(ProteinSketch {
                    signature,
                    moltype,
                    protein_ksize,
                    kmer_infos,
                    efficient_data: None, // Always None when deserializing
                })
            }
        }

        const FIELDS: &[&str] = &["signature", "moltype", "protein_ksize", "kmer_infos"];
        deserializer.deserialize_struct("ProteinSketch", FIELDS, ProteinSketchVisitor)
    }
}

// Represents a single protein's k-mer signature and hashval -> kmer info mapping
impl ProteinSketch {
    /// Create a new ProteinSketch with the given protein k-mer size
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

        let signature = StableSignature {
            location: String::new(),
            name: name.to_string(),
            md5sum: String::new(),
            minhash,
            moltype: moltype.to_string(),
            ksize: protein_ksize,
        };

        Ok(Self {
            signature,
            moltype: moltype.to_string(),
            protein_ksize,
            kmer_infos: HashMap::new(),
            efficient_data: None,
        })
    }

    /// Create a ProteinSketch from existing signature data
    /// This is useful for reconstructing signatures during load operations
    pub fn from_existing_data(
        signature: StableSignature,
        moltype: String,
        protein_ksize: u32,
        kmer_infos: HashMap<u64, KmerInfo>,
    ) -> Self {
        Self { signature, moltype, protein_ksize, kmer_infos, efficient_data: None }
    }

    /// Create a ProteinSketch from efficient storage data
    pub fn from_efficient_data(
        data: ProteinSketchStore,
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

        let signature = StableSignature {
            location: String::new(),
            name: data.name.clone(),
            md5sum,
            minhash,
            moltype: moltype.clone(),
            ksize: protein_ksize,
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
    pub fn to_efficient_data(&self, include_raw_sequence: bool) -> ProteinSketchStore {
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

        // Get encoded sequence if available
        let encoded_sequence = if let Some(ref data) = self.efficient_data {
            data.encoded_sequence.clone()
        } else {
            None
        };

        ProteinSketchStore::new(
            self.signature.name.clone(),
            mins,
            abunds,
            self.kmer_infos.clone(),
            raw_sequence,
            encoded_sequence,
        )
    }

    /// Convert to efficient storage format with pre-allocated sequence capacity
    pub fn to_efficient_data_with_capacity(&self, sequence_capacity: usize) -> ProteinSketchStore {
        let minhash = &self.signature.minhash;
        let mins = minhash.mins().to_vec();
        let abunds = minhash.abunds().map(|abunds| abunds.to_vec());

        ProteinSketchStore::with_sequence_capacity(
            self.signature.name.clone(),
            mins,
            abunds,
            self.kmer_infos.clone(),
            sequence_capacity,
        )
    }

    /// Set the efficient data (useful for performance optimization)
    pub fn set_efficient_data(&mut self, data: ProteinSketchStore) {
        self.efficient_data = Some(data);
    }

    /// Get the efficient data if available
    pub fn get_efficient_data(&self) -> Option<&ProteinSketchStore> {
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

    /// Get the encoded sequence if available
    pub fn get_encoded_sequence(&self) -> Option<&str> {
        self.efficient_data.as_ref()?.get_encoded_sequence()
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

    /// Get the underlying `StableSignature`
    pub fn into_signature(self) -> StableSignature {
        self.signature
    }

    /// Get a reference to the underlying `StableSignature`
    pub fn signature(&self) -> &StableSignature {
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

/// Efficient storage structure for protein sketches
/// Stores raw values to avoid serialization overhead and optionally includes raw sequences
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(crate = "serde")]
pub struct ProteinSketchStore {
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
    /// Encoded sequence (hp/dayhoff/protein) - None for protein encoding to save space
    pub encoded_sequence: Option<String>,
}

impl ProteinSketchStore {
    /// Create new signature data with optional raw sequence storage
    pub fn new(
        name: String,
        mins: Vec<u64>,
        abunds: Option<Vec<u64>>,
        kmer_infos: HashMap<u64, KmerInfo>,
        raw_sequence: Option<String>,
        encoded_sequence: Option<String>,
    ) -> Self {
        Self { name, mins, abunds, kmer_infos, raw_sequence, encoded_sequence }
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

        Self { name, mins, abunds, kmer_infos, raw_sequence, encoded_sequence: None }
    }

    /// Set the raw sequence (only if storage is enabled)
    pub fn set_raw_sequence(&mut self, sequence: String) {
        if self.raw_sequence.is_some() {
            self.raw_sequence = Some(sequence);
        }
    }

    /// Set the encoded sequence
    pub fn set_encoded_sequence(&mut self, sequence: String) {
        self.encoded_sequence = Some(sequence);
    }

    /// Get the raw sequence if stored
    pub fn get_raw_sequence(&self) -> Option<&str> {
        self.raw_sequence.as_deref()
    }

    /// Get the encoded sequence if stored
    pub fn get_encoded_sequence(&self) -> Option<&str> {
        self.encoded_sequence.as_deref()
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
            // Calculate size for the HashMap structure
            for (kmer, positions) in &kmer_info.original_kmer_to_position {
                size += kmer.len(); // original k-mer string
                size += positions.len() * 8; // positions
            }
        }

        // Raw sequence size if present
        if let Some(ref seq) = self.raw_sequence {
            size += seq.len();
        }

        size
    }
}
