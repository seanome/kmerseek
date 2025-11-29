use crate::encoding::get_hash_function_from_moltype;
use crate::kmer::KmerInfo;
use crate::signature::StableSignature;
use crate::SEED;
use serde::{Deserialize, Serialize};
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::collections::HashMap;

pub const PROTEIN_TO_MINHASH_RATIO: u32 = 3;

/// ProteinSketches contain Signatures and additional information around K-mer positions and their sequence data
#[derive(Debug, Clone)]
pub struct ProteinSketch {
    name: String,
    signature: StableSignature,
    moltype: String,
    protein_ksize: u32,
    scaled: u32,
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
        state.serialize_field("name", &self.name)?;
        state.serialize_field("signature", &self.signature)?;
        state.serialize_field("moltype", &self.moltype)?;
        state.serialize_field("protein_ksize", &self.protein_ksize)?;
        state.serialize_field("scaled", &self.scaled)?;
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
                let mut name = None;
                let mut signature = None;
                let mut moltype = None;
                let mut protein_ksize = None;
                let mut kmer_infos = None;
                let mut scaled = None;

                while let Some(key) = map.next_key()? {
                    match key {
                        "name" => {
                            if name.is_some() {
                                return Err(de::Error::duplicate_field("name"));
                            }
                            name = Some(map.next_value()?);
                        }
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
                        "scaled" => {
                            if scaled.is_some() {
                                return Err(de::Error::duplicate_field("scaled"));
                            }
                            scaled = Some(map.next_value()?);
                        }
                        "kmer_infos" => {
                            if kmer_infos.is_some() {
                                return Err(de::Error::duplicate_field("kmer_infos"));
                            }
                            kmer_infos = Some(map.next_value()?);
                        }
                        _ => {
                            map.next_value::<de::IgnoredAny>()?;
                        }
                    }
                }

                let name = name.ok_or_else(|| de::Error::missing_field("name"))?;
                let signature = signature.ok_or_else(|| de::Error::missing_field("signature"))?;
                let moltype = moltype.ok_or_else(|| de::Error::missing_field("moltype"))?;
                let protein_ksize =
                    protein_ksize.ok_or_else(|| de::Error::missing_field("protein_ksize"))?;
                let scaled = scaled.ok_or_else(|| de::Error::missing_field("scaled"))?;
                let kmer_infos =
                    kmer_infos.ok_or_else(|| de::Error::missing_field("kmer_infos"))?;

                Ok(ProteinSketch {
                    name,
                    signature,
                    moltype,
                    protein_ksize,
                    scaled,
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
    pub fn new(name: &str, protein_ksize: u32, scaled: u32, moltype: &str) -> anyhow::Result<Self> {
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
            name: name.to_string(),
            signature,
            moltype: moltype.to_string(),
            protein_ksize,
            scaled,
            kmer_infos: HashMap::new(),
            efficient_data: None,
        })
    }

    /// Create a ProteinSketch from existing signature data
    /// This is useful for reconstructing signatures during load operations
    pub fn from_existing_data(
        name: &str,
        signature: StableSignature,
        moltype: String,
        protein_ksize: u32,
        scaled: u32,
        kmer_infos: HashMap<u64, KmerInfo>,
    ) -> Self {
        Self {
            name: name.to_string(),
            signature,
            moltype,
            protein_ksize,
            scaled,
            kmer_infos,
            efficient_data: None,
        }
    }

    /// Create a ProteinSketch from efficient storage data
    pub fn from_efficient_data(
        data: ProteinSketchStore,
        moltype: String,
        protein_ksize: u32,
        scaled: u32,
    ) -> anyhow::Result<Self> {
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
        let name = data.name.clone();

        let signature = StableSignature {
            location: String::new(),
            name: name.clone(),
            md5sum,
            minhash,
            moltype: moltype.clone(),
            ksize: protein_ksize,
        };

        Ok(Self {
            name,
            signature,
            moltype,
            protein_ksize,
            scaled,
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
    pub fn get_moltype_sequence(&self) -> Option<&str> {
        self.efficient_data.as_ref()?.get_encoded_sequence()
    }

    /// Add a protein sequence to the signature
    pub fn add_protein(&mut self, sequence: &[u8]) -> anyhow::Result<()> {
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

    /// Get the minhash k-mer size
    pub fn scaled(&self) -> u32 {
        self.scaled
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

    // Check compatibility between another sketch
    pub fn is_compatible(&self, other: &ProteinSketch) -> bool {
        let same_moltype = self.moltype() == other.moltype();
        let same_scaled = self.scaled() == other.scaled();
        let same_ksize = self.protein_ksize() == other.protein_ksize();
        same_moltype && same_scaled && same_ksize
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
