use anyhow::Result;
use serde;
use serde::{Deserialize, Serialize};
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::SigStore;
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

use crate::encoding::{
    get_hash_function_from_moltype, get_moltype_from_hash_function,
    get_moltype_from_hash_function_string,
};

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
        // WHY: We serialize 8 fields: location, name, md5sum, mins, abunds, scaled, ksize, moltype.
        // The field count must match the actual number of fields serialized, otherwise serde
        // will fail during serialization. This matches the FIELDS constant used in deserialization.
        let mut state = serializer.serialize_struct("StableSignature", 8)?;
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
