use anyhow::Result;
use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};

use crate::index::ProteomeIndex;

fn get_hash_function_from_moltype(
    moltype: &str,
) -> Result<HashFunctions, std::result::Result<ProteomeIndex, anyhow::Error>> {
    let hash_function = match moltype {
        "protein" | "raw" => HashFunctions::Murmur64Protein,
        "hp" => HashFunctions::Murmur64Hp,
        "dayhoff" => HashFunctions::Murmur64Dayhoff,
        _ => {
            return Err(Err(anyhow::anyhow!(
                "Invalid moltype: {}, only 'protein', 'hp', or 'dayhoff' are supported",
                moltype
            )))
        }
    };
    Ok(hash_function)
}

fn get_encoding_fn_from_moltype(
    moltype: &str,
) -> Result<fn(u8) -> u8, std::result::Result<ProteomeIndex, anyhow::Error>> {
    let encoding_fn = match moltype {
        "protein" | "raw" => |b| b, // raw encoding - just returns the same amino acid
        "hp" => aa_to_hp,
        "dayhoff" => aa_to_dayhoff,
        _ => {
            return Err(Err(anyhow::anyhow!(
                "Invalid moltype: {}, only 'protein', 'hp', or 'dayhoff' are supported",
                moltype
            )))
        }
    };
    Ok(encoding_fn)
}
