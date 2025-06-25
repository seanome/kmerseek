use anyhow::Result;
use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};

/// Returns the appropriate `HashFunctions` variant for the given `moltype`.
///
/// # Arguments
///
/// * `moltype` - A string slice that specifies the molecule type. Supported values:
///   - `"protein"` or `"raw"` for standard protein encoding
///   - `"hp"` for hydrophobic/polar encoding
///   - `"dayhoff"` for Dayhoff encoding
///
/// # Returns
///
/// * `Ok(HashFunctions)` if the `moltype` is valid
/// * `Err(...)` if the `moltype` is unrecognized
pub fn get_hash_function_from_moltype(moltype: &str) -> Result<HashFunctions, anyhow::Error> {
    match moltype {
        "protein" | "raw" => Ok(HashFunctions::Murmur64Protein),
        "hp" => Ok(HashFunctions::Murmur64Hp),
        "dayhoff" => Ok(HashFunctions::Murmur64Dayhoff),
        _ => Err(anyhow::anyhow!(
            "Invalid moltype: {}, only 'protein', 'hp', or 'dayhoff' are supported",
            moltype
        )),
    }
}

/// Return an amino acid encoding function for a given `moltype` string.
///
/// # Arguments
/// * `moltype` - A string slice that specifies the molecule type. Supported values:
///   - `"protein"` or `"raw"` for standard protein encoding
///   - `"hp"` for hydrophobic/polar encoding
///   - `"dayhoff"` for Dayhoff encoding
///
/// # Returns
///
/// * `Ok(fn(u8) -> u8)` - A function that encodes an amino acid byte according to
///    the specified `moltype`.
/// * `Err(...)` - An error if the `moltype` is unrecognized.
pub fn get_encoding_fn_from_moltype(moltype: &str) -> Result<fn(u8) -> u8, anyhow::Error> {
    match moltype {
        "protein" | "raw" => Ok(|b| b),
        "hp" => Ok(aa_to_hp),
        "dayhoff" => Ok(aa_to_dayhoff),
        _ => Err(anyhow::anyhow!(
            "Invalid moltype: {}, only 'protein', 'hp', or 'dayhoff' are supported",
            moltype
        )),
    }
}

/// Process a k-mer to get its encoded version based on the specified moltype.
///
/// # Arguments
/// * `kmer` - A string slice containing the k-mer to encode
/// * `moltype` - A string slice that specifies the molecule type. Supported values:
///   - `"protein"` or `"raw"` for standard protein encoding
///   - `"hp"` for hydrophobic/polar encoding
///   - `"dayhoff"` for Dayhoff encoding
///
/// # Returns
/// * `Ok((String, String))` - A tuple containing (encoded_kmer, original_kmer)
/// * `Err(...)` - An error if the `moltype` is unrecognized
pub fn encode_kmer(kmer: &str, moltype: &str) -> Result<(String, String)> {
    let encoding_fn = get_encoding_fn_from_moltype(moltype)?;

    let mut encoded = String::with_capacity(kmer.len());
    let mut original = String::with_capacity(kmer.len());

    for &b in kmer.as_bytes() {
        encoded.push(encoding_fn(b) as char);
    }
    original.push_str(kmer);

    Ok((encoded, original))
}

/// Process a k-mer to get its encoded version based on the specified moltype and hash function.
///
/// # Arguments
/// * `kmer` - A string slice containing the k-mer to encode
/// * `encoding_fn` - A function that encodes an amino acid byte according to
///    the specified `moltype`.
///
/// # Returns
/// * `Ok((String, String))` - A tuple containing (encoded_kmer, original_kmer)
/// * `Err(...)` - An error if the `moltype` is unrecognized
pub fn encode_kmer_with_encoding_fn(
    kmer: &str,
    encoding_fn: fn(u8) -> u8,
) -> Result<(String, String)> {
    let mut encoded = String::with_capacity(kmer.len());
    let mut original = String::with_capacity(kmer.len());

    for &b in kmer.as_bytes() {
        encoded.push(encoding_fn(b) as char);
    }
    original.push_str(kmer);

    Ok((encoded, original))
}

#[cfg(test)]
mod tests {
    use super::*;
    use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};

    #[test]
    fn test_get_hash_function_from_moltype() {
        if let Ok(hf) = get_hash_function_from_moltype("protein") {
            assert_eq!(hf, HashFunctions::Murmur64Protein);
        } else {
            panic!("Expected HashFunctions::Murmur64Protein for 'protein'");
        }

        if let Ok(hf) = get_hash_function_from_moltype("raw") {
            assert_eq!(hf, HashFunctions::Murmur64Protein);
        } else {
            panic!("Expected HashFunctions::Murmur64Protein for 'raw'");
        }

        if let Ok(hf) = get_hash_function_from_moltype("hp") {
            assert_eq!(hf, HashFunctions::Murmur64Hp);
        } else {
            panic!("Expected HashFunctions::Murmur64Hp for 'hp'");
        }

        if let Ok(hf) = get_hash_function_from_moltype("dayhoff") {
            assert_eq!(hf, HashFunctions::Murmur64Dayhoff);
        } else {
            panic!("Expected HashFunctions::Murmur64Dayhoff for 'dayhoff'");
        }

        assert!(get_hash_function_from_moltype("xyz").is_err());
    }

    #[test]
    fn test_get_encoding_fn_from_moltype() {
        if let Ok(raw_fn) = get_encoding_fn_from_moltype("raw") {
            assert_eq!(raw_fn(b'A'), b'A');
        } else {
            panic!("Expected identity function for 'raw'");
        }

        if let Ok(protein_fn) = get_encoding_fn_from_moltype("protein") {
            assert_eq!(protein_fn(b'F'), b'F');
        } else {
            panic!("Expected identity function for 'protein'");
        }

        if let Ok(hp_fn) = get_encoding_fn_from_moltype("hp") {
            assert_eq!(hp_fn(b'A'), aa_to_hp(b'A'));
        } else {
            panic!("Expected aa_to_hp for 'hp'");
        }

        if let Ok(dayhoff_fn) = get_encoding_fn_from_moltype("dayhoff") {
            assert_eq!(dayhoff_fn(b'G'), aa_to_dayhoff(b'G'));
        } else {
            panic!("Expected aa_to_dayhoff for 'dayhoff'");
        }

        assert!(get_encoding_fn_from_moltype("xyz").is_err());
    }
}
