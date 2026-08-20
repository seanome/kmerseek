use anyhow::Result;
use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};

const MURMUR64PROTEIN: &str = "Murmur64Protein";
const MURMUR64DAYHOFF: &str = "Murmur64Dayhoff";
const MURMUR64HP: &str = "Murmur64Hp";

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
        // Custom HP alphabets pre-encode the sequence before hashing, so the
        // hash function sees an already-encoded h/p sequence and uses identity.
        s if s.starts_with("hp_") => Ok(HashFunctions::Murmur64Protein),
        _ => Err(anyhow::anyhow!(
            "Invalid moltype: {}. Supported values: 'protein', 'dayhoff', 'hp', \
             'hp_lehninger', 'hp_thomas_dill', 'hp_kyte_doolittle', \
             'hp_thomas_dill_no_c', 'hp_lehninger_c_nonpolar', 'hp_lehninger_hpc', 'hp_pbotc_1st_ed', \
             'hp_shuffled_control', 'hp_shuffled_control_1'..'hp_shuffled_control_10'",
            moltype
        )),
    }
}

pub fn get_moltype_from_hash_function_string(
    hash_function: String,
) -> Result<String, anyhow::Error> {
    match hash_function.as_str() {
        MURMUR64PROTEIN => Ok("protein".to_string()),
        MURMUR64DAYHOFF => Ok("dayhoff".to_string()),
        MURMUR64HP => Ok("hp".to_string()),
        _ => Err(anyhow::anyhow!(
            "Invalid hash function: {}, only 'Murmur64' with 'protein', 'dayhoff', or 'hp' are supported", hash_function
        ))
    }
}

pub fn get_moltype_from_hash_function(
    hash_function: HashFunctions,
) -> Result<String, anyhow::Error> {
    match hash_function {
        HashFunctions::Murmur64Protein => Ok("protein".to_string()),
        HashFunctions::Murmur64Dayhoff => Ok("dayhoff".to_string()),
        HashFunctions::Murmur64Hp => Ok("hp".to_string()),
        _ => Err(anyhow::anyhow!(
            "Invalid hash function: {}, only Sourmash HashFunctions::Murmur64 with 'protein', 'dayhoff', or 'hp' are supported", hash_function
        ))
    }
}

#[allow(clippy::doc_overindented_list_items)]
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
        // Custom HP alphabets pre-encode in HpAlphabet::table(); return identity here so
        // callers that only need a fn(u8)->u8 don't crash. process_kmers handles them separately.
        s if s.starts_with("hp_") => Ok(|b| b),
        _ => Err(anyhow::anyhow!(
            "Invalid moltype: {}, only 'protein', 'hp', or 'dayhoff' are supported",
            moltype
        )),
    }
}

/// Encode a sequence using the specified moltype.
///
/// This function applies the encoding function to each amino acid in the sequence,
/// producing an encoded sequence string. The encoded sequence is pre-allocated
/// with the exact capacity needed.
///
/// # Arguments
/// * `sequence` - The sequence to encode (can be a k-mer or full sequence)
/// * `moltype` - The molecule type encoding to use ("protein", "hp", or "dayhoff")
///
/// # Returns
/// * `Ok(String)` - The encoded sequence
/// * `Err(...)` - An error if the moltype is invalid
///
/// # Example
/// ```
/// use kmerseek::encoding::encode_by_moltype;
///
/// let sequence = "MKTAYIAKQR";
/// let encoded = encode_by_moltype(sequence, "hp").unwrap();
/// // encoded will be the HP-encoded version
/// ```
pub fn encode_by_moltype(sequence: &str, moltype: &str) -> Result<String> {
    let encoding_fn = get_encoding_fn_from_moltype(moltype)?;
    encode_with_fn(sequence, encoding_fn)
}

/// Encode a sequence into a molecular type using the provided encoding function.
///
/// This function applies the encoding function to each amino acid in the sequence,
/// producing an encoded sequence string. The encoded sequence is pre-allocated
/// with the exact capacity needed.
///
/// # Arguments
/// * `sequence` - The sequence to encode (can be a k-mer or full sequence)
/// * `encoding_fn` - A function that encodes an amino acid byte according to
///   the specified moltype
///
/// # Returns
/// * `Ok(String)` - The encoded sequence
///
/// # Example
/// ```
/// use kmerseek::encoding::encode_with_fn;
/// use sourmash::encodings::aa_to_hp;
///
/// let sequence = "MKTAYIAKQR";
/// let encoded = encode_with_fn(sequence, aa_to_hp).unwrap();
/// // encoded will be the HP-encoded version
/// ```
pub fn encode_with_fn(sequence: &str, encoding_fn: fn(u8) -> u8) -> Result<String> {
    // WHY: Pre-allocate with exact capacity since the encoded sequence length
    // will always equal the input sequence length. This avoids unnecessary
    // reallocations during encoding.
    let mut encoded = String::with_capacity(sequence.len());

    for &b in sequence.as_bytes() {
        encoded.push(encoding_fn(b) as char);
    }

    Ok(encoded)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::test_fixtures::TEST_KMER;
    use anyhow::Result;
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

    #[test]
    fn test_encode_by_moltype_protein() -> Result<()> {
        let encoded = encode_by_moltype(TEST_KMER, "protein")?;
        assert_eq!(encoded, TEST_KMER);
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_dayhoff() -> Result<()> {
        let encoded = encode_by_moltype(TEST_KMER, "dayhoff")?;
        assert_eq!(encoded, "eeeecbbeeec");
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_hp() -> Result<()> {
        let encoded = encode_by_moltype(TEST_KMER, "hp")?;
        assert_eq!(encoded, "hhhhphhhhhp");
        Ok(())
    }

    #[test]
    fn test_encode_with_fn_protein() -> Result<()> {
        let encoded = encode_with_fn(TEST_KMER, |b| b)?;
        assert_eq!(encoded, TEST_KMER);
        Ok(())
    }

    #[test]
    fn test_encode_with_fn_dayhoff() -> Result<()> {
        let encoded = encode_with_fn(TEST_KMER, aa_to_dayhoff)?;
        assert_eq!(encoded, "eeeecbbeeec");
        Ok(())
    }

    #[test]
    fn test_encode_with_fn_hp() -> Result<()> {
        let encoded = encode_with_fn(TEST_KMER, aa_to_hp)?;
        assert_eq!(encoded, "hhhhphhhhhp");
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_sequence() -> Result<()> {
        let sequence = "MKTAYIAKQR";
        let encoded = encode_by_moltype(sequence, "protein")?;
        assert_eq!(encoded, sequence);
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_sequence_hp() -> Result<()> {
        let sequence = "MKTAYIAKQR";
        let encoded = encode_by_moltype(sequence, "hp")?;
        // Verify the encoding produces the correct length and uses only 'h' and 'p'
        assert_eq!(encoded.len(), sequence.len());
        assert!(encoded.chars().all(|c| c == 'h' || c == 'p'));
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_sequence_dayhoff() -> Result<()> {
        let sequence = "MKTAYIAKQR";
        let encoded = encode_by_moltype(sequence, "dayhoff")?;
        // Dayhoff encoding should produce a different string
        assert_ne!(encoded, sequence);
        assert_eq!(encoded.len(), sequence.len());
        Ok(())
    }
}
