use anyhow::Result;
use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};

use crate::alphabets::{alphabet_table, canonical_moltype};

const MURMUR64PROTEIN: &str = "Murmur64Protein";
const MURMUR64DAYHOFF: &str = "Murmur64Dayhoff";
const MURMUR64HP: &str = "Murmur64Hp";

/// Returns the appropriate `HashFunctions` variant for the given `moltype`.
///
/// # Arguments
///
/// * `moltype` - A string slice that specifies the molecule type. Supported values:
///   - `"protein20"` for the full 20-letter alphabet
///   - `"dayhoff6"` for the Dayhoff alphabet
///   - `"hp_<name><size>"` for an HP table, e.g. `"hp_thomas_dill2"`
///   - `"<name><size>"` for a multi-letter reduced alphabet, e.g. `"sdm12"`
///
/// # Returns
///
/// * `Ok(HashFunctions)` if the `moltype` is valid
/// * `Err(...)` if the `moltype` is unrecognized
pub fn get_hash_function_from_moltype(moltype: &str) -> Result<HashFunctions, anyhow::Error> {
    match canonical_moltype(moltype) {
        "protein20" => Ok(HashFunctions::Murmur64Protein),
        // Dayhoff has no table of ours; sourmash encodes and hashes it. This arm has to come
        // before the table lookup below, which assumes a pre-encoded sequence.
        "dayhoff6" => Ok(HashFunctions::Murmur64Dayhoff),
        // Lehninger is sourmash's own aa_to_hp partition, so sourmash encodes and hashes it
        // and sketches are byte-identical to sourmash's.
        "hp_lehninger2" => Ok(HashFunctions::Murmur64Hp),
        // Every other alphabet pre-encodes through its own table, so the hash function sees
        // an already-encoded sequence and hashes it as protein.
        s if alphabet_table(s).is_some() => Ok(HashFunctions::Murmur64Protein),
        _ => Err(anyhow::anyhow!(
            "Invalid alphabet: {}. Supported values: 'protein20', 'dayhoff6', \
             'hp_lehninger2', 'hp_thomas_dill2', 'hp_kyte_doolittle2', \
             'hp_thomas_dill_no_c2', 'hp_lehninger_c_nonpolar2', 'hp_lehninger_hpc3', \
             'hp_pbotc_1st_ed2', 'hp_random_control2', \
             'hp_random_control2_1'..'hp_random_control2_10', \
             'gbmr4', 'wwmj5', 'gbmr7', 'sdm12', 'mmseqs12', 'wass14', 'hsdm17', \
             'uniprot18'",
            moltype
        )),
    }
}

pub fn get_moltype_from_hash_function_string(
    hash_function: String,
) -> Result<String, anyhow::Error> {
    match hash_function.as_str() {
        MURMUR64PROTEIN => Ok("protein20".to_string()),
        MURMUR64DAYHOFF => Ok("dayhoff6".to_string()),
        MURMUR64HP => Ok("hp_lehninger2".to_string()),
        _ => Err(anyhow::anyhow!(
            "Invalid hash function: {}, only 'Murmur64' with 'protein', 'dayhoff' or 'hp' \
             are supported",
            hash_function
        )),
    }
}

pub fn get_moltype_from_hash_function(
    hash_function: HashFunctions,
) -> Result<String, anyhow::Error> {
    match hash_function {
        HashFunctions::Murmur64Protein => Ok("protein20".to_string()),
        HashFunctions::Murmur64Dayhoff => Ok("dayhoff6".to_string()),
        HashFunctions::Murmur64Hp => Ok("hp_lehninger2".to_string()),
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
///   - `"protein20"` for the full 20-letter alphabet
///   - `"dayhoff6"` for the Dayhoff alphabet
///   - `"hp_<name><size>"` for an HP table, e.g. `"hp_thomas_dill2"`
///   - `"<name><size>"` for a multi-letter reduced alphabet, e.g. `"sdm12"`
///
/// # Returns
///
/// * `Ok(fn(u8) -> u8)` - A function that encodes an amino acid byte according to
///    the specified `moltype`.
/// * `Err(...)` - An error if the `moltype` is unrecognized.
pub fn get_encoding_fn_from_moltype(moltype: &str) -> Result<fn(u8) -> u8, anyhow::Error> {
    match canonical_moltype(moltype) {
        "protein20" => Ok(|b| b),
        // Before the table lookup below, for the same reason as in
        // get_hash_function_from_moltype: dayhoff is encoded by sourmash, not pre-encoded.
        "dayhoff6" => Ok(aa_to_dayhoff),
        // Encoded by sourmash, like dayhoff6 above.
        "hp_lehninger2" => Ok(aa_to_hp),
        // Table-backed alphabets cannot be expressed as a fn(u8) -> u8; return identity here
        // so callers that only need one don't crash. encode_by_moltype and add_protein apply
        // the table themselves.
        s if alphabet_table(s).is_some() => Ok(|b| b),
        _ => Err(anyhow::anyhow!(
            "Invalid alphabet: {}, only 'protein20', 'dayhoff6' and the reduced \
             alphabets are supported",
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
/// * `moltype` - The alphabet to encode with, e.g. "protein20", "dayhoff6", or any
///   reduced alphabet name
///
/// # Returns
/// * `Ok(String)` - The encoded sequence
/// * `Err(...)` - An error if the moltype is invalid
///
/// # Example
/// ```
/// use kmerseek::hash_functions::encode_by_moltype;
///
/// let sequence = "MKTAYIAKQR";
/// let encoded = encode_by_moltype(sequence, "hp_lehninger2").unwrap();
/// // encoded will be the HP-encoded version
/// ```
pub fn encode_by_moltype(sequence: &str, moltype: &str) -> Result<String> {
    // Table-backed alphabets cannot be expressed as a fn(u8) -> u8, so they are applied
    // here directly. Without this the HP family would encode to the identity.
    if let Some(table) = alphabet_table(moltype) {
        return Ok(sequence
            .bytes()
            .map(|b| {
                let upper = b.to_ascii_uppercase();
                table.get(&upper).copied().unwrap_or(upper) as char
            })
            .collect());
    }
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
/// use kmerseek::hash_functions::encode_with_fn;
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
        if let Ok(hf) = get_hash_function_from_moltype("protein20") {
            assert_eq!(hf, HashFunctions::Murmur64Protein);
        } else {
            panic!("Expected HashFunctions::Murmur64Protein for 'protein20'");
        }

        if let Ok(hf) = get_hash_function_from_moltype("protein20") {
            assert_eq!(hf, HashFunctions::Murmur64Protein);
        } else {
            panic!("Expected HashFunctions::Murmur64Protein for 'raw'");
        }

        // Lehninger is sourmash's own partition, so it keeps sourmash's HP hash function and
        // its sketches match sourmash's byte for byte.
        if let Ok(hf) = get_hash_function_from_moltype("hp_lehninger2") {
            assert_eq!(hf, HashFunctions::Murmur64Hp);
        } else {
            panic!("Expected HashFunctions::Murmur64Hp for 'hp_lehninger2'");
        }

        // Alphabets with no sourmash equivalent are pre-encoded, so they hash as protein.
        if let Ok(hf) = get_hash_function_from_moltype("sdm12") {
            assert_eq!(hf, HashFunctions::Murmur64Protein);
        } else {
            panic!("Expected HashFunctions::Murmur64Protein for 'sdm12'");
        }

        if let Ok(hf) = get_hash_function_from_moltype("dayhoff6") {
            assert_eq!(hf, HashFunctions::Murmur64Dayhoff);
        } else {
            panic!("Expected HashFunctions::Murmur64Dayhoff for 'dayhoff6'");
        }

        assert!(get_hash_function_from_moltype("xyz").is_err());
    }

    #[test]
    fn test_get_encoding_fn_from_moltype() {
        if let Ok(raw_fn) = get_encoding_fn_from_moltype("protein20") {
            assert_eq!(raw_fn(b'A'), b'A');
        } else {
            panic!("Expected identity function for 'raw'");
        }

        if let Ok(protein_fn) = get_encoding_fn_from_moltype("protein20") {
            assert_eq!(protein_fn(b'F'), b'F');
        } else {
            panic!("Expected identity function for 'protein20'");
        }

        // Lehninger is sourmash's own partition, so it encodes through sourmash's aa_to_hp
        // rather than a table of ours.
        if let Ok(hp_fn) = get_encoding_fn_from_moltype("hp_lehninger2") {
            assert_eq!(hp_fn(b'A'), aa_to_hp(b'A'));
        } else {
            panic!("Expected aa_to_hp for 'hp_lehninger2'");
        }

        // A table-backed alphabet cannot be a fn(u8) -> u8, so it falls back to the
        // identity here; encode_by_moltype applies the table instead.
        if let Ok(td_fn) = get_encoding_fn_from_moltype("hp_thomas_dill2") {
            assert_eq!(td_fn(b'A'), b'A');
        } else {
            panic!("Expected identity function for 'hp_thomas_dill2'");
        }

        if let Ok(dayhoff_fn) = get_encoding_fn_from_moltype("dayhoff6") {
            assert_eq!(dayhoff_fn(b'G'), aa_to_dayhoff(b'G'));
        } else {
            panic!("Expected aa_to_dayhoff for 'dayhoff6'");
        }

        assert!(get_encoding_fn_from_moltype("xyz").is_err());
    }

    #[test]
    fn test_encode_by_moltype_protein() -> Result<()> {
        let encoded = encode_by_moltype(TEST_KMER, "protein20")?;
        assert_eq!(encoded, TEST_KMER);
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_dayhoff() -> Result<()> {
        let encoded = encode_by_moltype(TEST_KMER, "dayhoff6")?;
        assert_eq!(encoded, "eeeecbbeeec");
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_hp() -> Result<()> {
        let encoded = encode_by_moltype(TEST_KMER, "hp_lehninger2")?;
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
        let encoded = encode_by_moltype(sequence, "protein20")?;
        assert_eq!(encoded, sequence);
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_sequence_hp() -> Result<()> {
        let sequence = "MKTAYIAKQR";
        let encoded = encode_by_moltype(sequence, "hp_lehninger2")?;
        // Verify the encoding produces the correct length and uses only 'h' and 'p'
        assert_eq!(encoded.len(), sequence.len());
        assert!(encoded.chars().all(|c| c == 'h' || c == 'p'));
        Ok(())
    }

    #[test]
    fn test_encode_by_moltype_sequence_dayhoff() -> Result<()> {
        let sequence = "MKTAYIAKQR";
        let encoded = encode_by_moltype(sequence, "dayhoff6")?;
        // Dayhoff encoding should produce a different string
        assert_ne!(encoded, sequence);
        assert_eq!(encoded.len(), sequence.len());
        Ok(())
    }

    /// Each hash function sourmash writes maps back to the kmerseek name for the same
    /// alphabet, which is what lets an index sourmash produced open under our naming.
    #[test]
    fn test_get_moltype_from_hash_function_each_supported_variant() -> Result<()> {
        assert_eq!(get_moltype_from_hash_function(HashFunctions::Murmur64Protein)?, "protein20");
        assert_eq!(get_moltype_from_hash_function(HashFunctions::Murmur64Dayhoff)?, "dayhoff6");
        assert_eq!(get_moltype_from_hash_function(HashFunctions::Murmur64Hp)?, "hp_lehninger2");
        Ok(())
    }

    /// DNA is a hash function sourmash writes that kmerseek has no alphabet for, so it
    /// has to be refused rather than read as some protein alphabet.
    #[test]
    fn test_get_moltype_from_hash_function_rejects_dna() {
        let err = get_moltype_from_hash_function(HashFunctions::Murmur64Dna)
            .expect_err("DNA has no protein alphabet to map to");
        assert_eq!(
            err.to_string(),
            "Invalid hash function: DNA, only Sourmash HashFunctions::Murmur64 with \
             'protein', 'dayhoff', or 'hp' are supported"
        );
    }

    /// Table-backed alphabets are applied by `encode_by_moltype` itself, since a table
    /// cannot be expressed as a `fn(u8) -> u8`. Without that branch they would fall
    /// through to the identity and encode to the input unchanged.
    #[test]
    fn test_encode_by_moltype_applies_the_reduced_table() -> Result<()> {
        // SDM12 writes each class as its first residue in lowercase, so KER is k,
        // TSQ is t, YF is y and LIVM is l.
        assert_eq!(encode_by_moltype("MKTAYIAKQR", "sdm12")?, "lktaylaktk");
        // The HP family reaches the same branch: h is ACFILMVWY, p is DEGHKNPQRST.
        assert_eq!(encode_by_moltype("MKTAYIAKQR", "hp_thomas_dill2")?, "hpphhhhppp");
        Ok(())
    }

    /// The table is keyed on uppercase residues, and anything it does not hold passes
    /// through unchanged so the encoded sequence keeps the length of its input.
    #[test]
    fn test_encode_by_moltype_uppercases_input_and_keeps_unknown_residues() -> Result<()> {
        assert_eq!(encode_by_moltype("mktayiakqr", "sdm12")?, "lktaylaktk");
        // X (any residue) and Z (Glx) are not in the 20-residue table.
        assert_eq!(encode_by_moltype("MXKZ", "sdm12")?, "lXkZ");
        Ok(())
    }
}
