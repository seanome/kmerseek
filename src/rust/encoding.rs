use anyhow::Result;
use sourmash::encodings::{aa_to_dayhoff, HashFunctions};
use std::collections::HashMap;

use crate::hp_alphabets::HpAlphabet;
use crate::reduced_alphabets::ReducedAlphabet;

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
    match moltype {
        "protein20" | "protein" | "raw" => Ok(HashFunctions::Murmur64Protein),
        // Dayhoff has no table of ours; sourmash encodes and hashes it. This arm has to come
        // before the table lookup below, which assumes a pre-encoded sequence.
        "dayhoff6" | "dayhoff" => Ok(HashFunctions::Murmur64Dayhoff),
        // Every other alphabet pre-encodes through its own table, so the hash function sees
        // an already-encoded sequence and hashes it as protein. This includes bare `hp`,
        // which is a spelling of hp_lehninger2: it must NOT map to Murmur64Hp, or sourmash
        // would run its own HP encoder over an already-encoded string, encoding twice.
        s if custom_alphabet_table(s).is_some() => Ok(HashFunctions::Murmur64Protein),
        _ => Err(anyhow::anyhow!(
            "Invalid alphabet: {}. Supported values: 'protein20', 'dayhoff6', \
             'hp_lehninger2', 'hp_thomas_dill2', 'hp_kyte_doolittle2', \
             'hp_thomas_dill_no_c2', 'hp_lehninger_c_nonpolar2', 'hp_lehninger_hpc3', \
             'hp_pbotc_1st_ed2', 'hp_random_control2', \
             'hp_random_control2_1'..'hp_random_control2_10', \
             'gbmr4', 'wwmj5', 'gbmr7', 'sdm12', 'mmseqs12', 'wass14', 'hsdm17', \
             'uniprot18' (older spellings -- 'protein', 'raw', 'hp', 'dayhoff', \
             'hp_<name>' without a class count, a 'reduced_' prefix, 'shuffled_control' \
             -- are still accepted)",
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
        // Murmur64Hp is sourmash's own HP encoding, which hashes lowercase h/p. kmerseek
        // hashes the same partition as uppercase H/P via hp_lehninger2, so a
        // signature written this way is not interchangeable with ours and is refused rather
        // than relabelled.
        MURMUR64HP => Err(anyhow::anyhow!(
            "Signature uses sourmash's built-in HP hash function (Murmur64Hp), whose hashes \
             are not compatible with kmerseek's hp_lehninger2. Re-sketch the input \
             with kmerseek."
        )),
        _ => Err(anyhow::anyhow!(
            "Invalid hash function: {}, only 'Murmur64' with 'protein' or 'dayhoff' are supported",
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
        // See get_moltype_from_hash_function_string: sourmash's HP hashes differ from ours.
        HashFunctions::Murmur64Hp => Err(anyhow::anyhow!(
            "Signature uses sourmash's built-in HP hash function (Murmur64Hp), whose hashes \
             are not compatible with kmerseek's hp_lehninger2. Re-sketch the input \
             with kmerseek."
        )),
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
    match moltype {
        "protein20" | "protein" | "raw" => Ok(|b| b),
        // Before the table lookup below, for the same reason as in
        // get_hash_function_from_moltype: dayhoff is encoded by sourmash, not pre-encoded.
        "dayhoff6" | "dayhoff" => Ok(aa_to_dayhoff),
        // Table-backed alphabets cannot be expressed as a fn(u8) -> u8; return identity here
        // so callers that only need one don't crash. encode_by_moltype and add_protein apply
        // the table themselves.
        s if custom_alphabet_table(s).is_some() => Ok(|b| b),
        _ => Err(anyhow::anyhow!(
            "Invalid alphabet: {}, only 'protein20', 'dayhoff6' and the reduced \
             alphabets are supported",
            moltype
        )),
    }
}

/// Residue-to-symbol table for the moltypes that pre-encode a sequence before hashing.
///
/// Covers both alphabet families: the two- and three-letter HP tables (`hp_*2`/`hp_*3`,
/// plus bare `hp`) and the multi-letter reduced alphabets (`sdm12`, `gbmr4`, ...). Returns
/// `None` for `protein20` and `dayhoff6`, which sourmash encodes itself.
pub fn custom_alphabet_table(moltype: &str) -> Option<&'static HashMap<u8, u8>> {
    if let Some(alphabet) = HpAlphabet::from_moltype(moltype) {
        return Some(alphabet.table());
    }
    ReducedAlphabet::from_moltype(moltype).map(|alphabet| alphabet.table())
}

/// Encode a sequence using the specified moltype.
///
/// This function applies the encoding function to each amino acid in the sequence,
/// producing an encoded sequence string. The encoded sequence is pre-allocated
/// with the exact capacity needed.
///
/// # Arguments
/// * `sequence` - The sequence to encode (can be a k-mer or full sequence)
/// * `moltype` - The molecule type encoding to use, e.g. "protein", "reduced_dayhoff6",
///   or any reduced alphabet name
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
/// let encoded = encode_by_moltype(sequence, "hp_lehninger2").unwrap();
/// // encoded will be the HP-encoded version
/// ```
pub fn encode_by_moltype(sequence: &str, moltype: &str) -> Result<String> {
    // Table-backed alphabets cannot be expressed as a fn(u8) -> u8, so they are applied
    // here directly. Without this the HP family would silently encode to the identity.
    if let Some(table) = custom_alphabet_table(moltype) {
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
        if let Ok(hf) = get_hash_function_from_moltype("protein20") {
            assert_eq!(hf, HashFunctions::Murmur64Protein);
        } else {
            panic!("Expected HashFunctions::Murmur64Protein for 'protein20'");
        }

        if let Ok(hf) = get_hash_function_from_moltype("raw") {
            assert_eq!(hf, HashFunctions::Murmur64Protein);
        } else {
            panic!("Expected HashFunctions::Murmur64Protein for 'raw'");
        }

        // `hp` is the pre-rename spelling of reduced_hp_lehninger2, which pre-encodes with
        // our own table and so hashes as protein. It is deliberately no longer Murmur64Hp.
        if let Ok(hf) = get_hash_function_from_moltype("hp") {
            assert_eq!(hf, HashFunctions::Murmur64Protein);
        } else {
            panic!("Expected HashFunctions::Murmur64Protein for 'hp'");
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
        if let Ok(raw_fn) = get_encoding_fn_from_moltype("raw") {
            assert_eq!(raw_fn(b'A'), b'A');
        } else {
            panic!("Expected identity function for 'raw'");
        }

        if let Ok(protein_fn) = get_encoding_fn_from_moltype("protein20") {
            assert_eq!(protein_fn(b'F'), b'F');
        } else {
            panic!("Expected identity function for 'protein20'");
        }

        // `hp` now names the table-backed Lehninger alphabet, and a table cannot be a
        // fn(u8) -> u8, so this returns the identity. Callers that want the encoding go
        // through encode_by_moltype, which consults the table
        // (test_encode_by_moltype_hp_uses_the_lehninger_table).
        if let Ok(hp_fn) = get_encoding_fn_from_moltype("hp") {
            assert_eq!(hp_fn(b'A'), b'A');
        } else {
            panic!("Expected identity function for 'hp'");
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
        let encoded = encode_by_moltype(TEST_KMER, "protein")?;
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
        let encoded = encode_by_moltype(sequence, "dayhoff6")?;
        // Dayhoff encoding should produce a different string
        assert_ne!(encoded, sequence);
        assert_eq!(encoded.len(), sequence.len());
        Ok(())
    }
}
