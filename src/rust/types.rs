use serde::{Deserialize, Serialize};
use std::fmt;

use crate::hp_alphabets::HpAlphabet;

/// A type-safe wrapper for k-mer sizes
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct KmerSize(pub u32);

impl KmerSize {
    /// Create a new k-mer size with validation
    pub fn new(size: u32) -> Result<Self, String> {
        if size == 0 {
            Err("K-mer size must be greater than 0".to_string())
        } else if size > 100 {
            Err("K-mer size too large (max 100)".to_string())
        } else {
            Ok(KmerSize(size))
        }
    }

    /// Get the raw value
    pub fn get(&self) -> u32 {
        self.0
    }

    /// Get as usize for indexing
    pub fn as_usize(&self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for KmerSize {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A type-safe wrapper for scaled values
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Scaled(pub u32);

impl Scaled {
    /// Create a new scaled value with validation
    ///
    /// Scaled represents the sampling rate (1/scaled), so:
    /// - scaled=1: take every k-mer (100% sampling)
    /// - scaled=2: take every 2nd k-mer (50% sampling)  
    /// - scaled=5: take every 5th k-mer (20% sampling)
    /// - scaled=10: take every 10th k-mer (10% sampling)
    ///
    /// For protein analysis, we typically want scaled ≤ 10 to ensure
    /// meaningful k-mer coverage. Higher values result in too sparse sampling.
    pub fn new(scaled: u32) -> Result<Self, String> {
        if scaled == 0 {
            Err("Scaled value must be greater than 0".to_string())
        } else if scaled > 10 {
            Err(format!(
                "Scaled value too large: {}. For protein analysis, scaled should be ≤ 10 to ensure meaningful k-mer coverage. \
                Higher values result in too sparse sampling (e.g., scaled=100 means only ~1 k-mer per 100-amino acid protein).",
                scaled
            ))
        } else {
            Ok(Scaled(scaled))
        }
    }

    /// Get the raw value
    pub fn get(&self) -> u32 {
        self.0
    }
}

impl fmt::Display for Scaled {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A type-safe wrapper for molecular types
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct MolType(pub String);

impl MolType {
    /// Create a new molecular type with validation.
    ///
    /// Pre-rename HP spellings (`hp_thomas_dill`) are accepted and normalized to their
    /// current form (`reduced_hp_thomas_dill2`). Normalizing here rather than at each call
    /// site means an index written before the rename compares equal to a query sketched
    /// after it: `find_matched_regions` asserts the two moltypes match, and the same
    /// alphabet under two spellings would otherwise abort the search.
    pub fn new(moltype: &str) -> Result<Self, String> {
        if let Some(alphabet) = HpAlphabet::from_moltype(moltype) {
            return Ok(MolType(alphabet.to_moltype()));
        }
        match moltype {
            // Dayhoff carries its class count too; the bare name is the pre-rename
            // spelling and normalizes to it.
            "dayhoff" => Ok(MolType("reduced_dayhoff6".to_string())),
            // `raw` is a synonym for `protein` in encoding.rs; both keep their own name so
            // generated filenames stay stable.
            "protein" | "raw" => Ok(MolType(moltype.to_string())),
            s if s.starts_with("reduced_") => Ok(MolType(moltype.to_string())),
            _ => Err(format!(
                "Invalid molecular type: {}. Must be one of: protein, reduced_dayhoff6, hp, \
                 reduced_hp_lehninger2, reduced_hp_thomas_dill2, reduced_hp_kyte_doolittle2, \
                 reduced_hp_thomas_dill_no_c2, reduced_hp_lehninger_c_nonpolar2, \
                 reduced_hp_lehninger_hpc3, reduced_hp_pbotc_1st_ed2, \
                 reduced_hp_shuffled_control2, reduced_gbmr4, reduced_wwmj5, reduced_gbmr7, \
                 reduced_sdm12, reduced_mmseqs12, reduced_wass14, reduced_hsdm17, \
                 reduced_uniprot18 (the pre-rename hp_<name> and dayhoff spellings are \
                 still accepted)",
                moltype
            )),
        }
    }

    /// Get the raw value
    pub fn get(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for MolType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A type-safe wrapper for hash values
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct HashValue(pub u64);

impl HashValue {
    /// Create a new hash value
    pub fn new(hash: u64) -> Self {
        HashValue(hash)
    }

    /// Get the raw value
    pub fn get(&self) -> u64 {
        self.0
    }
}

impl fmt::Display for HashValue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A type-safe wrapper for sequence positions
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Position(pub usize);

impl Position {
    /// Create a new position
    pub fn new(pos: usize) -> Self {
        Position(pos)
    }

    /// Get the raw value
    pub fn get(&self) -> usize {
        self.0
    }
}

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_size_validation() {
        assert!(KmerSize::new(0).is_err());
        assert!(KmerSize::new(1).is_ok());
        assert!(KmerSize::new(5).is_ok());
        assert!(KmerSize::new(100).is_ok());
        assert!(KmerSize::new(101).is_err());
    }

    #[test]
    fn test_scaled_validation() {
        assert!(Scaled::new(0).is_err());
        assert!(Scaled::new(1).is_ok());
        assert!(Scaled::new(2).is_ok());
        assert!(Scaled::new(5).is_ok());
        assert!(Scaled::new(10).is_ok());
        assert!(Scaled::new(11).is_err());
        assert!(Scaled::new(100).is_err());
        assert!(Scaled::new(1000).is_err());

        // Test error message for large values
        let error = Scaled::new(100).unwrap_err();
        assert!(error.contains("too large"));
        assert!(error.contains("≤ 10"));
        assert!(error.contains("sparse sampling"));
    }

    #[test]
    fn test_moltype_validation() {
        assert!(MolType::new("protein").is_ok());
        assert!(MolType::new("dayhoff").is_ok());
        assert!(MolType::new("hp").is_ok());
        assert!(MolType::new("invalid").is_err());
    }

    /// A moltype written before class counts were added to the HP names must normalize to
    /// the current spelling. `find_matched_regions` asserts query and target moltypes are
    /// equal, so an index built as `hp_thomas_dill` and a query sketched afterwards would
    /// otherwise abort the search even though both use the same alphabet.
    #[test]
    fn test_moltype_normalizes_legacy_hp_names() {
        let expected = [
            ("hp_lehninger", "reduced_hp_lehninger2"),
            ("hp_thomas_dill", "reduced_hp_thomas_dill2"),
            ("hp_kyte_doolittle", "reduced_hp_kyte_doolittle2"),
            ("hp_thomas_dill_no_c", "reduced_hp_thomas_dill_no_c2"),
            ("hp_lehninger_c_nonpolar", "reduced_hp_lehninger_c_nonpolar2"),
            ("hp_lehninger_hpc", "reduced_hp_lehninger_hpc3"),
            ("hp_pbotc_1st_ed", "reduced_hp_pbotc_1st_ed2"),
            ("hp_shuffled_control", "reduced_hp_shuffled_control2"),
            ("hp_shuffled_control_4", "reduced_hp_shuffled_control2_4"),
        ];

        for (legacy, current) in expected {
            assert_eq!(MolType::new(legacy).unwrap().get(), current, "{legacy}");
            // Already-current names pass through unchanged.
            assert_eq!(MolType::new(current).unwrap().get(), current, "{current}");
            assert_eq!(MolType::new(legacy).unwrap(), MolType::new(current).unwrap());
        }
    }

    /// `dayhoff` carries its class count now, so the bare name normalizes the way the
    /// pre-rename HP names do. Its hash function is unchanged, so existing dayhoff indexes
    /// keep matching (see `test_dayhoff_rename_preserves_hashes`).
    #[test]
    fn test_moltype_normalizes_legacy_dayhoff() {
        assert_eq!(MolType::new("dayhoff").unwrap().get(), "reduced_dayhoff6");
        assert_eq!(MolType::new("reduced_dayhoff6").unwrap().get(), "reduced_dayhoff6");
        assert_eq!(MolType::new("dayhoff").unwrap(), MolType::new("reduced_dayhoff6").unwrap());
    }

    /// `hp` is now the pre-rename spelling of the Lehninger 2-class alphabet.
    #[test]
    fn test_moltype_normalizes_builtin_hp_to_lehninger2() {
        assert_eq!(MolType::new("hp").unwrap().get(), "reduced_hp_lehninger2");
        assert_eq!(MolType::new("hp").unwrap(), MolType::new("reduced_hp_lehninger2").unwrap());
    }

    /// `protein` and its synonym `raw` are the full 20-letter alphabet and keep their names,
    /// as do the alphabets that already carry a class count.
    #[test]
    fn test_moltype_leaves_current_names_alone() {
        for moltype in ["protein", "raw", "reduced_sdm12", "reduced_gbmr4", "reduced_uniprot18"] {
            assert_eq!(MolType::new(moltype).unwrap().get(), moltype);
        }
    }
}
