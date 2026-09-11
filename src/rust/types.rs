use serde::{Deserialize, Serialize};
use std::fmt;

use crate::alphabets::{canonical_moltype, Alphabet};

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
    /// A k-mer is kept when its hash falls in the lowest `1/scaled` of the hash space
    /// (FracMinHash), so the same k-mer is kept or dropped in every sequence and the
    /// expected fraction kept is `1/scaled`. Which positions survive is decided by
    /// hash value, not by stride: scaled=2 keeps about half the k-mers, not every
    /// second one.
    ///
    /// Capped at 10. Beyond that a typical protein keeps too few k-mers for a matched
    /// region to be sampled at all.
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
    /// Create a molecular type, rejecting anything that is not a known alphabet.
    ///
    /// sourmash's `protein`, `dayhoff` and `hp` are accepted and stored under kmerseek's
    /// names for the same alphabets. Normalizing here rather than at each call site keeps a
    /// single spelling in indexes and results: `find_matched_regions` asserts query and
    /// target moltypes are equal, and two names for one alphabet would abort the search.
    pub fn new(moltype: &str) -> Result<Self, String> {
        if let Some(alphabet) = Alphabet::from_moltype(canonical_moltype(moltype)) {
            return Ok(MolType(alphabet.to_moltype().to_string()));
        }
        Err(format!(
            "Invalid molecular type: {}. Must be one of: protein20, dayhoff6, \
             hp_lehninger2, hp_thomas_dill2, hp_kyte_doolittle2, hp_thomas_dill_no_c2, \
             hp_lehninger_c_nonpolar2, hp_lehninger_hpc3, hp_pbotc_1st_ed2, \
             gbmr4, polarity4, wwmj5, gbmr7, funcgroups8, sdm12, \
             mmseqs12, wass14, hsdm17, uniprot18 (sourmash's protein, dayhoff and hp are \
             also read)",
            moltype
        ))
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
        assert!(MolType::new("protein20").is_ok());
        assert!(MolType::new("dayhoff6").is_ok());
        assert!(MolType::new("hp_lehninger2").is_ok());
        assert!(MolType::new("invalid").is_err());
    }

    /// `protein` and its synonym `raw` are the full 20-letter alphabet and keep their names,
    /// as do the alphabets that already carry a class count.
    #[test]
    fn test_moltype_leaves_current_names_alone() {
        for moltype in ["protein20", "dayhoff6", "hp_lehninger2", "sdm12", "gbmr4", "uniprot18"] {
            assert_eq!(MolType::new(moltype).unwrap().get(), moltype);
        }
    }

    /// kmerseek reads sourmash's three moltypes and stores them under its own names, so
    /// sourmash-labelled data stays usable and only one spelling reaches the rest of the
    /// code.
    #[test]
    fn test_moltype_reads_sourmash_names() {
        for (sourmash, kmerseek) in
            [("protein", "protein20"), ("dayhoff", "dayhoff6"), ("hp", "hp_lehninger2")]
        {
            assert_eq!(MolType::new(sourmash).unwrap().get(), kmerseek, "{sourmash}");
            assert_eq!(MolType::new(sourmash).unwrap(), MolType::new(kmerseek).unwrap());
        }
    }

    /// kmerseek's own earlier spellings were never sourmash's, so they stay rejected.
    #[test]
    fn test_moltype_rejects_earlier_kmerseek_spellings() {
        for moltype in ["raw", "hp_lehninger", "hp_thomas_dill", "reduced_sdm12"] {
            assert!(MolType::new(moltype).is_err(), "{moltype}");
        }
    }
}
