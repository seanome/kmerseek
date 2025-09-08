use serde::{Deserialize, Serialize};
use std::fmt;

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
    /// Create a new molecular type with validation
    pub fn new(moltype: &str) -> Result<Self, String> {
        match moltype {
            "protein" | "dayhoff" | "hp" => Ok(MolType(moltype.to_string())),
            _ => Err(format!(
                "Invalid molecular type: {}. Must be one of: protein, dayhoff, hp",
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
}
