use anyhow::{bail, Result};
use rand::prelude::*;
use std::collections::HashMap;

/// Standard amino acids and their properties
pub const STANDARD_AA: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y',
];

/// Represents amino acid ambiguity codes and their possible resolutions
#[derive(Debug)]
pub struct AminoAcidAmbiguity {
    /// Maps ambiguous amino acid codes to their possible resolutions
    replacements: HashMap<char, Vec<char>>,
}

impl AminoAcidAmbiguity {
    pub fn new() -> Self {
        let mut replacements = HashMap::new();

        // Standard amino acids map to themselves
        for aa in STANDARD_AA.iter() {
            replacements.insert(*aa, vec![*aa]);
        }

        // Add ambiguous amino acid codes
        replacements.insert('X', STANDARD_AA.to_vec()); // Unknown - any amino acid
        replacements.insert('B', vec!['D', 'N']); // Aspartic acid or Asparagine
        replacements.insert('Z', vec!['E', 'Q']); // Glutamic acid or Glutamine
        replacements.insert('J', vec!['I', 'L']); // Isoleucine or Leucine

        AminoAcidAmbiguity { replacements }
    }

    pub fn is_valid_aa(&self, aa: char) -> bool {
        self.replacements.contains_key(&aa)
    }

    pub fn resolve_ambiguity(&self, aa: char) -> char {
        if let Some(possible_aas) = self.replacements.get(&aa) {
            // If there's only one possibility, use it
            if possible_aas.len() == 1 {
                possible_aas[0]
            } else {
                // Otherwise randomly choose one of the possibilities
                let mut rng = rand::rng();
                *possible_aas.choose(&mut rng).unwrap_or(&aa)
            }
        } else {
            // If not found in replacements, return the original character
            aa
        }
    }

    /// Validates a protein sequence and returns an error if invalid characters are found
    pub fn validate_sequence(&self, sequence: &str) -> Result<()> {
        for (i, c) in sequence.chars().enumerate() {
            if !self.is_valid_aa(c) {
                bail!("Invalid amino acid '{}' found at position {}", c, i + 1);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_amino_acids() {
        let aa = AminoAcidAmbiguity::new();

        // Test standard amino acids
        for c in STANDARD_AA.iter() {
            assert!(aa.is_valid_aa(*c));
        }

        // Test ambiguous codes
        assert!(aa.is_valid_aa('X'));
        assert!(aa.is_valid_aa('B'));
        assert!(aa.is_valid_aa('Z'));
        assert!(aa.is_valid_aa('J'));

        // Test invalid characters
        assert!(!aa.is_valid_aa('1'));
        assert!(!aa.is_valid_aa('$'));
    }

    #[test]
    fn test_ambiguity_resolution() {
        let aa = AminoAcidAmbiguity::new();

        // Test standard amino acids (should return themselves)
        for c in STANDARD_AA.iter() {
            assert_eq!(aa.resolve_ambiguity(*c), *c);
        }

        // Test ambiguous codes (should return one of their possible resolutions)
        let b_resolution = aa.resolve_ambiguity('B');
        assert!(vec!['D', 'N'].contains(&b_resolution));

        let z_resolution = aa.resolve_ambiguity('Z');
        assert!(vec!['E', 'Q'].contains(&z_resolution));

        let j_resolution = aa.resolve_ambiguity('J');
        assert!(vec!['I', 'L'].contains(&j_resolution));
    }

    #[test]
    fn test_sequence_validation() {
        let aa = AminoAcidAmbiguity::new();

        // Test valid sequence
        assert!(aa.validate_sequence("ACDEFGHIKLMNPQRSTVWY").is_ok());
        assert!(aa.validate_sequence("ACDEFXBZJ").is_ok());

        // Test invalid sequence
        let result = aa.validate_sequence("ACDEF1GHIKL");
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Invalid amino acid '1'"));
    }
}
