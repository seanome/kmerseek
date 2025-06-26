use anyhow::{bail, Result};
use rand::prelude::*;
use std::borrow::Cow;
use std::collections::HashMap;

/// Standard amino acids and their properties
pub const STANDARD_AA: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y',
];

pub const SPECIAL_AA: [char; 4] = ['X', 'U', 'O', '*'];

/// Represents amino acid ambiguity codes and their possible resolutions
#[derive(Debug)]
pub struct AminoAcidAmbiguity {
    /// Maps ambiguous amino acid codes to their possible resolutions
    replacements: HashMap<char, Vec<char>>,
}

impl AminoAcidAmbiguity {
    pub fn new() -> Self {
        let mut replacements = HashMap::new();

        // Add ambiguous amino acid codes
        replacements.insert('B', vec!['D', 'N']); // Aspartic acid or Asparagine
        replacements.insert('Z', vec!['E', 'Q']); // Glutamic acid or Glutamine
        replacements.insert('J', vec!['I', 'L']); // Isoleucine or Leucine

        AminoAcidAmbiguity { replacements }
    }

    fn is_valid_aa(&self, aa: char) -> bool {
        // Put the most common (99.9%) case, of the 20 standard amino acids first
        STANDARD_AA.contains(&aa) || SPECIAL_AA.contains(&aa) || self.replacements.contains_key(&aa)
    }

    fn resolve_ambiguity(&self, aa: char) -> char {
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
    /// Stops reading at the first stop codon (*)
    pub fn validate_sequence(&self, sequence: &str) -> Result<()> {
        for (i, c) in sequence.chars().enumerate() {
            if c == '*' {
                // Stop codon - this is valid, but we stop reading here
                return Ok(());
            }
            if !self.is_valid_aa(c) {
                bail!("Invalid amino acid '{}' found at position {}", c, i + 1);
            }
        }
        Ok(())
    }

    /// Validates a protein sequence and resolves ambiguity if needed.
    /// Returns Ok(Cow<str>) with ambiguity resolved if needed, or an error if invalid characters are found.
    /// Stops processing at the first stop codon (*).
    pub fn validate_and_resolve<'a>(&self, sequence: &'a str) -> Result<Cow<'a, str>> {
        // First validate the sequence (this will stop at stop codon)
        self.validate_sequence(sequence)?;

        // Find the position of the first stop codon, if any
        let stop_pos = sequence.find('*');
        let sequence_to_process = if let Some(pos) = stop_pos {
            &sequence[..pos + 1] // Include the stop codon
        } else {
            sequence
        };

        let has_ambiguous = sequence_to_process.chars().any(|c| matches!(c, 'B' | 'Z' | 'J'));
        if has_ambiguous {
            // Only replace ambiguous characters if needed, since <0.01% of amino acids in UniProtKB are ambiguous
            // See: https://www.uniprot.org/uniprotkb/statistics#amino-acid-composition
            let resolved: String =
                sequence_to_process.chars().map(|c| self.resolve_ambiguity(c)).collect();
            Ok(Cow::Owned(resolved))
        } else {
            Ok(Cow::Borrowed(sequence_to_process))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::test_fixtures::{TEST_PROTEIN, TEST_PROTEIN_INVALID};

    #[test]
    fn test_valid_amino_acids() {
        let aa = AminoAcidAmbiguity::new();

        // Test standard amino acids
        for c in STANDARD_AA.iter() {
            assert!(aa.is_valid_aa(*c));
        }

        // Test ambiguous codes
        assert!(aa.is_valid_aa('B'));
        assert!(aa.is_valid_aa('Z'));
        assert!(aa.is_valid_aa('J'));

        // Test special amino acids
        assert!(aa.is_valid_aa('X'));
        assert!(aa.is_valid_aa('U'));
        assert!(aa.is_valid_aa('O'));
        assert!(aa.is_valid_aa('*'));

        // Test invalid characters
        assert!(!aa.is_valid_aa('1'));
        assert!(!aa.is_valid_aa('$'));
        assert!(!aa.is_valid_aa('@'));
    }

    #[test]
    fn test_ambiguity_resolution() {
        let aa = AminoAcidAmbiguity::new();

        // Test standard amino acids (should return themselves)
        for c in STANDARD_AA.iter() {
            assert_eq!(aa.resolve_ambiguity(*c), *c);
        }

        // Test special amino acids (should return themselves)
        assert_eq!(aa.resolve_ambiguity('X'), 'X');
        assert_eq!(aa.resolve_ambiguity('U'), 'U');
        assert_eq!(aa.resolve_ambiguity('O'), 'O');
        assert_eq!(aa.resolve_ambiguity('*'), '*');

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
        assert!(aa.validate_sequence(TEST_PROTEIN).is_ok());
        assert!(aa.validate_sequence("ACDEFGHIKLMNPQRSTVWY").is_ok());
        assert!(aa.validate_sequence("ACDEFXBZJ").is_ok());

        // Test sequences with special amino acids
        assert!(aa.validate_sequence("ACDEFXUO").is_ok());
        assert!(aa.validate_sequence("ACDEF*").is_ok());
        assert!(aa.validate_sequence("ACDEF*GHI").is_ok()); // Should stop at *

        // Test invalid sequence
        let result = aa.validate_sequence(TEST_PROTEIN_INVALID);
        assert!(result.is_err());
        // Contains invalid character '1'
        assert!(result.unwrap_err().to_string().contains("Invalid amino acid '1'"));
    }

    #[test]
    fn test_validate_and_resolve_with_stop_codon() {
        let aa = AminoAcidAmbiguity::new();

        // Test sequence with stop codon
        let result = aa.validate_and_resolve("ACDEF*GHI");
        assert!(result.is_ok());
        assert_eq!(result.unwrap().as_ref(), "ACDEF*");

        // Test sequence with stop codon and ambiguous amino acids
        let result = aa.validate_and_resolve("ACDEFB*GHI");
        assert!(result.is_ok());
        let resolved = result.unwrap();
        assert!(resolved.as_ref().starts_with("ACDEF"));
        assert!(resolved.as_ref().ends_with("*"));
        // The B should be resolved to either D or N
        let fifth_char = resolved.as_ref().chars().nth(5).unwrap();
        assert!(vec!['D', 'N'].contains(&fifth_char));
    }

    #[test]
    fn test_validate_and_resolve_special_amino_acids() {
        let aa = AminoAcidAmbiguity::new();

        // Test sequence with X, U, O
        let result = aa.validate_and_resolve("ACDEFXUO");
        assert!(result.is_ok());
        assert_eq!(result.unwrap().as_ref(), "ACDEFXUO");

        // Test sequence with X, U, O and ambiguous amino acids
        let result = aa.validate_and_resolve("ACDEFXBZJ");
        assert!(result.is_ok());
        let resolved = result.unwrap();
        assert!(resolved.as_ref().starts_with("ACDEFX"));
        // The B, Z, J should be resolved
        let seventh_char = resolved.as_ref().chars().nth(6).unwrap();
        let eighth_char = resolved.as_ref().chars().nth(7).unwrap();
        let ninth_char = resolved.as_ref().chars().nth(8).unwrap();
        assert!(vec!['D', 'N'].contains(&seventh_char));
        assert!(vec!['E', 'Q'].contains(&eighth_char));
        assert!(vec!['I', 'L'].contains(&ninth_char));
    }
}
