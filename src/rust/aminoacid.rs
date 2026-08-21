use std::borrow::Cow;

use crate::encoding::custom_alphabet_table;
use crate::errors::{IndexError, IndexResult};
use crate::hp_alphabets::HpAlphabet;

/// Standard amino acids and their properties
pub const STANDARD_AA: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y',
];

/// Codes that carry no residue identity and so can never be reduced: X (any residue) and
/// the stop codon.
pub const SPECIAL_AA: [char; 2] = ['X', '*'];

/// Ambiguity codes paired with the representative used to encode them.
///
/// Each code stands for one of two residues that land on the *same* side of every reduced
/// alphabet, so either representative yields an identical encoding and the choice is
/// lossless:
///   - B (Asx) = Asp or Asn — Dayhoff `c`, polar in every HP table
///   - J (Xle) = Ile or Leu — Dayhoff `e`, hydrophobic in every HP table
///   - Z (Glx) = Glu or Gln — Dayhoff `c`, polar in every HP table
pub const AMBIGUITY_CODES: [(char, char); 3] = [('B', 'D'), ('J', 'I'), ('Z', 'E')];

/// The two residues each ambiguity code stands for.
///
/// Substituting a fixed representative is lossless only in an alphabet that maps both
/// residues of a pair onto the same symbol, which is what `preserves_ambiguity_equivalence`
/// checks.
pub const AMBIGUITY_ALTERNATIVES: [(char, [char; 2]); 3] =
    [('B', ['D', 'N']), ('J', ['I', 'L']), ('Z', ['E', 'Q'])];

/// Whether both residues behind every ambiguity code encode identically under `moltype`.
///
/// Dayhoff and the built-in `hp` have no lookup table of ours (sourmash encodes them
/// internally) and are covered by
/// `test_ambiguity_alternatives_encode_identically_in_reduced_alphabets`, so they answer
/// `true`. Most of the multi-letter `reduced_*` alphabets answer `false`: SDM12 and HSDM17
/// give Asp and Asn their own classes, so reading B as Asp would assert a residue the
/// source never claimed.
fn preserves_ambiguity_equivalence(moltype: &str) -> bool {
    let Some(table) = custom_alphabet_table(moltype) else {
        return true;
    };
    AMBIGUITY_ALTERNATIVES
        .iter()
        .all(|(_, [first, second])| table.get(&(*first as u8)) == table.get(&(*second as u8)))
}

/// Non-canonical residues paired with their closest canonical analogue.
///
/// These are specific amino acids, not ambiguity codes, so they carry real chemistry that
/// would otherwise be discarded as unknown. U (Sec, selenocysteine) is cysteine with
/// selenium in place of sulfur; O (Pyl, pyrrolysine) is a lysine derivative. Each takes
/// whichever side its analogue takes in the active alphabet.
pub const NONCANONICAL_AA: [(char, char); 2] = [('U', 'C'), ('O', 'K')];

/// Resolves amino acid codes that are not one of the 20 canonical residues.
#[derive(Debug, Default)]
pub struct AminoAcidAmbiguity;

impl AminoAcidAmbiguity {
    pub fn new() -> Self {
        Self
    }

    fn is_valid_aa(&self, aa: char) -> bool {
        STANDARD_AA.contains(&aa) || SPECIAL_AA.contains(&aa) || Self::representative(aa).is_some()
    }

    /// The canonical residue used to encode `aa`, or `None` if `aa` needs no substitution.
    ///
    /// WHY deterministic: this previously drew at random from the alternatives, which made
    /// indexing non-reproducible — the same FASTA produced different k-mers, hashes and
    /// index contents on every run. The alternatives are interchangeable under every reduced
    /// alphabet, so a fixed representative is both reproducible and lossless there.
    fn representative(aa: char) -> Option<char> {
        AMBIGUITY_CODES
            .iter()
            .chain(NONCANONICAL_AA.iter())
            .find(|(code, _)| *code == aa)
            .map(|(_, representative)| *representative)
    }

    /// Validates a protein sequence and returns an error if invalid characters are found
    /// Stops reading at the first stop codon (*)
    pub fn validate_sequence(&self, sequence: &str) -> IndexResult<()> {
        for (i, c) in sequence.chars().enumerate() {
            if c == '*' {
                // Stop codon - this is valid, but we stop reading here
                return Ok(());
            }
            if !self.is_valid_aa(c) {
                return Err(IndexError::InvalidAminoAcid(c, i + 1));
            }
        }
        Ok(())
    }

    /// Validates a protein sequence, substituting representatives for non-canonical codes
    /// when `moltype` reduces the alphabet. Stops processing at the first stop codon (*).
    ///
    /// WHY only for alphabets that collapse each ambiguous pair: under Dayhoff or a named HP
    /// table a code like B encodes identically whether it is read as Asp or Asn, so
    /// substituting a representative is lossless. Three groups of moltypes are excluded:
    ///   - `protein`/`raw` keep the full 20-letter alphabet, so there is no such equivalence
    ///     to exploit — picking Asp would assert a residue the source never claimed.
    ///   - `reduced_hp_shuffled_control2[_1..10]` are HP tables too, but their partition is randomized
    ///     rather than biochemically derived, so the two alternatives can land on opposite
    ///     sides (see `test_shuffled_control_does_not_preserve_ambiguity_equivalence`). These
    ///     are excluded by name rather than by the equivalence check below, because a control
    ///     should not start substituting on the seeds where the shuffle happens to agree.
    ///   - the `reduced_*` alphabets that split an ambiguous pair across classes, which is
    ///     most of them (`preserves_ambiguity_equivalence`).
    ///
    /// All of them keep the original code and hash it as itself, consistent with how X is
    /// already handled.
    pub fn validate_and_resolve<'a>(
        &self,
        sequence: &'a str,
        moltype: &str,
    ) -> IndexResult<Cow<'a, str>> {
        let reduces_alphabet = !matches!(moltype, "protein" | "raw")
            && !matches!(
                HpAlphabet::from_moltype(moltype),
                Some(HpAlphabet::ShuffledControl | HpAlphabet::Shuffled(_))
            )
            && preserves_ambiguity_equivalence(moltype);

        // Validate first, recording where the kept region ends and where substitution first
        // becomes necessary. Almost every sequence needs neither (roughly 900 of SwissProt's
        // 207.6 M residues are non-canonical), so building an owned copy up front would
        // allocate and copy once per sequence only to discard it.
        let mut end = sequence.len();
        let mut first_substitution = None;
        for (offset, c) in sequence.char_indices() {
            if c == '*' {
                end = offset + c.len_utf8();
                break;
            }
            if !self.is_valid_aa(c) {
                return Err(IndexError::InvalidAminoAcid(c, offset + 1));
            }
            if first_substitution.is_none() && reduces_alphabet && Self::representative(c).is_some()
            {
                first_substitution = Some(offset);
            }
        }

        let kept = &sequence[..end];
        let Some(start) = first_substitution else {
            return Ok(if end == sequence.len() {
                Cow::Borrowed(sequence)
            } else {
                Cow::Owned(kept.to_string())
            });
        };

        // Reaching here means reduces_alphabet held, so every remaining code can be looked up
        // unconditionally. The stop codon has no representative and copies through.
        let mut result = String::with_capacity(kept.len());
        result.push_str(&kept[..start]);
        for c in kept[start..].chars() {
            result.push(Self::representative(c).unwrap_or(c));
        }
        Ok(Cow::Owned(result))
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
    fn test_representative_is_deterministic_and_exact() {
        // Canonical residues and identity-free codes need no substitution.
        for c in STANDARD_AA.iter() {
            assert_eq!(AminoAcidAmbiguity::representative(*c), None, "{c}");
        }
        assert_eq!(AminoAcidAmbiguity::representative('X'), None);
        assert_eq!(AminoAcidAmbiguity::representative('*'), None);

        // Ambiguity codes and non-canonical residues resolve to one fixed representative.
        assert_eq!(AminoAcidAmbiguity::representative('B'), Some('D'));
        assert_eq!(AminoAcidAmbiguity::representative('J'), Some('I'));
        assert_eq!(AminoAcidAmbiguity::representative('Z'), Some('E'));
        assert_eq!(AminoAcidAmbiguity::representative('U'), Some('C'));
        assert_eq!(AminoAcidAmbiguity::representative('O'), Some('K'));
    }

    /// The substitution must be lossless: both residues an ambiguity code stands for have to
    /// encode identically under every biochemically-derived alphabet, otherwise picking a
    /// representative would silently commit to one reading.
    ///
    /// `hp_shuffled_control` is deliberately exempt — it is a negative control whose partition
    /// is randomized, so it has no reason to respect biochemical equivalence, and is asserted
    /// separately below so that a real alphabet breaking equivalence still fails this test.
    #[test]
    fn test_ambiguity_alternatives_encode_identically_in_reduced_alphabets() {
        use crate::encoding::encode_by_moltype;
        use crate::hp_alphabets::HpAlphabet;

        for (code, [first, second]) in AMBIGUITY_ALTERNATIVES {
            let pair: String = [first, second].iter().collect();
            let pair = pair.as_str();
            for moltype in ["dayhoff", "hp"] {
                let encoded: Vec<String> = pair
                    .chars()
                    .map(|c| encode_by_moltype(&c.to_string(), moltype).unwrap())
                    .collect();
                assert_eq!(
                    encoded[0], encoded[1],
                    "{code}: {moltype} encodes {pair} differently ({encoded:?})"
                );
            }

            for alphabet in HpAlphabet::all_named() {
                if matches!(alphabet, HpAlphabet::ShuffledControl) {
                    continue;
                }
                let table = alphabet.table();
                let codes: Vec<u8> =
                    pair.bytes().map(|b| *table.get(&b).expect("canonical residue")).collect();
                assert_eq!(
                    codes[0],
                    codes[1],
                    "{code}: {} encodes {pair} differently",
                    alphabet.name()
                );
            }
        }
    }

    /// Documents the one alphabet where substituting a representative is *not* lossless.
    /// The shuffled control randomizes the partition, so D/N, I/L and E/Q land on opposite
    /// sides. Substitution there is still deterministic, which is what matters for a control.
    #[test]
    fn test_shuffled_control_does_not_preserve_ambiguity_equivalence() {
        use crate::hp_alphabets::HpAlphabet;

        let table = HpAlphabet::ShuffledControl.table();
        for (code, pair) in [('B', "DN"), ('J', "IL"), ('Z', "EQ")] {
            let codes: Vec<u8> = pair.bytes().map(|b| *table.get(&b).unwrap()).collect();
            assert_ne!(
                codes[0], codes[1],
                "{code}: shuffled control unexpectedly preserves {pair} equivalence"
            );
        }
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
        let error = result.unwrap_err();
        match error {
            IndexError::InvalidAminoAcid(c, pos) => {
                assert_eq!(c, '1');
                assert!(pos > 0);
            }
            _ => panic!("Expected InvalidAminoAcid error, got {:?}", error),
        }
    }

    #[test]
    fn test_validate_and_resolve_with_stop_codon() {
        let aa = AminoAcidAmbiguity::new();

        assert_eq!(aa.validate_and_resolve("ACDEF*GHI", "hp").unwrap().as_ref(), "ACDEF*");
        assert_eq!(aa.validate_and_resolve("ACDEFB*GHI", "hp").unwrap().as_ref(), "ACDEFD*");
    }

    #[test]
    fn test_validate_and_resolve_substitutes_in_reduced_alphabets() {
        let aa = AminoAcidAmbiguity::new();

        // B/J/Z take their representative; U/O take their canonical analogue; X is untouched.
        assert_eq!(
            aa.validate_and_resolve("ACDEFXBZJUO", "dayhoff").unwrap().as_ref(),
            "ACDEFXDEICK"
        );
    }

    #[test]
    fn test_validate_and_resolve_keeps_codes_verbatim_for_protein() {
        let aa = AminoAcidAmbiguity::new();

        // Under `protein` there is no equivalence to exploit, so nothing is substituted.
        let resolved = aa.validate_and_resolve("ACDEFXBZJUO", "protein").unwrap();
        assert_eq!(resolved.as_ref(), "ACDEFXBZJUO");
        assert!(matches!(resolved, Cow::Borrowed(_)), "unchanged input should not allocate");
    }

    /// `raw` is a synonym for `protein` in encoding.rs (get_hash_function_from_moltype and
    /// get_encoding_fn_from_moltype both treat them identically), so it must be exempt from
    /// substitution for the same reason `protein` is.
    #[test]
    fn test_validate_and_resolve_keeps_codes_verbatim_for_raw() {
        let aa = AminoAcidAmbiguity::new();

        let resolved = aa.validate_and_resolve("ACDEFXBZJUO", "raw").unwrap();
        assert_eq!(resolved.as_ref(), "ACDEFXBZJUO");
        assert!(matches!(resolved, Cow::Borrowed(_)), "unchanged input should not allocate");
    }

    /// The shuffled-control HP alphabets randomize the h/p partition, so a fixed representative
    /// is not lossless there (test_shuffled_control_does_not_preserve_ambiguity_equivalence).
    /// Both the base control and a seeded variant must keep codes verbatim, not substitute.
    #[test]
    fn test_validate_and_resolve_keeps_codes_verbatim_for_shuffled_control() {
        let aa = AminoAcidAmbiguity::new();

        for moltype in ["reduced_hp_shuffled_control2", "reduced_hp_shuffled_control2_1"] {
            let resolved = aa.validate_and_resolve("ACDEFXBZJUO", moltype).unwrap();
            assert_eq!(resolved.as_ref(), "ACDEFXBZJUO", "{moltype}");
            assert!(matches!(resolved, Cow::Borrowed(_)), "{moltype}: should not allocate");
        }
    }

    #[test]
    fn test_validate_and_resolve_is_deterministic() {
        let aa = AminoAcidAmbiguity::new();

        // Previously this drew at random, so repeated calls disagreed.
        let first = aa.validate_and_resolve("BZJUO", "hp").unwrap().into_owned();
        for _ in 0..10 {
            assert_eq!(aa.validate_and_resolve("BZJUO", "hp").unwrap().as_ref(), first);
        }
        assert_eq!(first, "DEICK");
    }

    #[test]
    fn test_validate_and_resolve_no_ambiguous() {
        let aa = AminoAcidAmbiguity::new();

        // Test sequence with no ambiguous amino acids
        let result = aa.validate_and_resolve("ACDEFGHIKLMNPQRSTVWY", "hp");
        assert!(result.is_ok());
        // Should return borrowed string (no allocation)
        match result.unwrap() {
            Cow::Borrowed(s) => assert_eq!(s, "ACDEFGHIKLMNPQRSTVWY"),
            Cow::Owned(_) => panic!("Expected borrowed string for non-ambiguous sequence"),
        }
    }

    /// Only GBMR4 and GBMR7 keep both residues of every ambiguous pair in one class. The
    /// other six reduced alphabets split at least one pair, so a fixed representative would
    /// commit to a reading the source never made. Spelled out per alphabet so that changing
    /// a cluster constant surfaces here instead of silently changing indexed hashes.
    #[test]
    fn test_which_reduced_alphabets_preserve_ambiguity_equivalence() {
        use crate::reduced_alphabets::ReducedAlphabet;

        let expected = [
            (ReducedAlphabet::Gbmr4, true),
            (ReducedAlphabet::Wwmj5, false),
            (ReducedAlphabet::Gbmr7, true),
            (ReducedAlphabet::Sdm12, false),
            (ReducedAlphabet::Mmseqs12, false),
            (ReducedAlphabet::Wass14, false),
            (ReducedAlphabet::Hsdm17, false),
            (ReducedAlphabet::Uniprot18, false),
        ];

        for (alphabet, preserves) in expected {
            assert_eq!(
                preserves_ambiguity_equivalence(&alphabet.to_moltype()),
                preserves,
                "{}",
                alphabet.name()
            );
        }
    }

    /// GBMR4 maps D/N, I/L and E/Q each onto a single class, so B/J/Z take their
    /// representative and U/O take their canonical analogue, exactly as under dayhoff.
    #[test]
    fn test_validate_and_resolve_substitutes_in_gbmr4() {
        let aa = AminoAcidAmbiguity::new();

        assert_eq!(
            aa.validate_and_resolve("ACDEFXBZJUO", "reduced_gbmr4").unwrap().as_ref(),
            "ACDEFXDEICK"
        );
    }

    /// SDM12 gives Asp and Asn separate classes, so nothing is substituted and B/J/Z/U/O
    /// hash as themselves, the same way X already does.
    #[test]
    fn test_validate_and_resolve_keeps_codes_verbatim_for_split_reduced_alphabets() {
        let aa = AminoAcidAmbiguity::new();

        for moltype in [
            "reduced_wwmj5",
            "reduced_sdm12",
            "reduced_mmseqs12",
            "reduced_wass14",
            "reduced_hsdm17",
            "reduced_uniprot18",
        ] {
            let resolved = aa.validate_and_resolve("ACDEFXBZJUO", moltype).unwrap();
            assert_eq!(resolved.as_ref(), "ACDEFXBZJUO", "{moltype}");
            assert!(matches!(resolved, Cow::Borrowed(_)), "{moltype}: should not allocate");
        }
    }
}
