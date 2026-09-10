use std::borrow::Cow;

use crate::alphabets::canonical_moltype;
use crate::errors::{IndexError, IndexResult};

/// Standard amino acids and their properties
pub const STANDARD_AA: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y',
];

/// Amino acids that carry no residue identity and so can never be reduced: X (any residue) and
/// the stop codon.
pub const SPECIAL_AA: [char; 2] = ['X', '*'];

/// The two residues each ambiguous residue stands for.
///
///   - `B` (Asx) is `D` (Asp, aspartate) or `N` (Asn, asparagine)
///   - `J` (Xle) is `I` (Ile, isoleucine) or `L` (Leu, leucine)
///   - `Z` (Glx) is `E` (Glu, glutamate) or `Q` (Gln, glutamine)
///
/// These appear when the source method could not tell the pair apart, most often because
/// Asn and Gln deamidate to Asp and Glu during acid hydrolysis.
///
/// An ambiguous residue is never resolved to one of the pair. Every k-mer window covering
/// it is indexed under *both* readings instead (`disambiguate_kmer`), so a search matches whichever
/// residue the query holds. Picking one would assert a residue the source never
/// claimed, and which reading is safe depends on the alphabet: SDM12 and HSDM17 give Asp
/// and Asn separate classes, so under them the two readings are different k-mers.
pub const AMBIGUITY_ALTERNATIVES: [(char, [char; 2]); 3] =
    [('B', ['D', 'N']), ('J', ['I', 'L']), ('Z', ['E', 'Q'])];

/// The residues `residue` stands for, or `None` if it is not an ambiguous residue.
fn alternatives(residue: u8) -> Option<[u8; 2]> {
    AMBIGUITY_ALTERNATIVES
        .iter()
        .find(|(ambiguous, _)| *ambiguous as u8 == residue)
        .map(|(_, [first, second])| [*first as u8, *second as u8])
}

/// Whether `residues` contains any ambiguous residue, and so needs disambiguating.
pub fn has_ambiguous_residues(residues: &[u8]) -> bool {
    residues.iter().any(|b| alternatives(*b).is_some())
}

/// Ceiling on how many ambiguous residues one k-mer may carry before it is dropped rather
/// than expanded.
///
/// This bounds memory; it expresses no view about which readings are worth keeping. `2^n`
/// readings would otherwise grow without bound on pathological input, such as a long run of
/// `B`. The densest window in Swiss-Prot 2026_03 holds 9 of them, so no real sequence comes
/// near this.
pub const MAX_AMBIGUOUS_RESIDUES_PER_KMER: usize = 16;

/// Disambiguate one k-mer into every reading its ambiguous residues allow: `D` and `N` for
/// `B`, `I` and `L` for `J`, and `E` and `Q` for `Z`. A k-mer with no ambiguous residue yields
/// itself; one carrying `n` of them yields all `2^n` readings.
///
/// Indexing every reading keeps a search matching whichever residue the query holds, and
/// keeps matching from depending on which reading was kept. Picking one would assert a
/// residue the source never claimed, and which reading is safe depends on the alphabet:
/// SDM12 and HSDM17 give Asp and Asn separate classes, so under them the readings are
/// different k-mers.
///
/// The expansion is affordable because ambiguous residues are rare and stay sparse within
/// window. Swiss-Prot 2026_03 holds 525 of them, 276 `B` and 249 `Z` with no `J` anywhere,
/// across 146 of its 575_748 sequences. The densest window holds 9, so the worst single
/// k-mer expands to 512 readings, and expansion grows the index by 0.0012% at k=4 and
/// 0.034% at k=30.
///
/// WHY bytes rather than `&str`: callers hash the result, and hashing reads bytes. Going
/// through `String` would add a UTF-8 validation per reading and a panic path for input
/// that validation has already ruled out.
pub fn disambiguate_kmer(kmer: &[u8]) -> Option<Vec<Vec<u8>>> {
    let ambiguous: Vec<(usize, [u8; 2])> = kmer
        .iter()
        .enumerate()
        .filter_map(|(position, residue)| alternatives(*residue).map(|pair| (position, pair)))
        .collect();

    if ambiguous.len() > MAX_AMBIGUOUS_RESIDUES_PER_KMER {
        return None;
    }

    let mut readings = vec![kmer.to_vec()];
    for (position, [first, second]) in ambiguous {
        let mut expanded = Vec::with_capacity(readings.len() * 2);
        for reading in readings {
            let mut with_second = reading.clone();
            with_second[position] = second;
            let mut with_first = reading;
            with_first[position] = first;
            expanded.push(with_first);
            expanded.push(with_second);
        }
        readings = expanded;
    }
    Some(readings)
}

/// Non-canonical residues paired with their closest canonical analogue.
///
/// These are specific amino acids, not ambiguous residues, so they carry real chemistry that
/// would otherwise be discarded as unknown. U (Sec, selenocysteine) is cysteine with
/// selenium in place of sulfur; O (Pyl, pyrrolysine) is a lysine derivative. Each takes
/// whichever side its analogue takes in the active alphabet.
pub const NONCANONICAL_AA: [(char, char); 2] = [('U', 'C'), ('O', 'K')];

/// Resolves amino acids that are not one of the 20 canonical residues.
#[derive(Debug, Default)]
pub struct AminoAcidAmbiguity;

impl AminoAcidAmbiguity {
    pub fn new() -> Self {
        Self
    }

    fn is_valid_aa(&self, aa: char) -> bool {
        STANDARD_AA.contains(&aa)
            || SPECIAL_AA.contains(&aa)
            || alternatives(aa as u8).is_some()
            || Self::representative(aa).is_some()
    }

    /// The canonical residue used to encode `aa`, or `None` if `aa` needs no substitution.
    ///
    /// WHY deterministic: this previously drew at random from the alternatives, which made
    /// indexing non-reproducible — the same FASTA produced different k-mers, hashes and
    /// index contents on every run. The alternatives are interchangeable under every reduced
    /// alphabet, so a fixed representative is both reproducible and lossless there.
    fn representative(aa: char) -> Option<char> {
        NONCANONICAL_AA
            .iter()
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

    /// Validates a protein sequence, substituting canonical analogues for the non-canonical
    /// residues U and O when `moltype` reduces the alphabet. Stops at the first stop codon.
    ///
    /// U (Sec) and O (Pyl) are specific residues, not ambiguities: U is cysteine with
    /// selenium for sulfur, O is a lysine derivative. Under a reduced alphabet each lands in
    /// its analogue's class anyway, so substituting is lossless and keeps them from being
    /// discarded as unknown. Under `protein20` nothing is substituted, since there is no
    /// class to fall into and rewriting U as C would assert a residue the source never had.
    ///
    /// The ambiguous residues B, J and Z are never resolved here. They stay in the sequence and
    /// every k-mer covering one is indexed under both readings; see `disambiguate_kmer`.
    pub fn validate_and_resolve<'a>(
        &self,
        sequence: &'a str,
        moltype: &str,
    ) -> IndexResult<Cow<'a, str>> {
        // Canonicalized like every other entry point: sourmash's `protein` names the same
        // alphabet as `protein20`, and the two must not resolve U and O differently.
        let reduces_alphabet = canonical_moltype(moltype) != "protein20";

        // Validate first, recording where the kept region ends and where substitution first
        // becomes necessary. Almost every sequence needs neither (9_066 of Swiss-Prot
        // 2026_03's 209.0 M residues are non-canonical: 8_181 X, 331 U, 276 B, 249 Z, 29 O
        // and no J at all), so building an owned copy up front would allocate and copy once
        // per sequence only to discard it.
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

        // Reaching here means reduces_alphabet held. Amino acids with no representative -- X, the
        // stop codon, and the ambiguous residues, which are expanded at k-mer time instead --
        // copy through unchanged.
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

        // Test ambiguous residues
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
        // Canonical residues, X and the stop codon need no substitution.
        for c in STANDARD_AA.iter() {
            assert_eq!(AminoAcidAmbiguity::representative(*c), None, "{c}");
        }
        assert_eq!(AminoAcidAmbiguity::representative('X'), None);
        assert_eq!(AminoAcidAmbiguity::representative('*'), None);

        // Non-canonical residues resolve to one fixed analogue.
        // B, J and Z are expanded at k-mer time, not resolved to one residue.
        assert_eq!(AminoAcidAmbiguity::representative('B'), None);
        assert_eq!(AminoAcidAmbiguity::representative('J'), None);
        assert_eq!(AminoAcidAmbiguity::representative('Z'), None);
        assert_eq!(AminoAcidAmbiguity::representative('U'), Some('C'));
        assert_eq!(AminoAcidAmbiguity::representative('O'), Some('K'));
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

        assert_eq!(
            aa.validate_and_resolve("ACDEF*GHI", "hp_lehninger2").unwrap().as_ref(),
            "ACDEF*"
        );
        // B survives: it is expanded at k-mer time, not resolved here.
        assert_eq!(
            aa.validate_and_resolve("ACDEFB*GHI", "hp_lehninger2").unwrap().as_ref(),
            "ACDEFB*"
        );
    }

    #[test]
    fn test_validate_and_resolve_substitutes_noncanonical_in_reduced_alphabets() {
        let aa = AminoAcidAmbiguity::new();

        // U/O take their canonical analogue; B/J/Z and X are untouched.
        assert_eq!(
            aa.validate_and_resolve("ACDEFXBZJUO", "dayhoff6").unwrap().as_ref(),
            "ACDEFXBZJCK"
        );
    }

    #[test]
    fn test_validate_and_resolve_keeps_ambiguous_residues_verbatim_for_protein() {
        let aa = AminoAcidAmbiguity::new();

        // Under protein20 there is no equivalence to exploit, so nothing is substituted.
        let resolved = aa.validate_and_resolve("ACDEFXBZJUO", "protein20").unwrap();
        assert_eq!(resolved.as_ref(), "ACDEFXBZJUO");
        assert!(matches!(resolved, Cow::Borrowed(_)), "unchanged input should not allocate");
    }

    #[test]
    fn test_validate_and_resolve_is_deterministic() {
        let aa = AminoAcidAmbiguity::new();

        // Previously this drew at random, so repeated calls disagreed.
        let first = aa.validate_and_resolve("BZJUO", "hp_lehninger2").unwrap().into_owned();
        for _ in 0..10 {
            assert_eq!(aa.validate_and_resolve("BZJUO", "hp_lehninger2").unwrap().as_ref(), first);
        }
        assert_eq!(first, "BZJCK");
    }

    #[test]
    fn test_validate_and_resolve_no_ambiguous() {
        let aa = AminoAcidAmbiguity::new();

        // Test sequence with no ambiguous amino acids
        let result = aa.validate_and_resolve("ACDEFGHIKLMNPQRSTVWY", "hp_lehninger2");
        assert!(result.is_ok());
        // Should return borrowed string (no allocation)
        match result.unwrap() {
            Cow::Borrowed(s) => assert_eq!(s, "ACDEFGHIKLMNPQRSTVWY"),
            Cow::Owned(_) => panic!("Expected borrowed string for non-ambiguous sequence"),
        }
    }

    /// Disambiguation readable as strings; the function itself works in bytes so that
    /// hashing does not pay for UTF-8 validation.
    fn readings(kmer: &str) -> Option<Vec<String>> {
        Some(
            disambiguate_kmer(kmer.as_bytes())?
                .into_iter()
                .map(|reading| String::from_utf8(reading).expect("input is ASCII"))
                .collect(),
        )
    }

    /// B stands for Asp or Asn, so a k-mer covering one is indexed under both readings.
    /// Picking a single residue would commit to a reading the source never made, and under
    /// SDM12 or HSDM17 -- where Asp and Asn are separate classes -- the two readings are
    /// different k-mers.
    #[test]
    fn test_disambiguation_yields_both_residues_of_each_code() {
        assert_eq!(readings("MKBTA").unwrap(), vec!["MKDTA", "MKNTA"]);
        assert_eq!(readings("MKJTA").unwrap(), vec!["MKITA", "MKLTA"]);
        assert_eq!(readings("MKZTA").unwrap(), vec!["MKETA", "MKQTA"]);
    }

    /// A window free of ambiguous residues yields itself, so disambiguation adds no k-mers in
    /// the common case.
    #[test]
    fn test_disambiguation_passes_through_unambiguous_kmers() {
        assert_eq!(readings("MKTAY").unwrap(), vec!["MKTAY"]);
        // X and the stop codon carry no residue identity, so they are not expanded either.
        assert_eq!(readings("MXT*A").unwrap(), vec!["MXT*A"]);
    }

    /// Human BCL-2 (UniProt P10415) residues 1-30. The tests below write Asp10 as the
    /// ambiguous residue that stands for it, B, so the real fragment is one of the two
    /// readings that comes back.
    const BCL2_1_30: &str = "MAHAGRTGYDNREIVMKYIHYKLSQRGYEW";

    /// Every ambiguous residue in the window expands, so `n` of them give `2^n` readings and the real
    /// fragment is among them whichever residue the query holds.
    #[test]
    fn test_every_code_in_the_kmer_expands() {
        // Residues 1-20 with Asp10 written as B: one ambiguous residue, two readings.
        assert_eq!(
            readings("MAHAGRTGYBNREIVMKYIH").unwrap(),
            vec![&BCL2_1_30[..20], "MAHAGRTGYNNREIVMKYIH"]
        );
        // The same window with Glu13 also written as Z: two ambiguous residues, so four readings.
        let two = readings("MAHAGRTGYBNRZIVMKYIH").unwrap();
        assert_eq!(two.len(), 4);
        assert!(two.contains(&BCL2_1_30[..20].to_string()), "{two:?}");
        // Length does not decide it either: the same two ambiguous residues across thirty residues.
        let longer = readings("MAHAGRTGYBNRZIVMKYIHYKLSQRGYEW").unwrap();
        assert_eq!(longer.len(), 4);
        assert!(longer.contains(&BCL2_1_30.to_string()), "{longer:?}");
    }

    /// The ceiling bounds memory on pathological input. A k-mer past it is dropped whole,
    /// rather than indexed under an arbitrary subset of its readings.
    #[test]
    fn test_kmer_past_the_ceiling_is_dropped() {
        let at_ceiling = "B".repeat(MAX_AMBIGUOUS_RESIDUES_PER_KMER);
        assert_eq!(readings(&at_ceiling).unwrap().len(), 1 << MAX_AMBIGUOUS_RESIDUES_PER_KMER);
        assert_eq!(readings(&"B".repeat(MAX_AMBIGUOUS_RESIDUES_PER_KMER + 1)), None);
    }

    /// Every reading must be a sequence over the canonical residues, since each stands for a
    /// k-mer the source could have held.
    #[test]
    fn test_disambiguated_readings_are_canonical() {
        let aa = AminoAcidAmbiguity::new();
        for reading in readings("MAHAGRTGYBNREIVMKYIHYKLSQRGYEW").unwrap() {
            assert!(aa.validate_sequence(&reading).is_ok(), "{reading}");
            for c in reading.chars() {
                assert!(STANDARD_AA.contains(&c), "{reading}: {c} is not canonical");
            }
        }
    }

    /// sourmash's names must resolve identically to kmerseek's for the same alphabet.
    /// Before canonicalizing here, `protein` counted as a reducing alphabet and rewrote U
    /// and O to C and K, while `protein20` left them alone: two names, two sequences, two
    /// sets of k-mers.
    #[test]
    fn test_validate_and_resolve_treats_sourmash_names_identically() {
        let aa = AminoAcidAmbiguity::new();

        for (sourmash, kmerseek) in
            [("protein", "protein20"), ("dayhoff", "dayhoff6"), ("hp", "hp_lehninger2")]
        {
            assert_eq!(
                aa.validate_and_resolve("ACDEFXBZJUO", sourmash).unwrap(),
                aa.validate_and_resolve("ACDEFXBZJUO", kmerseek).unwrap(),
                "{sourmash} and {kmerseek} should resolve the same"
            );
        }

        // The full alphabet keeps U and O; a reduced one takes their analogues.
        assert_eq!(
            aa.validate_and_resolve("ACDEFXBZJUO", "protein").unwrap().as_ref(),
            "ACDEFXBZJUO"
        );
        assert_eq!(
            aa.validate_and_resolve("ACDEFXBZJUO", "dayhoff").unwrap().as_ref(),
            "ACDEFXBZJCK"
        );
    }
}
