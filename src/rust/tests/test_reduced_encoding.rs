/// End-to-end tests for the multi-letter reduced alphabets in [`crate::alphabets`].
///
/// These share the pre-encoding path with the custom HP alphabets, so they are exposed to
/// the same two bugs that `test_hp_encoding` guards against: kmer_positions hashed in a
/// different case from the minhash, and encoded_sequence falling back to raw amino acids.
///
/// The test sequence is residues 121-181 of C. elegans CED-9 (UniProt P41958), taken
/// verbatim from `tests/testdata/fasta/ced9.fasta`, and the expected encodings below are
/// that fragment run through each paper's published partition.
#[cfg(test)]
mod tests {
    use crate::alphabets::ReducedAlphabet;
    use crate::sketch::ProteinSketch;

    /// CED-9 (P41958) residues 121-181, spanning the BH1 region.
    const CED9_FRAGMENT: &str = "TIFEKKHAENFETFCEQLLAVPRISFSLYQDVVRTVGNAQTDQCPMSYGRLIGLISFGGFV";

    /// Every hash in kmer_positions must also be in the minhash, and vice versa. A case
    /// mismatch between the two encoding loops empties matched regions without failing, so
    /// each alphabet is checked rather than a representative one.
    #[test]
    fn test_reduced_alphabets_kmer_positions_match_minhash() {
        let ksize = 10;

        for alphabet in ReducedAlphabet::all() {
            let moltype = alphabet.to_moltype();
            let sketch = ProteinSketch::from_protein_sequence(
                "ced9_fragment",
                CED9_FRAGMENT,
                ksize,
                1,
                &moltype,
            )
            .unwrap_or_else(|e| panic!("{moltype}: from_protein_sequence failed: {e}"));

            let minhash_set = sketch.mins_as_set();
            let kmer_positions = sketch.kmer_positions();

            assert!(!minhash_set.is_empty(), "{moltype}: minhash is empty");
            assert!(!kmer_positions.is_empty(), "{moltype}: kmer_positions is empty");

            for &hashval in kmer_positions.keys() {
                assert!(
                    minhash_set.contains(&hashval),
                    "{moltype}: kmer_positions hash {hashval} is not in the minhash"
                );
            }
            for &hashval in &minhash_set {
                assert!(
                    kmer_positions.contains_key(&hashval),
                    "{moltype}: minhash hash {hashval} is missing from kmer_positions"
                );
            }
        }
    }

    /// The stored encoded_sequence must be the fragment written in the reduced alphabet,
    /// not the raw amino acids. Spelled out for four alphabets spanning the size range so
    /// that a wrong cluster constant shows up as a diff on a real sequence.
    #[test]
    fn test_encoded_sequence_matches_published_partition() {
        let expected = [
            (
                ReducedAlphabet::Gbmr4,
                "ayyaaayaaayaayyaayyaypayayayyaayyaaygaaaaaaypyaygayygyyayggyy",
            ),
            (
                ReducedAlphabet::Sdm12,
                "tlykkkhaknyktycktllalpkltytlytdllktlgnattdtcpltygkllglltyggyl",
            ),
            (
                ReducedAlphabet::Hsdm17,
                "tlfkkkhaknfktfckqllalprlsfslyqdllrtlgnaqtdqcpmsygrllgllsfggfl",
            ),
            (
                ReducedAlphabet::Uniprot18,
                "tifekkhaenfetfceqhhaverisfshyqdvvrtvgnaqtdqcemsygrhighisfggfv",
            ),
        ];

        for (alphabet, want) in expected {
            let moltype = alphabet.to_moltype();
            let sketch = ProteinSketch::from_protein_sequence(
                "ced9_fragment",
                CED9_FRAGMENT,
                8,
                1,
                &moltype,
            )
            .unwrap();
            let got = sketch
                .get_moltype_sequence()
                .unwrap_or_else(|| panic!("{moltype}: encoded_sequence is None"));
            assert_eq!(got, want, "{moltype}");
        }
    }

    /// Every character of an encoded sequence has to be one of the alphabet's own symbols,
    /// so the class count in the name is what the index stores.
    #[test]
    fn test_encoded_sequence_uses_only_alphabet_symbols() {
        for alphabet in ReducedAlphabet::all() {
            let moltype = alphabet.to_moltype();
            let symbols: Vec<char> = alphabet
                .clusters()
                .iter()
                .map(|c| c.chars().next().unwrap().to_ascii_lowercase())
                .collect();

            let sketch = ProteinSketch::from_protein_sequence(
                "ced9_fragment",
                CED9_FRAGMENT,
                8,
                1,
                &moltype,
            )
            .unwrap();
            let encoded = sketch.get_moltype_sequence().unwrap();

            for ch in encoded.chars() {
                assert!(
                    symbols.contains(&ch),
                    "{moltype}: encoded char {ch:?} is not a class symbol"
                );
            }
        }
    }

    /// A coarser alphabet must collapse the fragment at least as hard as any finer one.
    /// Checked over every ordered pair rather than a sliding chain. A table wired to the
    /// wrong alphabet could otherwise sit between two neighbours that agree with each other.
    #[test]
    fn test_coarser_alphabets_collapse_at_least_as_much() {
        use std::collections::HashSet;

        let distinct = |alphabet: &ReducedAlphabet| -> usize {
            let sketch = ProteinSketch::from_protein_sequence(
                "ced9_fragment",
                CED9_FRAGMENT,
                8,
                1,
                &alphabet.to_moltype(),
            )
            .unwrap();
            sketch.get_moltype_sequence().unwrap().chars().collect::<HashSet<char>>().len()
        };

        for coarse in ReducedAlphabet::all() {
            // No alphabet can use more symbols on this fragment than it has classes.
            assert!(
                distinct(coarse) <= coarse.size(),
                "{}: {} distinct symbols exceeds its {} classes",
                coarse.name(),
                distinct(coarse),
                coarse.size()
            );

            for fine in ReducedAlphabet::all() {
                if coarse.size() >= fine.size() {
                    continue;
                }
                assert!(
                    distinct(coarse) <= distinct(fine),
                    "{} ({} classes) uses {} symbols, more than {} ({} classes) at {}",
                    coarse.name(),
                    coarse.size(),
                    distinct(coarse),
                    fine.name(),
                    fine.size(),
                    distinct(fine)
                );
            }
        }
    }

    /// A sequence searched against itself must yield matched regions under every alphabet.
    #[test]
    fn test_find_matched_regions_self_hit() {
        use crate::search::find_matched_regions;

        for alphabet in ReducedAlphabet::all() {
            let moltype = alphabet.to_moltype();
            let sketch = ProteinSketch::from_protein_sequence(
                "ced9_fragment",
                CED9_FRAGMENT,
                8,
                1,
                &moltype,
            )
            .unwrap();

            let shared_hashes = sketch.mins_as_set();
            assert!(!shared_hashes.is_empty(), "{moltype}: self-hit shares no hashes");

            let regions = find_matched_regions(&sketch, &sketch, &shared_hashes);
            assert!(!regions.is_empty(), "{moltype}: self-hit produced no matched regions");
        }
    }

    /// Two different alphabets must not produce the same sketch, otherwise a run labelled
    /// SDM12 could silently be HSDM17. CED-9 residues 121-181 contain Q and D, which
    /// HSDM17 keeps apart but SDM12 merges into its TSQ and acidic classes.
    #[test]
    fn test_alphabets_produce_distinct_sketches() {
        let sketches: Vec<(String, std::collections::HashSet<u64>)> = ReducedAlphabet::all()
            .iter()
            .map(|alphabet| {
                let moltype = alphabet.to_moltype();
                let sketch = ProteinSketch::from_protein_sequence(
                    "ced9_fragment",
                    CED9_FRAGMENT,
                    8,
                    1,
                    &moltype,
                )
                .unwrap();
                (moltype, sketch.mins_as_set())
            })
            .collect();

        for (i, (left_name, left)) in sketches.iter().enumerate() {
            for (right_name, right) in &sketches[i + 1..] {
                assert_ne!(left, right, "{left_name} and {right_name} produced identical sketches");
            }
        }
    }
}
