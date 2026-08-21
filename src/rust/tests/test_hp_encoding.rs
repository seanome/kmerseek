/// Regression tests for HP-alphabet encoding bugs:
///
/// Bug 1: kmer_positions were computed with lowercase h/p hashes that never matched
///        the uppercase H/P hashes stored by sourmash's ReadingFrame::new_protein.
///
/// Bug 2: encoded_sequence was storing raw amino acids instead of HP-encoded h/p codes,
///        because encode_by_moltype returned the identity function for custom HP alphabets.

#[cfg(test)]
mod tests {
    use crate::sketch::ProteinSketch;
    use crate::SEED;
    use sourmash::_hash_murmur;

    /// For a custom HP alphabet, every hash in kmer_positions must also be in the minhash.
    #[test]
    fn test_custom_hp_kmer_positions_match_minhash() {
        let seq = "NSQLAGKRILVTQADTFMGPTLCEVFAEMGNTLSGFLNYCSFNLNLQTLRHYVLAKILNKH";
        let ksize = 10;
        let sketch = ProteinSketch::from_protein_sequence(
            "test_seq",
            seq,
            ksize,
            1,
            "reduced_hp_kyte_doolittle2",
        )
        .unwrap();

        let minhash_set = sketch.mins_as_set();
        assert!(!minhash_set.is_empty(), "minhash should not be empty");

        let kmer_positions = sketch.kmer_positions();
        assert!(
            !kmer_positions.is_empty(),
            "kmer_positions should not be empty for a sequence with k={ksize}"
        );

        // Every hash in kmer_positions must be in the minhash (Bug 1 regression).
        for &hashval in kmer_positions.keys() {
            assert!(
                minhash_set.contains(&hashval),
                "hash {hashval} in kmer_positions is not in minhash — uppercase mismatch (Bug 1)"
            );
        }

        // Every minhash entry must appear in kmer_positions (completeness check).
        for &hashval in &minhash_set {
            assert!(
                kmer_positions.contains_key(&hashval),
                "minhash hash {hashval} missing from kmer_positions"
            );
        }
    }

    /// The encoded_sequence for a custom HP alphabet must consist solely of 'h' and 'p'
    /// characters (or pass-through bytes for unknown amino acids). Lowercase matches the
    /// built-in hp/dayhoff encodings, which sourmash emits lowercase.
    #[test]
    fn test_custom_hp_encoded_sequence_contains_hp_chars() {
        let seq = "MKTAYIAKQRFLVS";
        let sketch = ProteinSketch::from_protein_sequence(
            "test_hp_enc",
            seq,
            8,
            1,
            "reduced_hp_kyte_doolittle2",
        )
        .unwrap();

        let enc = sketch
            .get_moltype_sequence()
            .expect("encoded_sequence should be Some for custom HP alphabet (Bug 2)");

        assert_eq!(enc.len(), seq.len(), "encoded sequence should have same length as input");

        for ch in enc.chars() {
            assert!(
                ch == 'h' || ch == 'p',
                "encoded_sequence char {ch:?} is not 'h' or 'p' — raw amino acid stored instead (Bug 2)"
            );
        }
    }

    /// The encoded_sequence must differ from raw_sequence for custom HP alphabets.
    #[test]
    fn test_custom_hp_encoded_sequence_differs_from_raw() {
        let seq = "MKTAYIAKQRFLVS";
        let sketch = ProteinSketch::from_protein_sequence(
            "test_hp_diff",
            seq,
            8,
            1,
            "reduced_hp_kyte_doolittle2",
        )
        .unwrap();

        let raw = sketch.get_raw_sequence().expect("raw_sequence should be Some");
        let enc = sketch.get_moltype_sequence().expect("encoded_sequence should be Some");

        assert_ne!(
            raw, enc,
            "encoded_sequence must differ from raw_sequence for a custom HP alphabet (Bug 2)"
        );
    }

    /// Hashing uppercase H/P bytes must match what sourmash's ReadingFrame stores.
    /// This verifies the upstream invariant that our fix depends on.
    #[test]
    fn test_uppercase_hp_hash_matches_sourmash() {
        // sourmash uppercases h/p to H/P before hashing — confirm hash of uppercase differs
        // from lowercase so that the fix (uppercasing) is actually necessary.
        let lower: Vec<u8> = b"ppphpphpphpppppppphppp".to_vec();
        let upper: Vec<u8> = b"PPPHPPHPPHPPPPPPPPHPPP".to_vec();
        let hash_lower = _hash_murmur(&lower, SEED);
        let hash_upper = _hash_murmur(&upper, SEED);
        assert_ne!(
            hash_lower, hash_upper,
            "uppercase and lowercase HP k-mers must hash differently; if equal the fix is unnecessary"
        );
    }

    /// End-to-end: two identical sequences with a custom HP alphabet must produce
    /// at least one matched region (find_matched_regions succeeds).
    #[test]
    fn test_custom_hp_find_matched_regions_self_hit() {
        use crate::search::find_matched_regions;

        let seq = "MKTAYIAKQRFLVSNSQLAGKRILVTQAD";

        let sketch = ProteinSketch::from_protein_sequence(
            "self_test",
            seq,
            8,
            1,
            "reduced_hp_kyte_doolittle2",
        )
        .unwrap();

        let shared_hashes = sketch.mins_as_set();
        assert!(!shared_hashes.is_empty(), "self-hit must share hashes");

        let regions = find_matched_regions(&sketch, &sketch, &shared_hashes);

        assert!(
            !regions.is_empty(),
            "find_matched_regions should return at least one region for a self-hit with custom HP alphabet"
        );
    }

    /// All custom HP alphabets must produce kmer_positions that match the minhash.
    /// Catches per-alphabet regressions of Bug 1.
    #[test]
    fn test_all_custom_hp_alphabets_kmer_positions_match_minhash() {
        let seq = "MKTAYIAKQRFLVSNSQLAGKRILVTQADTFMGPTLCEVFAEMG";
        let ksize = 8;

        let alphabets = [
            "reduced_hp_kyte_doolittle2",
            "reduced_hp_thomas_dill2",
            "reduced_hp_lehninger2",
            "reduced_hp_thomas_dill_no_c2",
            "reduced_hp_lehninger_c_nonpolar2",
            "reduced_hp_lehninger_hpc3",
            "reduced_hp_pbotc_1st_ed2",
        ];

        for moltype in alphabets {
            let sketch = ProteinSketch::from_protein_sequence("test", seq, ksize, 1, moltype)
                .unwrap_or_else(|e| panic!("{moltype}: from_protein_sequence failed: {e}"));

            let minhash_set = sketch.mins_as_set();
            let kmer_positions = sketch.kmer_positions();

            assert!(!minhash_set.is_empty(), "{moltype}: minhash is empty");
            assert!(!kmer_positions.is_empty(), "{moltype}: kmer_positions is empty (Bug 1)");

            for &h in kmer_positions.keys() {
                assert!(
                    minhash_set.contains(&h),
                    "{moltype}: kmer_positions hash {h} not in minhash (Bug 1 regression)"
                );
            }
            for &h in &minhash_set {
                assert!(
                    kmer_positions.contains_key(&h),
                    "{moltype}: minhash hash {h} missing from kmer_positions"
                );
            }

            let enc = sketch
                .get_moltype_sequence()
                .unwrap_or_else(|| panic!("{moltype}: encoded_sequence is None (Bug 2)"));
            for ch in enc.chars() {
                assert!(
                    ch == 'h' || ch == 'p' || ch == 'c',
                    "{moltype}: encoded_sequence char {ch:?} is not h/p/c (Bug 2 regression)"
                );
            }
        }
    }

    /// For standard encodings (dayhoff, hp), kmer_positions hashes must match the minhash.
    /// Guards against accidentally uppercasing dayhoff/hp codes in the else branch (which
    /// would produce hash mismatches since sourmash hashes those as lowercase).
    #[test]
    fn test_standard_encodings_kmer_positions_match_minhash() {
        let seq = "MKTAYIAKQRFLVSNSQLAGKRILVTQAD";
        let ksize = 8;

        for moltype in ["hp", "dayhoff"] {
            let sketch = ProteinSketch::from_protein_sequence("test", seq, ksize, 1, moltype)
                .unwrap_or_else(|e| panic!("{moltype}: from_protein_sequence failed: {e}"));

            let minhash_set = sketch.mins_as_set();
            let kmer_positions = sketch.kmer_positions();

            assert!(!minhash_set.is_empty(), "{moltype}: minhash is empty");
            assert!(!kmer_positions.is_empty(), "{moltype}: kmer_positions is empty");

            for &h in kmer_positions.keys() {
                assert!(
                    minhash_set.contains(&h),
                    "{moltype}: kmer_positions hash {h} not in minhash — uppercase applied to {moltype} codes by mistake"
                );
            }
            for &h in &minhash_set {
                assert!(
                    kmer_positions.contains_key(&h),
                    "{moltype}: minhash hash {h} missing from kmer_positions"
                );
            }
        }
    }

    /// find_matched_regions with an empty intersection must return an empty vec without panic.
    /// This guards against any accidental assumption that intersection is non-empty.
    #[test]
    fn test_find_matched_regions_empty_intersection_returns_empty() {
        use crate::search::find_matched_regions;
        use std::collections::HashSet;

        let seq = "MKTAYIAKQRFLVSNSQLAGKRILVTQAD";
        let sketch =
            ProteinSketch::from_protein_sequence("test", seq, 8, 1, "reduced_hp_kyte_doolittle2")
                .unwrap();

        let empty: HashSet<u64> = HashSet::new();
        let regions = find_matched_regions(&sketch, &sketch, &empty);
        assert!(regions.is_empty(), "empty intersection must yield no matched regions");
    }

    /// find_matched_regions with a disjoint query/target (no shared sequence) must not panic.
    /// The soft-skip (`if query_moltype_seq != target_moltype_seq`) guards hash collisions
    /// where two different HP sequences share the same hash. Simulate by using two completely
    /// different proteins that share zero HP k-mers.
    #[test]
    fn test_find_matched_regions_disjoint_sequences_returns_empty() {
        use crate::search::find_matched_regions;

        let seq_a = "MKTAYIAKQRFLVSNSQLAGKRILVTQAD";
        let seq_b = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"; // all-cys: likely different HP pattern

        let sketch_a =
            ProteinSketch::from_protein_sequence("a", seq_a, 8, 1, "reduced_hp_kyte_doolittle2")
                .unwrap();
        let sketch_b =
            ProteinSketch::from_protein_sequence("b", seq_b, 8, 1, "reduced_hp_kyte_doolittle2")
                .unwrap();

        let shared = sketch_a
            .mins_as_set()
            .intersection(&sketch_b.mins_as_set())
            .cloned()
            .collect::<std::collections::HashSet<_>>();

        // Must not panic even if shared hashes map to different HP subsequences.
        let _regions = find_matched_regions(&sketch_a, &sketch_b, &shared);
    }

    /// An index built before the HP rename stores `hp_<name>`; a query sketched afterwards
    /// uses `reduced_hp_<name>2`. Both must resolve to the same table and therefore the same
    /// hashes, otherwise every pre-rename index would silently stop matching.
    #[test]
    fn test_legacy_and_current_moltype_names_sketch_identically() {
        use crate::hp_alphabets::HpAlphabet;

        let seq = "TIFEKKHAENFETFCEQLLAVPRISFSLYQDVVRTVGNAQTDQCPMSYGRLIGLISFGGFV";

        for alphabet in HpAlphabet::all_named() {
            let legacy = alphabet.legacy_moltype();
            let current = alphabet.to_moltype();

            let from_legacy =
                ProteinSketch::from_protein_sequence("x", seq, 8, 1, &legacy).unwrap();
            let from_current =
                ProteinSketch::from_protein_sequence("x", seq, 8, 1, &current).unwrap();

            assert_eq!(
                from_legacy.mins_as_set(),
                from_current.mins_as_set(),
                "{legacy} and {current} produced different hashes"
            );
            // The legacy spelling is normalized away, so nothing downstream sees two names.
            assert_eq!(from_legacy.moltype().get(), current, "{legacy} was not normalized");
        }
    }

    /// Renaming `dayhoff` to `reduced_dayhoff6` must not change a single hash: both names
    /// resolve to sourmash's Murmur64Dayhoff and its own encoder, so existing dayhoff
    /// indexes stay searchable. This is the guarantee that `hp` cannot make, because the
    /// built-in HP path hashes lowercase h/p while the custom tables hash uppercase H/P.
    #[test]
    fn test_dayhoff_rename_preserves_hashes() {
        let seq = "TIFEKKHAENFETFCEQLLAVPRISFSLYQDVVRTVGNAQTDQCPMSYGRLIGLISFGGFV";

        let legacy = ProteinSketch::from_protein_sequence("x", seq, 8, 1, "dayhoff").unwrap();
        let current =
            ProteinSketch::from_protein_sequence("x", seq, 8, 1, "reduced_dayhoff6").unwrap();

        assert_eq!(legacy.mins_as_set(), current.mins_as_set());
        assert_eq!(legacy.get_moltype_sequence(), current.get_moltype_sequence());
        assert_eq!(legacy.moltype().get(), "reduced_dayhoff6");
    }

    /// The built-in `hp` and our Lehninger table agree on the partition but not on the
    /// bytes they hash: sourmash encodes `hp` to lowercase h/p, while a pre-encoded custom
    /// table is uppercased to H/P before hashing. They are therefore separate moltypes, and
    /// merging their names would silently invalidate every index built under the other one.
    #[test]
    fn test_builtin_hp_and_lehninger2_share_a_partition_but_not_hashes() {
        let seq = "TIFEKKHAENFETFCEQLLAVPRISFSLYQDVVRTVGNAQTDQCPMSYGRLIGLISFGGFV";

        let builtin = ProteinSketch::from_protein_sequence("x", seq, 8, 1, "hp").unwrap();
        let custom =
            ProteinSketch::from_protein_sequence("x", seq, 8, 1, "reduced_hp_lehninger2").unwrap();

        // Same partition: the encoded sequences are character-for-character equal.
        assert_eq!(builtin.get_moltype_sequence(), custom.get_moltype_sequence());

        // Different hashed bytes: no k-mer hash is shared.
        let (b, c) = (builtin.mins_as_set(), custom.mins_as_set());
        assert_eq!(b.len(), c.len());
        assert_eq!(b.intersection(&c).count(), 0, "hp and lehninger2 unexpectedly share hashes");
    }
}
