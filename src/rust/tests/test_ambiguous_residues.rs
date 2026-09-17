//! Searching with a protein dense in ambiguous residues, under every alphabet.
//!
//! Topi pancreatic ribonuclease (UniProt P00659, RNAS1_DAMKO) was sequenced by Edman
//! degradation, and 22 of its 124 residues came out as B (Asp or Asn) or Z (Glu or Gln):
//! 12 B and 10 Z, never more than 17 residues apart. It is the query here, against the 125
//! UniProtKB ribonucleases it was downloaded with.
//!
//! Two of those targets make the expansion checkable by hand. Goat RNase (P67926,
//! RNAS1_CAPHI) is the same 124 residues and agrees with topi at every ambiguous position,
//! differing at one real position: residue 103 is K in topi and E in goat. Bovine RNase A
//! (P61823, RNAS1_BOVIN) carries a 26-residue signal peptide and then a mature chain that
//! differs from topi at four real positions (3, 19, 37 and 103) and, again, agrees at every
//! ambiguous one. So a shared k-mer count equal to "every window not covering a real
//! mismatch" means every B and Z was expanded to the residue the homolog holds.
//!
//! Each alphabet is run at the k where a k-mer carries as many bits as protein20 at k=10.

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use tempfile::TempDir;

    use crate::alphabets::Alphabet;
    use crate::index::ProteomeIndex;
    use crate::search::{ProteinSearcher, SearchFilters, SearchResult, DEFAULT_BATCH_SIZE};
    use crate::sketch::ProteinSketch;

    const RIBONUCLEASES_FASTA_GZ: &str =
        "tests/testdata/fasta/ribonuclease_125_entries_uniprotkb_2026_09_17.fasta.gz";
    const RNAS1_DAMKO_FASTA: &str = "tests/testdata/fasta/rnas1_damko_P00659.fasta";

    /// The k-mer size at which `alphabet` carries as many bits per k-mer as protein20 at
    /// k=10. A k-mer over N classes carries k·log2(N) bits, and 10·log2(20) is 43.2, so a
    /// two-class alphabet needs k=43 and dayhoff6 needs k=17. Rounded to the nearest k.
    fn equivalent_ksize(alphabet: Alphabet) -> u32 {
        let bits = 10.0 * 20f64.log2();
        (bits / (alphabet.size() as f64).log2()).round() as u32
    }

    /// What indexing and searching the topi query gives under one alphabet.
    struct Expected {
        moltype: &'static str,
        ksize: u32,
        /// Readings sketched, summed over windows. A window holding n residues that are
        /// still ambiguous after encoding contributes 2^n. Under an alphabet that merges
        /// Asp with Asn and Glu with Gln, no residue of this protein is ambiguous once
        /// encoded, so every window contributes one reading and none reaches
        /// `MAX_AMBIGUOUS_RESIDUES_PER_KMER`.
        readings: usize,
        /// Distinct hashes those readings produce, equal to `readings` for every alphabet
        /// since no two windows of the query encode the same.
        query_hashes: usize,
        /// Query hashes also found in bovine RNase A.
        shared_with_bovine: usize,
    }

    /// One row per alphabet, in `Alphabet::all()` order.
    ///
    /// Windows per k: 115 at k=10, 114 at k=11, 113 at k=12, 111 at k=14, 110 at k=15,
    /// 108 at k=17, 106 at k=19, 103 at k=22, 98 at k=27 and 82 at k=43. Where the alphabet
    /// keeps D/N and E/Q apart, the readings are the sum of 2^n over windows holding n
    /// ambiguous residues: at most 4 per window at k=10, giving 541; 617 at k=11; 709 at
    /// k=12; 933 at k=14; 1984 at k=19; 3250 at k=22. Where it merges both pairs (dayhoff6,
    /// gbmr4, gbmr7, mmseqs12 and the HP family), the readings are the windows. The densest
    /// window, 12 ambiguous residues at k=43, is under the two-class alphabets, where it is
    /// not ambiguous at all.
    ///
    /// The bovine count is the number of query hashes also found in bovine RNase A, which is
    /// the windows that cover none of the four real mismatches, less any the alphabet merges
    /// away. protein20, uniprot18, wass14 and hsdm17 keep all four apart. gbmr4 puts S, T,
    /// A, Q, K and N in one class, so it sees no mismatch at all and every one of its 103
    /// windows is shared. The HP alphabets merge three of the four (S/T, Q/K and K/N are
    /// each polar/polar) and keep only residue 19's S/A.
    const EXPECTED: &[Expected] = &[
        Expected {
            moltype: "protein20",
            ksize: 10,
            readings: 541,
            query_hashes: 541,
            shared_with_bovine: 82,
        },
        Expected {
            moltype: "dayhoff6",
            ksize: 17,
            readings: 108,
            query_hashes: 108,
            shared_with_bovine: 74,
        },
        Expected {
            moltype: "hp_lehninger2",
            ksize: 43,
            readings: 82,
            query_hashes: 82,
            shared_with_bovine: 63,
        },
        Expected {
            moltype: "hp_thomas_dill2",
            ksize: 43,
            readings: 82,
            query_hashes: 82,
            shared_with_bovine: 63,
        },
        Expected {
            moltype: "hp_kyte_doolittle2",
            ksize: 43,
            readings: 82,
            query_hashes: 82,
            shared_with_bovine: 63,
        },
        Expected {
            moltype: "hp_thomas_dill_no_c2",
            ksize: 43,
            readings: 82,
            query_hashes: 82,
            shared_with_bovine: 63,
        },
        Expected {
            moltype: "hp_lehninger_c_nonpolar2",
            ksize: 43,
            readings: 82,
            query_hashes: 82,
            shared_with_bovine: 63,
        },
        Expected {
            moltype: "hp_lehninger_hpc3",
            ksize: 27,
            readings: 98,
            query_hashes: 98,
            shared_with_bovine: 79,
        },
        Expected {
            moltype: "hp_pbotc_1st_ed2",
            ksize: 43,
            readings: 82,
            query_hashes: 82,
            shared_with_bovine: 63,
        },
        Expected {
            moltype: "gbmr4",
            ksize: 22,
            readings: 103,
            query_hashes: 103,
            shared_with_bovine: 103,
        },
        Expected {
            moltype: "polarity4",
            ksize: 22,
            readings: 3250,
            query_hashes: 3250,
            shared_with_bovine: 44,
        },
        Expected {
            moltype: "wwmj5",
            ksize: 19,
            readings: 1984,
            query_hashes: 1984,
            shared_with_bovine: 87,
        },
        Expected {
            moltype: "gbmr7",
            ksize: 15,
            readings: 110,
            query_hashes: 110,
            shared_with_bovine: 77,
        },
        Expected {
            moltype: "funcgroups8",
            ksize: 14,
            readings: 933,
            query_hashes: 933,
            shared_with_bovine: 69,
        },
        Expected {
            moltype: "sdm12",
            ksize: 12,
            readings: 709,
            query_hashes: 709,
            shared_with_bovine: 77,
        },
        Expected {
            moltype: "mmseqs12",
            ksize: 12,
            readings: 113,
            query_hashes: 113,
            shared_with_bovine: 89,
        },
        Expected {
            moltype: "wass14",
            ksize: 11,
            readings: 617,
            query_hashes: 617,
            shared_with_bovine: 78,
        },
        Expected {
            moltype: "hsdm17",
            ksize: 11,
            readings: 617,
            query_hashes: 617,
            shared_with_bovine: 78,
        },
        Expected {
            moltype: "uniprot18",
            ksize: 10,
            readings: 541,
            query_hashes: 541,
            shared_with_bovine: 82,
        },
    ];

    /// Index the 125 ribonucleases and sketch the topi query under one alphabet, then
    /// search. Returns the query sketch and every hit.
    fn search_topi(
        temp_dir: &TempDir,
        moltype: &str,
        ksize: u32,
    ) -> Result<(ProteinSketch, Vec<SearchResult>)> {
        let targets = ProteomeIndex::new(temp_dir.path().join("targets"), ksize, 1, moltype, true)?;
        targets.process_fasta(RIBONUCLEASES_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(targets)?;

        let query_index =
            ProteomeIndex::new(temp_dir.path().join("query"), ksize, 1, moltype, true)?;
        query_index.process_fasta(RNAS1_DAMKO_FASTA, 0, DEFAULT_BATCH_SIZE)?;
        query_index.load_state()?;
        let queries: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();
        assert_eq!(queries.len(), 1, "{moltype}: the query file holds one protein");

        let results = searcher.search(&queries, &SearchFilters::default())?;
        Ok((queries.into_iter().next().unwrap(), results))
    }

    fn shared_with<'a>(results: &'a [SearchResult], entry_name: &str) -> &'a SearchResult {
        results
            .iter()
            .find(|r| r.target_name.contains(entry_name))
            .unwrap_or_else(|| panic!("{entry_name} should be a hit"))
    }

    /// Every alphabet, at its protein20-k=10-equivalent k, expands all 22 ambiguous
    /// residues: the query sketches the expected number of readings, and the shared k-mer
    /// count with bovine RNase A is the windows that hold no real mismatch.
    #[test]
    fn test_every_alphabet_expands_the_22_ambiguous_residues_of_topi_ribonuclease() -> Result<()> {
        let moltypes: Vec<&str> = EXPECTED.iter().map(|e| e.moltype).collect();
        let all: Vec<&str> = Alphabet::all().iter().map(Alphabet::to_moltype).collect();
        assert_eq!(moltypes, all, "every alphabet must have a row here");

        for expected in EXPECTED {
            let moltype = expected.moltype;
            let alphabet = Alphabet::from_moltype(moltype).unwrap();
            assert_eq!(equivalent_ksize(alphabet), expected.ksize, "{moltype}");

            let temp_dir = TempDir::new()?;
            let (query, results) = search_topi(&temp_dir, moltype, expected.ksize)?;

            let readings: usize = query.kmer_positions().values().map(Vec::len).sum();
            assert_eq!(readings, expected.readings, "{moltype}: readings sketched");
            assert_eq!(query.mins_as_set().len(), expected.query_hashes, "{moltype}: query hashes");

            let bovine = shared_with(&results, "RNAS1_BOVIN");
            assert_eq!(
                bovine.n_intersecting_hashes, expected.shared_with_bovine,
                "{moltype}: k-mers shared with bovine RNase A"
            );
        }
        Ok(())
    }

    /// The spans of a hit's matched regions, as `(start, end)` in query coordinates,
    /// zero-based and end-exclusive, in order along the query.
    fn region_spans(hit: &SearchResult) -> Vec<(u32, u32)> {
        let mut spans: Vec<_> = hit.matched_regions.iter().map(|r| (r.start, r.end)).collect();
        spans.sort_unstable();
        spans
    }

    /// A matched region runs through an ambiguous residue instead of stopping at it. Every
    /// region below covers several B and Z, and each break is a real mismatch the alphabet
    /// can see. Bovine differs from topi at residues 3 (S/T), 19 (S/A), 37 (Q/K) and 103
    /// (K/N): sdm12 merges S with T and so breaks at the other three, giving four regions;
    /// every HP alphabet merges all but S/A, giving one region from residue 20 to the end.
    /// Goat differs only at residue 103 (K/E), which sdm12 merges, so under sdm12 goat is
    /// one region covering the whole protein.
    #[test]
    fn test_matched_regions_run_through_ambiguous_residues() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let (_, results) = search_topi(&temp_dir, "sdm12", 12)?;
        let bovine = shared_with(&results, "RNAS1_BOVIN");
        assert_eq!(region_spans(bovine), vec![(0, 18), (19, 36), (37, 102), (103, 124)]);
        let goat = shared_with(&results, "RNAS1_CAPHI");
        assert_eq!(region_spans(goat), vec![(0, 124)]);

        let temp_dir = TempDir::new()?;
        let (_, results) = search_topi(&temp_dir, "hp_lehninger2", 43)?;
        let bovine = shared_with(&results, "RNAS1_BOVIN");
        assert_eq!(region_spans(bovine), vec![(19, 124)]);
        assert_eq!(bovine.matched_regions[0].n_shared, 63);
        Ok(())
    }

    /// Under protein20 at k=10, goat RNase shares the 105 windows that do not cover
    /// residue 103, its one real difference from topi: 115 windows minus the 10 that cover
    /// it. Every other window covers at least one B or Z and matches only because that
    /// residue was expanded to the one goat holds.
    #[test]
    fn test_goat_ribonuclease_matches_every_window_off_its_one_real_mismatch() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let (query, results) = search_topi(&temp_dir, "protein20", 10)?;
        assert_eq!(query.kmer_positions().len(), 541);

        let goat = shared_with(&results, "RNAS1_CAPHI");
        assert_eq!(goat.n_intersecting_hashes, 105);
        assert_eq!(goat.ksize, 10);
        assert_eq!(goat.moltype, "protein20");
        Ok(())
    }
}
