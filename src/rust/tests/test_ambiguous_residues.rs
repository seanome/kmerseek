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
        /// Readings sketched, summed over windows: a window holding n ambiguous residues
        /// contributes 2^n, and a window holding more than
        /// `MAX_AMBIGUOUS_RESIDUES_PER_KMER` (10) is dropped.
        readings: usize,
        /// Distinct hashes those readings produce. Equal to `readings` when the alphabet
        /// keeps D/N and E/Q apart, and equal to the number of windows kept when it merges
        /// both pairs, since then every reading of a window hashes the same.
        query_hashes: usize,
        /// Query hashes also found in bovine RNase A.
        shared_with_bovine: usize,
    }

    /// One row per alphabet, in `Alphabet::all()` order.
    ///
    /// Windows and readings per k: k=10 has 115 windows and at most 4 ambiguous residues
    /// per window, giving 541 readings; k=11 has 114 windows and 617 readings; k=12, 113
    /// and 709; k=14, 111 and 933; k=15, 110 and 1083; k=17, 108 and 1452; k=19, 106 and
    /// 1984; k=22, 103 and 3250; k=27, 98 and 6916. At k=43 the 11 windows holding 11 or 12
    /// ambiguous residues are dropped, leaving 71 of 82, with 28032 readings.
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
            readings: 1452,
            query_hashes: 108,
            shared_with_bovine: 74,
        },
        Expected {
            moltype: "hp_lehninger2",
            ksize: 43,
            readings: 28032,
            query_hashes: 71,
            shared_with_bovine: 52,
        },
        Expected {
            moltype: "hp_thomas_dill2",
            ksize: 43,
            readings: 28032,
            query_hashes: 71,
            shared_with_bovine: 52,
        },
        Expected {
            moltype: "hp_kyte_doolittle2",
            ksize: 43,
            readings: 28032,
            query_hashes: 71,
            shared_with_bovine: 52,
        },
        Expected {
            moltype: "hp_thomas_dill_no_c2",
            ksize: 43,
            readings: 28032,
            query_hashes: 71,
            shared_with_bovine: 52,
        },
        Expected {
            moltype: "hp_lehninger_c_nonpolar2",
            ksize: 43,
            readings: 28032,
            query_hashes: 71,
            shared_with_bovine: 52,
        },
        Expected {
            moltype: "hp_lehninger_hpc3",
            ksize: 27,
            readings: 6916,
            query_hashes: 98,
            shared_with_bovine: 79,
        },
        Expected {
            moltype: "hp_pbotc_1st_ed2",
            ksize: 43,
            readings: 28032,
            query_hashes: 71,
            shared_with_bovine: 52,
        },
        Expected {
            moltype: "gbmr4",
            ksize: 22,
            readings: 3250,
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
            readings: 1083,
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
            readings: 709,
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
