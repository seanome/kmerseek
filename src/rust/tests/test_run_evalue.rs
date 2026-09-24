//! Checks of `region_run_evalue` (`evalue::run_evalue`) against real searches.

use std::io::Write;
use std::path::Path;

use tempfile::TempDir;

use crate::evalue::{run_evalue, SplitMix64};
use crate::index::ProteomeIndex;
use crate::search::{ProteinSearcher, SearchFilters, SearchResult, DEFAULT_BATCH_SIZE};
use crate::sketch::ProteinSketch;
use crate::tests::test_fixtures::{TEST_DECOYS_2MER_GZ, TEST_FASTA_GZ};

/// Amino-acid composition of UniProtKB/Swiss-Prot release 2026_03, in percent
/// (web.expasy.org/docs/relnotes/relstat.html).
const SWISSPROT_2026_03_PERCENT: [(u8, f64); 20] = [
    (b'A', 8.25),
    (b'R', 5.52),
    (b'N', 4.06),
    (b'D', 5.46),
    (b'C', 1.38),
    (b'E', 6.71),
    (b'Q', 3.93),
    (b'G', 7.07),
    (b'H', 2.27),
    (b'I', 5.90),
    (b'L', 9.64),
    (b'K', 5.79),
    (b'M', 2.41),
    (b'F', 3.86),
    (b'P', 4.75),
    (b'S', 6.66),
    (b'T', 5.36),
    (b'W', 1.10),
    (b'Y', 2.92),
    (b'V', 6.85),
];

/// A protein of `len` residues, each drawn on its own from the Swiss-Prot composition.
fn random_protein(rng: &mut SplitMix64, len: usize) -> String {
    let total: f64 = SWISSPROT_2026_03_PERCENT.iter().map(|&(_, p)| p).sum();
    (0..len)
        .map(|_| {
            let mut u = (rng.next_u64() >> 11) as f64 / (1u64 << 53) as f64 * total;
            for &(residue, p) in &SWISSPROT_2026_03_PERCENT {
                if u < p {
                    return residue as char;
                }
                u -= p;
            }
            'V'
        })
        .collect()
}

fn write_random_fasta(path: &Path, prefix: &str, n: usize, len: usize, rng: &mut SplitMix64) {
    let mut file = std::fs::File::create(path).unwrap();
    for i in 0..n {
        writeln!(file, ">{prefix}{i}\n{}", random_protein(rng, len)).unwrap();
    }
}

/// Every sketch of a FASTA, indexed at `moltype` and `ksize`, with sequences stored.
fn sketches(dir: &Path, fasta: &str, moltype: &str, ksize: u32) -> Vec<ProteinSketch> {
    let index = ProteomeIndex::new(dir, ksize, 1, moltype, true).unwrap();
    index.process_fasta(fasta, 0, DEFAULT_BATCH_SIZE).unwrap();
    index.load_state().unwrap();
    let mut out: Vec<_> = index.get_signatures().iter().map(|e| e.value().clone()).collect();
    out.sort_by(|a, b| a.signature().name.cmp(&b.signature().name));
    out
}

/// Exact search of every query in `query_fasta` against `target_fasta`, with every filter
/// open, and the target sketches for recounting compositions.
fn exact_search(
    query_fasta: &str,
    target_fasta: &str,
    moltype: &str,
    ksize: u32,
) -> (Vec<SearchResult>, Vec<ProteinSketch>, Vec<ProteinSketch>) {
    let dir = TempDir::new().unwrap();
    let targets = sketches(&dir.path().join("t_all"), target_fasta, moltype, ksize);
    let index = ProteomeIndex::new(dir.path().join("t"), ksize, 1, moltype, true).unwrap();
    index.process_fasta(target_fasta, 0, DEFAULT_BATCH_SIZE).unwrap();
    let searcher = ProteinSearcher::new(index).unwrap();
    let queries = sketches(&dir.path().join("q"), query_fasta, moltype, ksize);
    let results = searcher.search(&queries, &SearchFilters::default()).unwrap();
    (results, queries, targets)
}

/// Chance that two positions, one from each encoded sequence, share a class, counted here
/// rather than taken from the search.
fn pr_same(q: &[u8], t: &[u8]) -> f64 {
    let fraction = |s: &[u8], c: u8| s.iter().filter(|&&b| b == c).count() as f64 / s.len() as f64;
    let mut classes: Vec<u8> = q.iter().chain(t).copied().collect();
    classes.sort_unstable();
    classes.dedup();
    classes.iter().map(|&c| fraction(q, c) * fraction(t, c)).sum()
}

fn encoded(sketch: &ProteinSketch) -> &[u8] {
    sketch.get_moltype_sequence().unwrap().as_bytes()
}

/// Every (query, target) pair of encoded sequences.
fn pairs<'a>(
    queries: &'a [ProteinSketch],
    targets: &'a [ProteinSketch],
) -> impl Iterator<Item = (&'a [u8], &'a [u8])> {
    queries.iter().flat_map(move |q| targets.iter().map(move |t| (encoded(q), encoded(t))))
}

/// Length of every maximal run of equal letters on any diagonal of `q` against `t`, found
/// by walking every diagonal: what the search should report, independent of how it finds
/// regions.
fn run_lengths(q: &[u8], t: &[u8]) -> Vec<u32> {
    let mut lengths = Vec::new();
    for offset in -(q.len() as isize - 1)..t.len() as isize {
        let q_start = (-offset).max(0) as usize;
        let t_start = offset.max(0) as usize;
        let mut run = 0;
        for (&a, &b) in q[q_start..].iter().zip(&t[t_start..]) {
            if a == b {
                run += 1;
            } else if run > 0 {
                lengths.push(run);
                run = 0;
            }
        }
        if run > 0 {
            lengths.push(run);
        }
    }
    lengths
}

/// For each L0: regions the search reported with run_length >= L0, and the sum of E_run at
/// L0 over every query-target pair, including the pairs that share no k-mer at all.
fn counted_and_predicted(moltype: &str, ksize: u32, l0s: &[u32]) -> Vec<(u32, usize, f64)> {
    const N_QUERIES: usize = 10;
    const N_TARGETS: usize = 40;
    const LENGTH: usize = 1000;
    let dir = TempDir::new().unwrap();
    let (q_fasta, t_fasta) = (dir.path().join("q.fasta"), dir.path().join("t.fasta"));
    let mut rng = SplitMix64::new(20260924);
    write_random_fasta(&q_fasta, "query", N_QUERIES, LENGTH, &mut rng);
    write_random_fasta(&t_fasta, "target", N_TARGETS, LENGTH, &mut rng);
    let (results, queries, targets) =
        exact_search(q_fasta.to_str().unwrap(), t_fasta.to_str().unwrap(), moltype, ksize);
    let scanned_runs: Vec<u32> =
        pairs(&queries, &targets).flat_map(|(q, t)| run_lengths(q, t)).collect();
    l0s.iter()
        .map(|&l0| {
            let counted = results
                .iter()
                .flat_map(|r| &r.matched_regions)
                .filter(|region| region.run_length >= l0)
                .count();
            let predicted: f64 = pairs(&queries, &targets)
                .map(|(q, t)| run_evalue(pr_same(q, t), l0, q.len() as f64, t.len() as f64))
                .sum();
            let scanned = scanned_runs.iter().filter(|&&length| length >= l0).count();
            assert_eq!(counted, scanned, "{moltype} L0={l0}: the search missed or split runs");
            (l0, counted, predicted)
        })
        .collect()
}

/// Test 2 of the run E-value: on sequences whose letters are drawn independently, the
/// number of runs at least L0 long that an exact search reports matches the sum of E_run
/// over all pairs, within 3 Poisson standard deviations. Ten 1,000-residue queries against
/// forty 1,000-residue targets. The counts also equal a walk along every diagonal
/// (`run_lengths`), so the search neither misses nor splits a run.
///
/// The formula counts m x n places a run can start, but a run of L0 cannot start in the last
/// L0 - 1 positions of either sequence, so it predicts (1 - (L0 - 1) / 1000)^2 too many:
/// 2.2% at L0 = 12, 3.8% at L0 = 20. Where the counts are in the tens of thousands the
/// Poisson spread (0.4%) is smaller than that, and the counts fall outside 3 standard
/// deviations for this reason alone. The assertions use the L0 where the spread is 3.5% to
/// 13%, larger than the edge effect. The shorter L0 are listed for the record.
///
/// | alphabet      | k  | L0 | counted | predicted | spread | asserted |
/// |---------------|----|----|---------|-----------|--------|----------|
/// | hp_lehninger2 | 12 | 12 | 49,236  | 50,273.9  | 0.4%   | no       |
/// | hp_lehninger2 | 12 | 16 | 3,044   | 3,176.2   | 1.8%   | no       |
/// | hp_lehninger2 | 12 | 18 | 768     | 798.4     | 3.5%   | yes      |
/// | hp_lehninger2 | 12 | 19 | 392     | 400.3     | 5.0%   | yes      |
/// | hp_lehninger2 | 12 | 20 | 198     | 200.7     | 7.1%   | yes      |
/// | gbmr4         | 10 | 10 | 30,495  | 30,978.6  | 0.6%   | no       |
/// | gbmr4         | 10 | 13 | 2,039   | 2,132.4   | 2.2%   | no       |
/// | gbmr4         | 10 | 15 | 330     | 358.7     | 5.3%   | yes      |
/// | gbmr4         | 10 | 16 | 154     | 147.2     | 8.2%   | yes      |
/// | gbmr4         | 10 | 17 | 60      | 60.4      | 12.9%  | yes      |
///
/// The counts are pinned; the sequences are seeded. This checks the code against the
/// formula. It says nothing about real proteins, whose classes repeat in patterns: see
/// `region_run_evalue_on_2mer_shuffled_decoys`.
#[test]
fn region_run_evalue_matches_runs_between_random_sequences() {
    let cases = [
        ("hp_lehninger2", 12, [(18, 768), (19, 392), (20, 198)]),
        ("gbmr4", 10, [(15, 330), (16, 154), (17, 60)]),
    ];
    for (moltype, ksize, expected) in cases {
        let l0s: Vec<u32> = expected.iter().map(|&(l0, _)| l0).collect();
        let measured = counted_and_predicted(moltype, ksize, &l0s);
        for ((l0, counted, predicted), (_, expected_count)) in measured.into_iter().zip(expected) {
            let sd = predicted.sqrt();
            assert!(sd / predicted < 0.15, "{moltype} L0={l0}: predicted only {predicted:.1}");
            assert!(
                (counted as f64 - predicted).abs() <= 3.0 * sd,
                "{moltype} L0={l0}: counted {counted}, predicted {predicted:.1}"
            );
            assert_eq!(counted, expected_count, "{moltype} L0={l0}");
        }
    }
}

/// Query-target pairs whose best region has region_run_evalue at most `x`.
fn pairs_at_or_below(
    results: &[SearchResult],
    x: f64,
    evalue: impl Fn(&SearchResult) -> f64,
) -> usize {
    results.iter().filter(|r| evalue(r) <= x).count()
}

fn best_run_evalue(r: &SearchResult) -> f64 {
    r.matched_regions.iter().filter_map(|region| region.run_evalue).fold(f64::INFINITY, f64::min)
}

fn best_poisson_evalue(r: &SearchResult) -> f64 {
    r.matched_regions.iter().map(|region| region.poisson_evalue).fold(f64::INFINITY, f64::min)
}

/// Test 3 of the run E-value, on PR #79's decoys: the 25 BCL-2-like proteins against 500
/// sequences shuffled with their 2-mer counts kept, exact search at hp_lehninger2. A
/// calibrated E-value gives about 25 query-target pairs with a best region at E <= 1 (one
/// per query) and 250 at E <= 10. Pinned as measured, nothing tuned:
///
/// | k  | E <= | region_run_evalue | region_poisson_evalue | calibrated |
/// |----|------|-------------------|-----------------------|------------|
/// | 15 | 1    | 6                 | 2,995                 | 25         |
/// | 15 | 10   | 165               | 4,165                 | 250        |
/// | 12 | 1    | 6                 | 2,351                 | 25         |
/// | 12 | 10   | 175               | 4,858                 | 250        |
///
/// The run E-value calls fewer decoy pairs than a calibrated one would, so on these decoys
/// it errs toward too large. The Poisson E-value calls 94 to 120 times too many at E <= 1.
#[test]
fn region_run_evalue_on_2mer_shuffled_decoys() {
    // (k, pairs, run E-value at <= 1 and <= 10, Poisson E-value at <= 1 and <= 10)
    let expected = [(15, 8762, [6, 165], [2995, 4165]), (12, 12440, [6, 175], [2351, 4858])];
    for (ksize, n_pairs, run, poisson) in expected {
        let (results, queries, targets) =
            exact_search(TEST_FASTA_GZ, TEST_DECOYS_2MER_GZ, "hp_lehninger2", ksize);
        assert_eq!((queries.len(), targets.len(), results.len()), (25, 500, n_pairs));
        let count =
            |e: fn(&SearchResult) -> f64| [1.0, 10.0].map(|x| pairs_at_or_below(&results, x, e));
        assert_eq!(count(best_run_evalue), run, "k={ksize}");
        assert_eq!(count(best_poisson_evalue), poisson, "k={ksize}");
    }
}
