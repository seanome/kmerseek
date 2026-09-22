use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::path::Path;

use anyhow::Result;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use statrs::distribution::{DiscreteCDF, Poisson};

use crate::aminoacid::encoded_residues_agree;
use crate::errors::{IndexError, IndexResult};
use crate::hash_functions::residue_encoder;
use crate::index::{ProteomeIndex, SearchCache};
use crate::karlin_altschul::{
    fit_scores, fit_scores_with_reference, make_decoy, survival_counts, DecoyNull, KaCalibration,
    SplitMix64, BIN_WIDTH, MIN_FIT_POINTS,
};
use crate::significance;
use crate::sketch::ProteinSketch;
use crate::types::MolType;

/// Default progress interval for FASTA processing (log progress every N sequences)
pub const DEFAULT_PROGRESS_INTERVAL: u32 = 1000;

/// Default batch size for FASTA processing (process N sequences per batch)
pub const DEFAULT_BATCH_SIZE: usize = 1000;

/// Result-level filters applied while a search is running, so that results failing the
/// filters are never allocated into the results `Vec` in the first place (as opposed to
/// building the full unfiltered `Vec` and then filtering it down afterward).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SearchFilters {
    /// Minimum containment score (query-side) required to keep a match.
    pub threshold: f64,
    /// Minimum number of shared k-mers required to keep a match.
    pub min_shared_kmers: usize,
    /// Maximum whole-query Poisson p-value required to keep a match.
    pub max_query_pvalue: f64,
    /// Minimum region-scoped score required to keep a match, applied to the best region in
    /// the pair. Bigger means more surprising (see `MatchedRegion::poisson_score`); a
    /// heuristic cutoff on a ranking score, not a statistically calibrated significance
    /// threshold.
    pub min_region_score: f64,
    /// Drop a target whose sketch is the query's own (same md5). On for all-vs-all searches
    /// of one index against itself, where every query would otherwise hit itself; off when a
    /// FASTA is searched against an index, where the query's identical entry is a real hit.
    pub skip_self_matches: bool,
}

impl Default for SearchFilters {
    /// Accepts every result `compare()` produces.
    ///
    /// `max_query_pvalue` must be `f64::INFINITY`, not 1.0: `query_poisson_pvalue` is 1.0
    /// whenever there's no database frequency context (e.g. `set_query_frequencies` was never
    /// called), which is common, so a cap of 1.0 combined with `compare()`'s strict-less-than
    /// keep check would wrongly reject those results.
    ///
    /// `min_region_score` must be `f64::NEG_INFINITY`, not 0.0, for the mirror-image reason:
    /// `poisson_score` is 0.0 whenever there's no evidence at all (`poisson_pvalue` of 1.0), so
    /// a floor of 0.0 combined with the strict-greater-than keep check would wrongly reject
    /// those results too.
    fn default() -> Self {
        Self {
            threshold: 0.0,
            min_shared_kmers: 0,
            max_query_pvalue: f64::INFINITY,
            min_region_score: f64::NEG_INFINITY,
            skip_self_matches: false,
        }
    }
}

/// How far past its exact seed a matched region may grow, and at what cost per mismatch.
///
/// A region from `find_matched_regions` is a maximal exact run in the encoded alphabet: one
/// class flip ends it. Between remote homologs the HP pattern is conserved per column far
/// better than any 23-residue stretch of it is conserved exactly (2024-kmerseek-analysis
/// notebooks 230 and 232), so an exact run is a seed, not the match. `extend_regions` walks
/// outward from each seed along the stored encoded sequences, +1 per agreeing position and
/// `-mismatch_penalty` per disagreeing one, and stops when the running score has fallen
/// `xdrop` below its best (Altschul's X-drop). Two seeds on one diagonal whose extensions
/// meet become one region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtensionParams {
    /// Score subtracted per disagreeing encoded position. Must be positive.
    pub mismatch_penalty: f64,
    /// Extension stops once the running score is this far below its best so far.
    pub xdrop: f64,
    /// Karlin-Altschul K for the E-value, fitted on calibration searches of the index
    /// (`KaCalibration`), or set by hand with `--ka-k`.
    pub ka_k: f64,
    /// Multiplier on every pair's closed-form lambda: the fitted lambda over the
    /// closed-form lambda at the database's own composition (`KaCalibration::lambda_scale`).
    /// 1.0 keeps the closed form, which assumes independent positions.
    pub ka_lambda_scale: f64,
    /// Chain colinear regions at most this many residues apart (on the query) into one
    /// region scored with Karlin-Altschul sum statistics. 0 leaves every region on its own.
    /// See `chain_regions`.
    pub chain_max_gap: u32,
    /// How far apart in diagonal two chained regions may sit, i.e. the largest net indel a
    /// chain tolerates between members. 0 chains only on one diagonal.
    pub chain_max_shift: u32,
}

/// Targets whose encoded sequences give the database's match probability during a fit.
const COMPOSITION_SAMPLE: usize = 200;

/// What `ProteinSearcher::calibrate_ka` found.
#[derive(Debug, Clone, PartialEq)]
pub struct KaCalibrationReport {
    /// The fit, or None when the calibration queries gave too few usable score bins.
    pub fitted: Option<KaCalibration>,
    /// Chance that two positions drawn from the sampled database sequences share a class:
    /// the `a` of the database against itself.
    pub match_probability: f64,
    pub n_queries: usize,
    pub n_regions: usize,
    /// Reference queries searched and the regions they produced: `n_queries` x
    /// `reference_shuffles` shuffles of the same sequences, and the chance curve the fit
    /// is read against. Zero under a null that has no reference.
    pub n_reference_queries: usize,
    pub n_reference_regions: usize,
    /// Residues in the calibration queries and k-mers in the database: the L and N the
    /// fit's K is read against.
    pub query_residues: u64,
    pub database_kmers: u64,
    pub lambda_analytic: f64,
    /// Regions with score >= each half-nat bin of x = lambda_region S, for the calibration
    /// queries (`survival`) and, under the database null, for the same queries shuffled
    /// (`reference_survival`, empty otherwise). Kept whether or not the fit succeeded, so
    /// a refused fit can still show the histogram it refused. Same shape as
    /// `KaCalibration::survival`.
    pub survival: Vec<(i64, u64)>,
    pub reference_survival: Vec<(i64, u64)>,
}

/// What a calibration run is asked for; see `ProteinSearcher::calibrate_ka`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KaCalibrationSettings {
    pub mismatch_penalty: f64,
    pub xdrop: f64,
    /// What the calibration queries are.
    pub null: DecoyNull,
    /// For `DecoyNull::Database`: what the reference queries are.
    pub reference: DecoyNull,
    /// For `DecoyNull::Database`: how many times each reference query is shuffled. The
    /// reference only has to hold `MIN_BIN_COUNT` regions per score bin for that bin to
    /// be usable, and above about 40 bits per seed (k x log2 classes) a chance k-mer
    /// match is rare enough that one shuffle per query does not reach it however many
    /// queries are searched: in the 0.4 dark-set store 21 of 140 indexes refused their
    /// fit with 6 to 2_655 chance regions in total, protein20 at k >= 9 and uniprot18 at
    /// k >= 10 among them, against a median 1.7M for the indexes that fitted. Shuffling
    /// each query n times multiplies the chance curve by n and leaves the real curve
    /// alone; the fit reads a ratio whose baseline is a fitted line, so the constant
    /// ln n it adds is absorbed by that line's intercept. 1 is the old behaviour.
    pub reference_shuffles: usize,
    pub n_queries: usize,
    pub seed: u64,
}

/// The lambda scale and K a search runs with.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KaParams {
    pub k: f64,
    pub lambda_scale: f64,
}

/// Where the lambda and K in use came from; see `ProteinSearcher::resolve_ka`.
#[derive(Debug, Clone, PartialEq)]
pub enum KaSource {
    /// `--ka-k`, with the closed-form lambda per region.
    Flag,
    Index(KaCalibration),
    Fitted(KaCalibration),
}

impl Display for KaSource {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            KaSource::Flag => write!(f, "--ka-k, closed-form lambda"),
            KaSource::Index(c) | KaSource::Fitted(c) => write!(
                f,
                "{}: {} {} queries, {} regions; slope {:.3} per nat of lambda_region S (1 = closed form holds; closed form {:.3} at the database's match probability {:.3}), K {:.4}, fit on x {:.1}..{:.1}{}, rms {:.3}{}",
                if matches!(self, KaSource::Index(_)) { "stored in the index" } else { "fitted now" },
                c.n_queries,
                c.null,
                c.n_regions,
                c.slope,
                c.lambda_analytic,
                c.match_probability,
                c.k,
                c.x_range().0,
                c.x_range().1,
                c.bend_score.map_or(String::new(), |b| format!(
                    ", {} x {:.1}",
                    if c.null == DecoyNull::Database { "relatives from" } else { "counts rise above the line from" },
                    b as f64 * c.bin_width
                )),
                c.rms_residual,
                c.reference_lambda.map_or(String::new(), |l| format!(
                    "; {} reference slope {:.3} over the same bins",
                    c.reference.map_or("shuffled".to_string(), |r| r.to_string()),
                    l / c.bin_width
                ))
            ),
        }
    }
}

/// Karlin-Altschul sum statistic for `r` segments whose normalised scores add to `t`:
/// P(T_r >= t) ~ e^-t t^(r-1) / (r! (r-1)!) (Karlin & Altschul 1993, PNAS 90:5873). Each
/// normalised score is lambda S_i - ln(K m n) and must be positive for the approximation to
/// hold; a non-positive sum returns 1.
pub fn karlin_altschul_sum_p(t: f64, r: u32) -> f64 {
    if t <= 0.0 || r == 0 {
        return 1.0;
    }
    let r_f = r as f64;
    let ln_fact = |x: f64| libm_lgamma(x + 1.0);
    let ln_p = -t + (r_f - 1.0) * t.ln() - ln_fact(r_f) - ln_fact(r_f - 1.0);
    ln_p.exp().min(1.0)
}

/// ln Gamma, Lanczos approximation; enough for the factorials of small chain lengths.
fn libm_lgamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln()
            - libm_lgamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + G + 0.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Class frequencies of an encoded sequence, keyed by byte. Gaps and unknowns count too,
/// since a match against them is also a match in the run.
fn class_composition(encoded: &[u8]) -> HashMap<u8, f64> {
    let mut counts: HashMap<u8, f64> = HashMap::new();
    for &b in encoded {
        *counts.entry(b).or_insert(0.0) += 1.0;
    }
    let n = encoded.len().max(1) as f64;
    counts.values_mut().for_each(|v| *v /= n);
    counts
}

/// The chance match rate `u`: the probability that one position drawn from `p` and one
/// drawn from `q` carry the same class.
fn match_probability(p: &HashMap<u8, f64>, q: &HashMap<u8, f64>) -> f64 {
    p.iter().map(|(b, pb)| pb * q.get(b).copied().unwrap_or(0.0)).sum()
}

/// The slice of `encoded` a region covers, clamped to the sequence. An out-of-range span
/// gives an empty slice rather than panicking.
fn region_span(encoded: &[u8], start: u32, end: u32) -> &[u8] {
    let lo = (start as usize).min(encoded.len());
    let hi = (end as usize).clamp(lo, encoded.len());
    &encoded[lo..hi]
}

/// The chance match rate `u` of one region: the two spans the region covers are counted
/// on their own, not as part of the proteins they sit in.
///
/// This is the whole point of scoring a region rather than a pair. An RS domain
/// (RNRDRDHKRRHRSRSRSRS...) inside an ordinary protein matches another polar-rich stretch
/// at three positions in four for free, but the two proteins around it look ordinary, so
/// the pair's `u` is ordinary and the free matches get scored as evidence. Counted over
/// the spans themselves, `u` lands past the boundary and `karlin_altschul_lambda` returns
/// 0: not assessable.
///
/// The cost is a noisier `u`. A 40-residue span gives 40 draws per sequence, so the
/// estimate wobbles where a whole protein's would not, and the wobble runs both ways.
/// Erring toward the boundary is the safe direction, since it withholds an E-value rather
/// than inventing one.
fn region_match_probability(q_span: &[u8], t_span: &[u8]) -> f64 {
    match_probability(&class_composition(q_span), &class_composition(t_span))
}

/// The Karlin-Altschul lambda for +1 / -penalty scoring when a random pair of positions
/// matches with probability `u`: the positive root of u e^x + (1-u) e^(-penalty x) = 1.
///
/// A positive root exists only when the expected score u - penalty (1-u) is negative,
/// i.e. u < penalty / (1 + penalty). Above that, agreement is what these two compositions
/// do by default and no run of it is surprising: returns 0. The left side is convex with
/// value 1 at x = 0 and slope u - penalty (1-u) there, so bisection on [0, hi] with hi
/// pushed out until f(hi) > 1 is safe.
///
/// `u` is what the E-value explainer writes for Pr(match | unrelated); the calibration
/// report calls the database-wide version `match_probability`.
pub fn karlin_altschul_lambda(u: f64, penalty: f64) -> f64 {
    let one_minus_u = 1.0 - u;
    if u <= 0.0 || u.is_nan() || u - penalty * one_minus_u >= 0.0 {
        return 0.0;
    }
    let f = |x: f64| u * x.exp() + one_minus_u * (-penalty * x).exp();
    let mut hi = 1.0;
    while f(hi) <= 1.0 {
        hi *= 2.0;
        if hi > 1e6 {
            return 0.0;
        }
    }
    let (mut lo, mut hi) = (0.0, hi);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid) > 1.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let root = 0.5 * (lo + hi);
    // At the boundary u = penalty / (1 + penalty) the root is 0 up to rounding; a lambda of
    // 1e-8 would make every E-value ~ K m n, which is the same "no evidence" answer.
    if root < 1e-6 {
        0.0
    } else {
        root
    }
}

impl SearchFilters {
    /// A pair is kept when either scope clears its cap. Both are not required.
    ///
    /// Requiring both would bring back the problem this PR fixes: a real sub-protein domain
    /// match diluted into insignificance by the rest of the protein. BCL2/CED9 at k=15 has a
    /// whole-query p-value of 0.99 and a region score of 3.16 (p=0.0007). Requiring both to
    /// pass would discard the sub-protein domain match that region scoring exists to surface.
    fn scopes_pass(&self, query_pvalue: f64, best_region_score: Option<f64>) -> bool {
        let query_passes = query_pvalue < self.max_query_pvalue;
        let region_passes = best_region_score.is_some_and(|score| score > self.min_region_score);
        query_passes || region_passes
    }
}

/// CSV-friendly version of SearchResult with matched region information
/// WHY: Each matched region gets its own row in the CSV, with all SearchResult similarity
/// metrics repeated for each region. This makes it easy to analyze individual matched regions
/// while still having access to the overall similarity metrics. Every CSV row must have matched
/// region data - SearchResults without matched regions are not included in the CSV output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResultCsv {
    // SearchResult fields
    pub query_name: String,
    pub query_md5: String,
    pub target_name: String,
    pub target_md5: String,
    pub containment: f64,
    pub n_intersecting_hashes: usize,
    pub ksize: u32,
    pub scaled: u32,
    pub moltype: String,
    /// Whether low-complexity (homopolymer) k-mers were removed from both the
    /// index and the query sketches for this search. Recorded per row, like
    /// ksize/scaled/moltype, so a results file is self-describing.
    pub remove_low_complexity: bool,
    pub jaccard: f64,
    pub max_containment: f64,
    pub average_abund: f64,
    pub median_abund: f64,
    pub std_abund: f64,
    pub containment_target_in_query: f64,
    pub f_weighted_target_in_query: f64,
    pub query_tfidf: f64,
    pub mean_matched_kmer_freq: f64,
    pub sum_matched_kmer_freq: f64,
    pub query_expected_shared_kmers: f64,
    pub query_enrichment: f64,
    /// Joint k-mer frequency: Σ freq_query[h] * freq_target[h] for h in intersection.
    /// Dot product of query and target frequency vectors over shared k-mers.
    /// 0.0 when query frequencies are not available (one-pass mode).
    pub joint_kmer_freq: f64,
    /// Poisson p-value: P(X ≥ n_intersecting_hashes | λ = query_expected_shared_kmers).
    /// 1.0 when query_expected_shared_kmers is unavailable (no database context).
    pub query_poisson_pvalue: f64,
    // How many other things this result was tested alongside: candidate region placements,
    // targets, k-mers, and queries. Reported as separate numbers rather than multiplied into
    // a p-value, so no reported statistic changes depending on batch composition. See
    // SearchResult for what each one counts and how to combine them.
    pub region_search_space: usize,
    pub db_n_targets: usize,
    pub db_n_kmers: usize,
    pub run_n_queries: usize,
    // MatchedRegion fields (always present - every CSV row has a matched region)
    pub region_start: u32,
    pub region_end: u32,
    pub region_subseq: String,
    pub target_start: u32,
    pub target_end: u32,
    pub target_subseq: String,
    pub moltype_seq: String,
    pub region_length: u32,
    pub region_n_shared_kmers: u32,
    pub region_expected_shared_kmers: f64,
    pub region_poisson_score: f64,
    /// Raw Poisson survival probability behind region_poisson_score (see
    /// MatchedRegion::tail_probability). Reported so downstream tools that need a probability
    /// (e.g. Benjamini-Hochberg correction) don't have to invert the -log10 transform.
    pub region_tail_probability: f64,
    pub region_enrichment: f64,
    /// Sum of IDF over the k-mers inside this region (see MatchedRegion::tfidf).
    pub region_tfidf: f64,
    /// region_tfidf divided by region_n_shared_kmers (see MatchedRegion::mean_idf).
    pub region_mean_idf: f64,
    /// Encoded positions inside the region where query and target disagree. Zero unless the
    /// search ran with `--extend-mismatch-penalty`.
    pub region_n_mismatches: u32,
    /// Karlin-Altschul bit score of the region (see MatchedRegion::ka_bits). 0 without
    /// `--extend-mismatch-penalty`.
    /// The region's own chance match rate u (see MatchedRegion::ka_u). 0 without
    /// `--extend-mismatch-penalty`.
    pub region_ka_u: f64,
    /// The region's Karlin-Altschul lambda, nats per unit of score (see
    /// MatchedRegion::ka_lambda). **0 means the region is not assessable**: its own
    /// composition matches at or past the boundary penalty / (1 + penalty), so
    /// `region_evalue` is infinity for want of a null, not because the score was poor.
    pub region_ka_lambda: f64,
    pub region_ka_bits: f64,
    /// E-value of the region against the searched database (see MatchedRegion::evalue).
    pub region_evalue: f64,
    /// Number of extended regions chained into this row (see MatchedRegion::n_chained).
    pub region_n_chained: u32,
}

impl SearchResultCsv {
    /// Create a CSV row from a SearchResult and a MatchedRegion
    /// WHY: Every CSV row must have matched region data. This method combines the SearchResult
    /// similarity metrics with a specific matched region to create one CSV row. Each SearchResult
    /// will produce multiple CSV rows (one per matched region), with all similarity metrics
    /// repeated for each region.
    pub fn from_result_and_region(
        result: &SearchResult,
        region: &MatchedRegion,
        remove_low_complexity: bool,
    ) -> Self {
        // See the matching debug_assert in ProteinSearcher::compare: a region shorter than
        // ksize should never exist.
        debug_assert!(
            region.length >= result.ksize,
            "region shorter than ksize: length={}, ksize={}",
            region.length,
            result.ksize
        );
        Self {
            query_name: result.query_name.clone(),
            query_md5: result.query_md5.clone(),
            target_name: result.target_name.clone(),
            target_md5: result.target_md5.clone(),
            containment: result.containment,
            n_intersecting_hashes: result.n_intersecting_hashes,
            ksize: result.ksize,
            scaled: result.scaled,
            moltype: result.moltype.clone(),
            remove_low_complexity,
            jaccard: result.jaccard,
            max_containment: result.max_containment,
            average_abund: result.average_abund,
            median_abund: result.median_abund,
            std_abund: result.std_abund,
            containment_target_in_query: result.containment_target_in_query,
            f_weighted_target_in_query: result.f_weighted_target_in_query,
            query_tfidf: result.query_tfidf,
            mean_matched_kmer_freq: result.mean_matched_kmer_freq,
            sum_matched_kmer_freq: result.sum_matched_kmer_freq,
            query_expected_shared_kmers: result.query_expected_shared_kmers,
            query_enrichment: result.query_enrichment,
            joint_kmer_freq: result.joint_kmer_freq,
            query_poisson_pvalue: result.query_poisson_pvalue,
            region_search_space: result.region_search_space,
            db_n_targets: result.db_n_targets,
            db_n_kmers: result.db_n_kmers,
            run_n_queries: result.run_n_queries,
            region_start: region.start,
            region_end: region.end,
            region_subseq: region.subseq.clone(),
            target_start: region.target_start,
            target_end: region.target_end,
            target_subseq: region.target_subseq.clone(),
            moltype_seq: region.moltype_seq.clone(),
            region_length: region.length,
            region_n_shared_kmers: region.n_shared,
            region_expected_shared_kmers: region.expected_shared_kmers,
            region_poisson_score: region.poisson_score,
            region_tail_probability: region.tail_probability,
            region_enrichment: region.enrichment,
            region_tfidf: region.tfidf,
            region_mean_idf: region.mean_idf,
            region_n_mismatches: region.n_mismatches,
            region_ka_u: region.ka_u,
            region_ka_lambda: region.ka_lambda,
            region_ka_bits: region.ka_bits,
            region_evalue: region.evalue,
            region_n_chained: region.n_chained,
        }
    }
}

/// Search result for a single query-target pair
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResult {
    /// Query sequence name
    pub query_name: String,

    /// Query sequence MD5 hash
    pub query_md5: String,

    /// Target sequence name
    pub target_name: String,

    /// Target sequence MD5 hash
    pub target_md5: String,

    /// Containment score (intersection / query_size)
    pub containment: f64,

    /// Number of intersecting k-mers
    pub n_intersecting_hashes: usize,

    /// K-mer size used
    pub ksize: u32,

    /// Scaled factor used
    pub scaled: u32,

    /// Molecular type (hp, dayhoff, protein)
    pub moltype: String,

    /// Jaccard similarity (intersection / union)
    pub jaccard: f64,

    /// Maximum containment (max of query->target and target->query containment)
    pub max_containment: f64,

    /// Average abundance of intersecting k-mers
    pub average_abund: f64,

    /// Median abundance of intersecting k-mers
    pub median_abund: f64,

    /// Standard deviation of abundance of intersecting k-mers
    pub std_abund: f64,

    /// Containment of target in query
    pub containment_target_in_query: f64,

    /// Weighted fraction of target in query
    pub f_weighted_target_in_query: f64,

    /// TF-IDF score for the query signature against the target database
    pub query_tfidf: f64,

    /// Mean frequency of matched k-mers in the target database: mean(freq_target[h]/N) over intersection.
    /// Higher = matched k-mers are common in the target DB (less discriminative).
    pub mean_matched_kmer_freq: f64,

    /// Sum of target-DB frequencies for matched k-mers: Σ freq_target[h]/N over intersection.
    /// Higher = matched k-mers are collectively more common in the target DB.
    pub sum_matched_kmer_freq: f64,

    /// Expected number of shared k-mers by chance: Σ freq_target[h]/N over ALL query k-mers.
    /// Uses this query's specific k-mers, so it is query-dependent.
    pub query_expected_shared_kmers: f64,

    /// Fold-enrichment: n_intersecting_hashes / query_expected_shared_kmers.
    /// Higher = more k-mers matched than expected by chance given target DB composition.
    pub query_enrichment: f64,

    /// Joint k-mer frequency: Σ freq_query[h] * freq_target[h] for h in intersection.
    /// Dot product of query and target frequency vectors over shared k-mers (two-pass;
    /// 0.0 if query frequencies not provided via set_query_frequencies()).
    pub joint_kmer_freq: f64,

    /// Poisson p-value: P(X ≥ n_intersecting_hashes | λ = query_expected_shared_kmers).
    /// 1.0 when query_expected_shared_kmers is unavailable (no database context).
    pub query_poisson_pvalue: f64,

    /// Number of positions in this query where a region could have started:
    /// `query_length - ksize + 1`. Regions are chosen after the fact (the best-looking
    /// gapless run of shared k-mers is kept), so a p-value computed on one region looks better
    /// than it should unless it is corrected by how many candidate positions it was chosen
    /// from. This is that count. 0 when raw sequences are not stored.
    pub region_search_space: usize,

    /// Number of target signatures searched.
    pub db_n_targets: usize,

    /// Total k-mer occurrences across the database: the sum, over every distinct k-mer hash,
    /// of how many signatures contain that hash. Used as a stand-in for total residue count,
    /// which the index does not currently store.
    pub db_n_kmers: usize,

    /// Number of queries in this search run. Reported for anyone who wants to correct for
    /// having tested many queries in one run (a family-wise error rate correction). Not
    /// folded into any p-value here: a hit's reported significance must not change depending
    /// on what other queries happened to run alongside it in the same invocation.
    pub run_n_queries: usize,

    /// 1 or more regions of 1+ k-mers overlapping between query and target
    pub matched_regions: Vec<MatchedRegion>,
}

/// A region of k-mer overlap between the query and target sequences.
/// We use u16 (up to 65,535) for the integer indexing as the largest protein as of
/// Nov 2025 is PKZILLA-1 which is 45,212 amino acids long
/// Source: https://en.wikipedia.org/wiki/Prymnesin-1
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatchedRegion {
    /// Query sequence name
    pub query_name: String,

    /// Start position of this region in the query sequence
    pub start: u32,

    /// End position of this region in the query sequence
    pub end: u32,

    /// Query subsequence (stitched k-mers)
    pub subseq: String,

    /// Target sequence name
    pub target_name: String,

    /// Target sequence start position
    pub target_start: u32,

    /// Target sequence end position
    pub target_end: u32,

    /// Target subsequence (stitched k-mers)
    pub target_subseq: String,

    /// Encoded sequence (hp/dayhoff/protein encoding)
    pub moltype_seq: String,

    // One of "protein20", "dayhoff6", or "hp"
    pub moltype: MolType,

    /// Length of the match
    pub length: u32,

    /// Exact shared k-mers inside the region that are in the sketch. For an exact region at
    /// scaled=1 this is `length - ksize + 1`. Two things make it smaller: at scaled>1 only
    /// sampled k-mers are in the sketch, and extension (see `extend_regions`) adds residues
    /// that were not shared k-mers. Either way the region spans the whole match and only the
    /// shared, sampled k-mers count as observations. The Poisson test compares this against
    /// `expected_shared_kmers`, which is summed over the same k-mers, so the two stay on the
    /// same footing.
    pub n_shared: u32,

    /// Encoded positions inside the region where query and target disagree. Zero unless the
    /// region was extended with a mismatch penalty.
    pub n_mismatches: u32,

    /// Expected number of shared k-mers by chance within this region: for every query k-mer
    /// whose start position falls inside this region, sum how often that k-mer's hash appears
    /// across the database, divided by the number of signatures in the database. 0.0 without
    /// DB context.
    ///
    /// The frequencies used are averaged over the whole database, not this region's own local
    /// composition. A region sitting in an unusual stretch of the protein (for example, an
    /// unusually hydrophobic stretch) looks more enriched than it should, because it is
    /// compared against the database average rather than against similar local sequence.
    pub expected_shared_kmers: f64,

    /// -log10 of the Poisson survival-function probability for this region: the probability
    /// of seeing at least `n_shared` shared k-mers if matches happened at random, given the
    /// rate `expected_shared_kmers`, then negative-log10-transformed so bigger means more
    /// surprising. A probability of 1.0 (no evidence at all) is a score of 0.0; a probability
    /// of 0.05 is a score of about 1.3; a probability of 0.0007 is a score of about 3.16.
    /// Uses the region's own k-mer count instead of the whole protein's. 0.0 without DB
    /// context.
    ///
    /// Named `poisson_score`, not `poisson_pvalue`, on purpose: treat it as a heuristic score
    /// for ranking candidate regions against each other, not as a calibrated probability. Two
    /// problems keep the underlying probability from being a real p-value, and the -log10
    /// transform changes neither of them:
    ///
    /// 1. `n_shared` is not an independent observation. It is `region_length - ksize + 1`,
    ///    arithmetic on the region's own length, and the region's length is exactly what
    ///    `find_matched_regions` chose by keeping the longest gapless run of shared k-mers. The
    ///    test is being applied to the same quantity that defined the region, which is close to
    ///    circular.
    ///
    /// 2. The k-mers being counted overlap by `ksize - 1` residues, so they are not independent
    ///    trials the way the Poisson model assumes. Five overlapping k-mers spanning a single
    ///    19-residue stretch are closer to one piece of evidence, observed five times, than to
    ///    five separate pieces of evidence. Treating them as independent understates how likely
    ///    a run this long is to appear by chance.
    ///
    /// To account for having picked the best-looking window out of many candidates, convert
    /// back to a probability and multiply by `region_search_space` (how many positions a
    /// region could have started at) and `db_n_targets` (how many targets were searched):
    /// `10f64.powf(-poisson_score) * region_search_space * db_n_targets`. That correction still
    /// does not fix either problem above.
    ///
    /// A properly calibrated version of this statistic would model the length of the longest
    /// gapless run directly (an extreme-value distribution, the same kind of model behind
    /// BLAST's E-values), account for the k-mer overlap, and be checked empirically against
    /// a decoy database. None of that is implemented here. Use this score to prioritize which
    /// regions to look at first, not to make a significance claim about any single region.
    pub poisson_score: f64,

    /// The raw Poisson survival-function probability `poisson_score` was computed from, before
    /// the -log10 transform: `10f64.powf(-poisson_score)`. Reported alongside `poisson_score`
    /// so downstream code that needs a probability (for example Benjamini-Hochberg FDR
    /// correction, which multiplies and ranks p-values directly) doesn't have to reconstruct
    /// one by undoing the log. Not named `poisson_pvalue`: it carries the same two structural
    /// problems documented on `poisson_score` and is not a calibrated p-value either. 1.0 (no
    /// evidence) without DB context.
    pub tail_probability: f64,

    /// Fold-enrichment scoped to this region: n_shared / expected_shared_kmers. 0.0 without DB
    /// context or when expected_shared_kmers is 0.
    pub enrichment: f64,

    /// TF-IDF scoped to this region: the sum of `ln(N / freq_target(h))` over every query
    /// k-mer whose start position falls inside the region, with term frequency fixed at 1,
    /// the same weighting `SearchResult::query_tfidf` applies to the whole query. Every
    /// k-mer in a region is shared with the target by construction, so this is the summed
    /// rarity of the k-mers that make up the match. 0.0 without DB context.
    ///
    /// This is not independent evidence from `expected_shared_kmers`: both are built from the
    /// same per-position `freq_target(h) / N`, one summing it linearly and the other summing
    /// its negative log. It also grows with region length the same way `n_shared` does, so
    /// it inherits the length circularity described on `poisson_score`.
    pub tfidf: f64,

    /// `tfidf / n_shared`: the average rarity of one k-mer in this region, so a short run of
    /// rare k-mers and a long run of common ones can be told apart without the length term.
    /// Divides by the same `n_shared` the Poisson test uses. 0.0 without DB context.
    pub mean_idf: f64,
    /// The chance match rate `u` of this region: how often a query position and a target
    /// position drawn from the two spans this region covers carry the same class. Counted
    /// on the spans alone, not on the whole proteins (see `region_match_probability`), so
    /// a polar-rich stretch inside an ordinary protein is measured as the polar-rich
    /// stretch it is. For a chain, counted over the chained span. 0.0 without extension.
    pub ka_u: f64,

    /// This region's Karlin-Altschul lambda, in nats per unit of score, with
    /// `ExtensionParams::ka_lambda_scale` already applied: the positive root of
    /// u e^x + (1-u) e^(-penalty x) = 1 at this region's own `ka_u`
    /// (Karlin & Altschul 1990; composition-based after Schaffer et al. 2001).
    ///
    /// **0.0 means the region is not assessable**, not that it scored badly. There is no
    /// positive root once `ka_u` reaches penalty / (1 + penalty), which is what two
    /// hydrophobic runs or two low-complexity stretches look like: matching is what those
    /// two compositions do by chance, so no length of agreement is evidence and no
    /// E-value exists. `ka_bits` is then 0.0 and `evalue` is infinity, the same values
    /// they take with no extension at all; this field is what tells the two apart.
    pub ka_lambda: f64,

    /// Karlin-Altschul bit score of the region as an ungapped alignment in the encoded
    /// alphabet: (`ka_lambda` * S - ln K) / ln 2, where S = matches - penalty * mismatches
    /// over the region. 0.0 when the region is not assessable (`ka_lambda` is 0.0), which
    /// is the property that makes this the ranking statistic for extended regions rather
    /// than the Poisson count, which sees a transmembrane helix against any other as a
    /// long exact run. Requires an extension penalty (`ExtensionParams`), since the
    /// score's mismatch term is the penalty; 0.0 otherwise.
    pub ka_bits: f64,

    /// How many extended regions this row is a chain of (see `chain_regions`). 1 for a region
    /// that stands alone, which is every region unless `--chain-max-gap` is set.
    pub n_chained: u32,

    /// E-value for `ka_bits` against the searched database:
    /// K * m * n * exp(-`ka_lambda` * S), with m the query length and n the database's
    /// residue count (`db_n_kmers` stands in for it). K is `ExtensionParams::ka_k`, which
    /// has to be calibrated on decoys for the alphabet and penalty in use. Infinity
    /// without extension or DB context, and infinity when the region is not assessable -
    /// read `ka_lambda` to tell those apart.
    pub evalue: f64,
}

/// P(X >= observed | lambda) via the Poisson survival function, 1 - CDF(observed - 1).
///
/// Returns 1.0 (no evidence) rather than erroring when there is no usable null: a lambda of
/// zero or an observation of zero leaves nothing to be surprised by, and `Poisson::new` rejects
/// non-positive or non-finite rates.
fn poisson_survival(observed: u32, lambda: f64) -> f64 {
    if observed == 0 || lambda <= 0.0 {
        return 1.0;
    }
    match Poisson::new(lambda) {
        Ok(dist) => (1.0 - dist.cdf((observed - 1) as u64)).max(0.0),
        Err(_) => 1.0,
    }
}

/// -log10 of a probability in (0, 1], so bigger means more surprising. 0.0 at pvalue = 1.0
/// (no evidence). Floors the input at `f64::MIN_POSITIVE` before taking the log, so a
/// probability that underflows to exactly 0.0 (an extreme tail value past f64 precision)
/// produces a large finite score instead of infinity.
fn neg_log10_score(pvalue: f64) -> f64 {
    -pvalue.max(f64::MIN_POSITIVE).log10()
}

/// Observed over expected. 0.0 when there is no expectation to divide by, which reads as
/// "not computable" rather than the +inf the division would produce.
fn fold_enrichment(observed: u32, expected: f64) -> f64 {
    if expected > 0.0 {
        observed as f64 / expected
    } else {
        0.0
    }
}

/// Expected shared k-mers by chance within one region:
///
/// ```text
/// lambda = sum, over every query k-mer whose start position falls inside the region,
///          of freq_target(h) / N
/// ```
///
/// `freq_target(h)` is how many target signatures contain a k-mer with hash `h`, and `N` is the
/// total number of target signatures. Note there is no division by the region's length here:
/// each term in the sum is already a per-k-mer database frequency, and the sum has one term per
/// k-mer position inside the region. This is the same formula as
/// `ProteinSearcher::calculate_expected_shared_kmers`, restricted to k-mer positions inside the
/// region instead of the whole query.
///
/// The window is `[start, end - ksize + 1)`, not the region's full span: a k-mer belongs to the
/// region only if it fits entirely inside. That is what makes the count equal
/// `length - ksize + 1` at scaled=1.
///
/// `prefix` is `PreparedQuery::position_prefix`: `freq_target(h)/N` already summed by position,
/// one entry per query, built once regardless of how many targets or regions it is looked up
/// for. Handing it `PreparedQuery::idf_prefix` instead gives the region's TF-IDF over the
/// same window. This turns the lookup into a difference of two prefix sums, O(1), instead of
/// the O(query k-mer count) rescan that computing lambda from scratch for every region on
/// every target would otherwise cost.
fn region_expectation(prefix: &[f64], start: u32, end: u32, ksize: usize) -> f64 {
    let last_index = prefix.len() - 1;
    let window_start = (start as usize).min(last_index);
    // A k-mer at p covers [p, p + ksize), so it fits inside [start, end) only when
    // p <= end - ksize. A span shorter than ksize holds no whole k-mer at all - saturating
    // here would wrongly admit position 0.
    let window_end = match (end as usize).checked_sub(ksize) {
        Some(last_start) => (last_start + 1).min(last_index),
        None => window_start,
    };
    prefix[window_end] - prefix[window_start]
}

impl Display for MatchedRegion {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Query Name: {}", self.query_name)?;
        writeln!(f, "Match Name: {}", self.target_name)?;
        writeln!(f, "query: {} ({}-{})", self.subseq, self.start, self.end)?;
        writeln!(f, "alpha: {}", self.moltype_seq)?;
        write!(f, "match: {} ({}-{})", self.target_subseq, self.target_start, self.target_end)
    }
}

/// Search statistics for TF-IDF and probability calculations
#[derive(Debug, Clone)]
pub struct SearchStats {
    /// Total number of signatures in the database
    pub total_signatures: usize,
    /// IDF values for each k-mer hash
    pub idf: HashMap<u64, f64>,
    /// Frequency of each k-mer hash across all signatures
    pub kmer_frequencies: HashMap<u64, usize>,
}

/// Pre-computed query data for efficient batch searching
///
/// WHY: This struct encapsulates all the query-specific data that needs to be computed once
/// per query and reused across many target comparisons. Instead of passing 7 separate parameters
/// to `compare()`, we bundle them into a single struct. This makes the API cleaner and reduces
/// the chance of errors from passing incorrect parameters. The struct is designed for performance:
/// it pre-computes expensive operations (like extracting mins as a HashSet and calculating TF-IDF)
/// so they're only done once per query, not once per query-target pair.
pub struct PreparedQuery<'a> {
    /// Reference to the query sketch
    pub sketch: &'a ProteinSketch,
    /// Pre-computed minhash values as a HashSet for efficient intersection calculations
    pub mins: HashSet<u64>,
    /// Pre-computed TF-IDF score for the query
    pub tfidf: f64,
    /// Prefix sums of target-DB k-mer frequency by query position, indexed 0..=max_position.
    /// `position_prefix[p]` is the sum of freq_target(h)/N over every k-mer whose start
    /// position is < p. Lets `region_expectation` answer any region window in O(1) instead of
    /// rescanning every one of the query's k-mers per region per target (see
    /// `ProteinSearcher::build_position_prefix`).
    pub position_prefix: Vec<f64>,
    /// Same layout as `position_prefix`, but summing IDF (`ln(N / freq_target(h))`) instead of
    /// `freq_target(h)/N`, so a region's TF-IDF is the same O(1) prefix difference (see
    /// `ProteinSearcher::build_idf_prefix`).
    pub idf_prefix: Vec<f64>,
}

impl SearchStats {
    /// Search statistics from the structures `ProteomeIndex::load_search_cache` returns.
    fn from_cache(total_signatures: usize, kmer_frequencies: HashMap<u64, usize>) -> Self {
        let idf: HashMap<u64, f64> = kmer_frequencies
            .iter()
            .map(|(&kmer, &freq)| (kmer, (total_signatures as f64 / freq as f64).ln()))
            .collect();
        Self { total_signatures, idf, kmer_frequencies }
    }
}

/// Protein signature searcher
pub struct ProteinSearcher {
    index: ProteomeIndex,
    stats: SearchStats,
    /// Ordered list of target MD5 keys, indexed by u32 for the inverted index
    target_list: Vec<String>,
    /// Inverted k-mer index: kmer_hash → Vec of target indices into target_list.
    /// Enables candidate pre-filtering so search_one only compares against targets
    /// that share at least one k-mer with the query, instead of all targets.
    inverted_index: HashMap<u64, Vec<u32>>,
    /// Lazy signature cache: populated on demand from RocksDB during search_one().
    ///
    /// WHY: When on-demand loading is used (fast startup path), the same target protein
    /// may be a candidate for many queries. Without a cache, each query would trigger a
    /// separate RocksDB get() for the same signature. The DashMap cache is thread-safe,
    /// allowing concurrent reads from rayon's par_iter() in search_one(). Memory grows
    /// only as candidates are accessed — hot targets are cached, cold ones are never loaded.
    sig_cache: DashMap<String, ProteinSketch>,
    /// Other entry names stored under each target md5 (see `ProteomeIndex::aliases`). A hit
    /// on the md5 is reported once more under each of them.
    aliases: HashMap<String, Vec<String>>,
    /// K-mer frequencies across the query proteome, set via set_query_frequencies() before
    /// searching. None = single-pass mode; joint_kmer_freq will be 0.0 for all results.
    query_kmer_frequencies: Option<HashMap<u64, usize>>,
    /// Total number of query sequences used to build query_kmer_frequencies.
    total_queries: usize,
    /// Total k-mer occurrences across the database: the sum, over every distinct k-mer hash,
    /// of how many signatures contain that hash. Computed once here and reused for every
    /// result, rather than per comparison, since summing it fresh would cost O(unique k-mers)
    /// on every candidate pair.
    db_n_kmers: usize,
    /// Seed extension, set via `set_extension()`. None keeps every region an exact run.
    extension: Option<ExtensionParams>,
}

impl ProteinSearcher {
    /// Extend every matched region past its exact seed with these parameters. Off by default.
    pub fn set_extension(&mut self, params: Option<ExtensionParams>) {
        self.extension = params;
    }

    /// Fit lambda and K of this index for one mismatch penalty and X-drop by searching
    /// `n_queries` calibration queries built from the index's own sequences (`DecoyNull`)
    /// and reading the slope and intercept off ln(regions with score S) against S
    /// (`fit_scores`), the homolog excess cut off.
    ///
    /// Filters are wide open (one shared k-mer, no p-value cap) so the curve holds every
    /// region the procedure can produce; a filter a search adds only removes regions and
    /// makes its E-values conservative. A query's hits on its own database entry are left
    /// out. The searcher's extension setting is restored afterwards.
    pub fn calibrate_ka(
        &mut self,
        settings: KaCalibrationSettings,
    ) -> IndexResult<KaCalibrationReport> {
        let previous = self.extension;
        self.extension = Some(ExtensionParams {
            mismatch_penalty: settings.mismatch_penalty,
            xdrop: settings.xdrop,
            ka_k: 1.0,
            ka_lambda_scale: 1.0,
            chain_max_gap: 0,
            chain_max_shift: 0,
        });
        let report = self.run_calibration(settings);
        self.extension = previous;
        report
    }

    fn run_calibration(&self, settings: KaCalibrationSettings) -> IndexResult<KaCalibrationReport> {
        let KaCalibrationSettings {
            mismatch_penalty,
            xdrop,
            null,
            reference,
            reference_shuffles,
            n_queries,
            seed,
        } = settings;
        let (queries, match_probability) = self.calibration_queries(null, n_queries, 1, seed)?;
        let scores = self.calibration_scores(&queries);
        // The database null is censored against the same queries shuffled (`reference`),
        // each of them `reference_shuffles` times.
        let (fit, reference_survival, n_reference_queries, n_reference_regions) =
            if null == DecoyNull::Database {
                let (shuffled, _) =
                    self.calibration_queries(reference, n_queries, reference_shuffles, seed)?;
                let reference_scores = self.calibration_scores(&shuffled);
                (
                    fit_scores_with_reference(&scores, &reference_scores),
                    survival_counts(&reference_scores),
                    shuffled.len(),
                    reference_scores.len(),
                )
            } else {
                (fit_scores(&scores), Vec::new(), 0, 0)
            };
        let survival = survival_counts(&scores);
        let n_queries = queries.len();
        let query_residues: u64 =
            queries.iter().map(|(_, q)| q.get_raw_sequence().map_or(0, |r| r.len() as u64)).sum();
        let database_kmers = self.db_n_kmers as u64;
        let lambda_analytic = karlin_altschul_lambda(match_probability, mismatch_penalty);
        let fitted = fit.map(|fit| KaCalibration {
            mismatch_penalty,
            xdrop,
            null,
            seed,
            n_queries,
            query_residues,
            database_kmers,
            n_regions: scores.len(),
            match_probability,
            lambda_analytic,
            // The fit ran on x / BIN_WIDTH, so its slope per bin is slope-per-nat x width.
            slope: fit.lambda / BIN_WIDTH,
            // ln(regions in the bin at x) = ln(K L N (1 - e^(-slope w))) - slope x, and
            // slope w is the fit's slope per bin.
            k: fit.ln_intercept.exp()
                / (query_residues as f64 * database_kmers as f64 * (1.0 - (-fit.lambda).exp())),
            bin_width: BIN_WIDTH,
            score_lo: fit.score_lo,
            score_hi: fit.score_hi,
            bend_score: fit.bend_score,
            rms_residual: fit.rms_residual,
            survival: fit.survival,
            reference_survival: fit.reference_survival,
            reference_lambda: fit.reference_lambda,
            reference: (null == DecoyNull::Database).then_some(reference),
        });
        Ok(KaCalibrationReport {
            fitted,
            match_probability,
            n_queries,
            n_regions: scores.len(),
            n_reference_queries,
            n_reference_regions,
            query_residues,
            database_kmers,
            lambda_analytic,
            survival,
            reference_survival,
        })
    }

    /// The calibration queries, each with the md5 of the database entry it came from, and
    /// the database's match probability read off at least `COMPOSITION_SAMPLE` targets.
    ///
    /// `target_list` comes out of a DashMap, whose order changes from run to run, so the
    /// sample is drawn from the md5s in sorted order to make the fit reproducible.
    /// `repeats` decoys are made from each sequence picked, not one, so the reference
    /// curve can be filled in without searching a different part of the database than the
    /// real curve did: the same sequences, shuffled `repeats` times over.
    fn calibration_queries(
        &self,
        null: DecoyNull,
        n_queries: usize,
        repeats: usize,
        seed: u64,
    ) -> IndexResult<(Vec<(String, ProteinSketch)>, f64)> {
        let mut md5s: Vec<&String> = self.target_list.iter().collect();
        md5s.sort_unstable();
        let mut rng = SplitMix64::new(seed);
        let picks = rng.sample_indices(md5s.len(), n_queries.max(COMPOSITION_SAMPLE));
        let wanted = n_queries.saturating_mul(repeats.max(1));
        let mut queries = Vec::with_capacity(wanted);
        let mut class_counts: HashMap<u8, f64> = HashMap::new();
        for idx in picks {
            let md5 = md5s[idx].clone();
            let Some(target) = self.target_sketch(&md5)? else { continue };
            if let Some(encoded) = target.get_class_sequence() {
                for b in encoded.bytes() {
                    *class_counts.entry(b).or_insert(0.0) += 1.0;
                }
            }
            for _ in 0..repeats.max(1) {
                if queries.len() >= wanted {
                    break;
                }
                if let Some(query) = self.decoy_from_target(&target, null, &mut rng)? {
                    queries.push((md5.clone(), query));
                }
            }
        }
        let total: f64 = class_counts.values().sum::<f64>().max(1.0);
        let match_probability = class_counts.values().map(|c| (c / total).powi(2)).sum();
        Ok((queries, match_probability))
    }

    /// Normalised score x = lambda_region S of every region the calibration queries produce,
    /// in units of `BIN_WIDTH`, a query's own database entry excluded. The searcher runs
    /// with K = 1 and lambda scale 1 during calibration, so a region's `ka_bits` x ln 2 is
    /// exactly lambda S with the closed-form lambda of that region's own two spans.
    ///
    /// A region whose spans sit past the composition boundary has lambda 0, so `ka_bits`
    /// is 0 and it lands in `karlin_altschul::DEGENERATE_BIN`, which `tail_bins` already
    /// bars from being the peak. Those regions are counted, not fitted.
    fn calibration_scores(&self, queries: &[(String, ProteinSketch)]) -> Vec<f64> {
        let filters = SearchFilters {
            threshold: 0.0,
            min_shared_kmers: 1,
            max_query_pvalue: f64::INFINITY,
            min_region_score: f64::NEG_INFINITY,
            // A query's hits on its own database entry are dropped by md5 below, which
            // also covers the aliases #63 reports under the same md5.
            skip_self_matches: false,
        };
        queries
            .par_iter()
            .flat_map_iter(|(source_md5, query)| {
                self.search_one(query, &filters, queries.len())
                    .into_iter()
                    .filter(|r| &r.target_md5 != source_md5)
                    .flat_map(|r| {
                        r.matched_regions
                            .into_iter()
                            .map(|region| region.ka_bits * std::f64::consts::LN_2 / BIN_WIDTH)
                    })
                    .collect::<Vec<f64>>()
            })
            .collect()
    }

    /// The calibration query for one target, built the way a search query would be, or
    /// None when the index did not store raw sequences.
    fn decoy_from_target(
        &self,
        target: &ProteinSketch,
        null: DecoyNull,
        rng: &mut SplitMix64,
    ) -> IndexResult<Option<ProteinSketch>> {
        let Some(raw) = target.get_raw_sequence() else { return Ok(None) };
        let decoy_seq = make_decoy(raw, null, rng);
        let mut decoy = ProteinSketch::new(
            &format!("{null}_{}", target.signature().name),
            self.index.ksize(),
            self.index.scaled(),
            self.index.moltype(),
        )?;
        decoy.set_remove_low_complexity(self.index.remove_low_complexity());
        decoy.add_protein(&decoy_seq, true)?;
        Ok(Some(decoy))
    }

    /// A target by md5 from memory, the search cache, or RocksDB, cached for later searches.
    fn target_sketch(&self, md5: &str) -> IndexResult<Option<ProteinSketch>> {
        if let Some(entry) = self.index.get_signatures().get(md5) {
            return Ok(Some(entry.value().clone()));
        }
        if let Some(entry) = self.sig_cache.get(md5) {
            return Ok(Some(entry.value().clone()));
        }
        let Some(target) = self.index.get_signature_by_md5(md5)? else { return Ok(None) };
        self.sig_cache.insert(md5.to_string(), target.clone());
        Ok(Some(target))
    }

    /// The lambda scale and K a search should use for this penalty and X-drop, and where
    /// they came from, in order of preference: `explicit` (`--ka-k`, with the closed-form
    /// lambda); a fit stored in the index for this pair, whatever null it used; a fresh fit
    /// on `settings.n_queries` calibration queries if that is nonzero. Otherwise an error: an
    /// E-value without a fit for its own index is not printed.
    pub fn resolve_ka(
        &mut self,
        explicit: Option<f64>,
        settings: KaCalibrationSettings,
    ) -> IndexResult<(KaParams, KaSource)> {
        let KaCalibrationSettings { mismatch_penalty, xdrop, .. } = settings;
        if let Some(k) = explicit {
            return Ok((KaParams { k, lambda_scale: 1.0 }, KaSource::Flag));
        }
        if let Some(stored) = self.index.ka_calibration(mismatch_penalty, xdrop)? {
            let params = KaParams { k: stored.k, lambda_scale: stored.lambda_scale() };
            return Ok((params, KaSource::Index(stored)));
        }
        let report = self.calibrate_ka(settings)?;
        match report.fitted {
            Some(fit) => {
                let params = KaParams { k: fit.k, lambda_scale: fit.lambda_scale() };
                Ok((params, KaSource::Fitted(fit)))
            }
            None => Err(anyhow::anyhow!(
                "no Karlin-Altschul fit for penalty {mismatch_penalty}, X-drop {xdrop}: {} \
                 calibration queries gave {} regions, fewer than the {MIN_FIT_POINTS} score \
                 bins of 30 regions the fit needs. Rebuild the index with --ka-queries set \
                 higher, search with --ka-queries, or pass --ka-k.",
                report.n_queries,
                report.n_regions
            )
            .into()),
        }
    }

    /// Create a searcher over an index built in this process.
    ///
    /// Finalizes the index so its inverted index is on disk, then reads the search
    /// structures back exactly as `load` would for a saved index. Sketches the index
    /// holds in memory (from `store_signatures` or `load`) are used directly; the rest
    /// are read on demand.
    pub fn new(index: ProteomeIndex) -> IndexResult<Self> {
        index.finalize()?;
        let cache = index.load_search_cache()?.ok_or(IndexError::NoSavedState)?;
        Self::from_cache(index, cache)
    }

    /// Load a searcher from a saved index.
    ///
    /// Opens the database without loading any signatures into memory. Signatures are
    /// loaded on demand during search via `get_signature_by_md5()`.
    pub fn load<P: AsRef<Path>>(path: P) -> IndexResult<Self> {
        let index = ProteomeIndex::open_for_search(&path)?;
        let cache = index.load_search_cache()?.ok_or(IndexError::NoSavedState)?;
        let (targets, kmers) = (cache.target_list.len(), cache.inverted_index.len());
        eprintln!("Loaded search cache: {targets} targets, {kmers} k-mers indexed");
        Self::from_cache(index, cache)
    }

    fn from_cache(index: ProteomeIndex, cache: SearchCache) -> IndexResult<Self> {
        let stats = SearchStats::from_cache(cache.target_list.len(), cache.kmer_frequencies);
        let db_n_kmers = stats.kmer_frequencies.values().sum();
        let aliases = index.aliases()?;
        Ok(Self {
            index,
            stats,
            target_list: cache.target_list,
            inverted_index: cache.inverted_index,
            sig_cache: DashMap::new(),
            aliases,
            query_kmer_frequencies: None,
            total_queries: 0,
            db_n_kmers,
            extension: None,
        })
    }

    /// `result` again under every other entry name stored with its target's sketch. Those
    /// entries have the same k-mer set, so every statistic and region is the same; only the
    /// name differs.
    fn alias_results(&self, result: &SearchResult) -> Vec<SearchResult> {
        let Some(names) = self.aliases.get(&result.target_md5) else {
            return Vec::new();
        };
        names
            .iter()
            .map(|name| {
                let mut copy = result.clone();
                copy.target_name = name.clone();
                for region in &mut copy.matched_regions {
                    region.target_name = name.clone();
                }
                copy
            })
            .collect()
    }

    /// Prepare a query for efficient batch searching
    ///
    /// WHY: This method pre-computes expensive operations (extracting mins as HashSet, calculating
    /// TF-IDF) that would otherwise be repeated for each target comparison. By doing this once
    /// per query, we improve performance for batch searches (1 query vs many targets). This follows
    /// the idiomatic Rust pattern of preparing data once and reusing it, rather than recomputing
    /// it repeatedly.
    ///
    /// # Arguments
    /// * `query` - The query sketch to prepare
    ///
    /// # Returns
    /// A `PreparedQuery` struct containing pre-computed query data
    pub fn prepare_query<'a>(&self, query: &'a ProteinSketch) -> PreparedQuery<'a> {
        PreparedQuery {
            sketch: query,
            mins: query.mins_as_set(),
            tfidf: self.calculate_tfidf(query),
            position_prefix: self.build_position_prefix(query),
            idf_prefix: self.build_idf_prefix(query),
        }
    }

    /// Builds the per-position frequency prefix sums consumed by `region_expectation`. Runs
    /// once per query in `prepare_query`, not once per region per target: `search_one` calls
    /// `compare` once per candidate target sharing the query, and each `compare` call rescopes
    /// the Poisson test to every matched region, so without this the query's full k-mer set
    /// would be rescanned target-count x region-count times instead of once.
    fn build_position_prefix(&self, query: &ProteinSketch) -> Vec<f64> {
        let total_signatures = self.stats.total_signatures as f64;
        Self::build_prefix(query, |hashval| {
            self.stats.kmer_frequencies.get(&hashval).copied().unwrap_or(1) as f64
                / total_signatures
        })
    }

    /// Per-position IDF prefix sums, the region-scoped counterpart of `calculate_tfidf`. A
    /// hash the database has never seen contributes 0, matching `calculate_tfidf`, which
    /// skips such hashes. Same once-per-query reasoning as `build_position_prefix`.
    fn build_idf_prefix(&self, query: &ProteinSketch) -> Vec<f64> {
        Self::build_prefix(query, |hashval| self.stats.idf.get(&hashval).copied().unwrap_or(0.0))
    }

    /// Lays `per_hash(h)` out by k-mer start position and returns its prefix sums, so any
    /// window `[a, b)` of positions sums to `prefix[b] - prefix[a]`.
    ///
    /// A position holding an ambiguous residue (B, J or Z) is recorded under every reading's
    /// hash, so a position can carry several hashes. Their values are added, which keeps the
    /// whole-query total equal to the sums `calculate_expected_shared_kmers` and
    /// `calculate_tfidf` take over every query hash.
    ///
    /// Sized to the highest k-mer start position actually present in `query.kmer_positions()`,
    /// not to the raw sequence length, since raw sequences are only optionally stored.
    fn build_prefix(query: &ProteinSketch, per_hash: impl Fn(u64) -> f64) -> Vec<f64> {
        let n_positions = query
            .kmer_positions()
            .values()
            .flat_map(|positions| positions.iter().copied())
            .max()
            .map_or(0, |max_pos| max_pos + 1);

        let mut position_value = vec![0.0; n_positions];
        for (&hashval, positions) in query.kmer_positions() {
            let value = per_hash(hashval);
            for &p in positions {
                position_value[p] += value;
            }
        }

        let mut prefix = Vec::with_capacity(n_positions + 1);
        prefix.push(0.0);
        let mut running = 0.0;
        for value in position_value {
            running += value;
            prefix.push(running);
        }
        prefix
    }

    /// Comprehensive search method that calculates all metrics including TF-IDF and overlap probability
    ///
    /// This is the single, idiomatic search method that replaces search_single, search_multiple,
    /// and search_with_kmer_extraction. It performs parallel processing for multiple queries
    /// and calculates all similarity metrics in one pass for efficiency.
    ///
    /// # Arguments
    /// * `queries` - Slice of query signatures to search against the database
    /// * `filters` - Result-level filters applied while searching; failing results are never
    ///   allocated into the returned `Vec` (see `SearchFilters`)
    ///
    /// # Returns
    /// Vector of SearchResult containing all similarity metrics, sorted by containment score
    #[must_use = "search results should be used to process query matches"]
    pub fn search(
        &self,
        queries: &[ProteinSketch],
        filters: &SearchFilters,
    ) -> Result<Vec<SearchResult>> {
        // Create progress bar for tracking query processing
        let progress = ProgressBar::new(queries.len() as u64);
        progress.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} queries ({percent}%) | ETA: {eta}",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        progress.set_message("Searching...");

        let total_queries = queries.len();

        // Perform parallel search across all queries using the inverted index
        let all_results: Vec<SearchResult> = queries
            .par_iter()
            .flat_map(|query| {
                // search_one uses the inverted index to find candidates
                let results = self.search_one(query, filters, total_queries);

                // Update progress after processing each query
                progress.inc(1);

                results
            })
            .collect();

        progress.finish_with_message("Search complete");

        // Sort by containment score (descending) - this is the primary ranking metric
        let mut sorted_results = all_results;
        sorted_results.sort_by(|a, b| {
            b.containment.partial_cmp(&a.containment).unwrap_or(std::cmp::Ordering::Equal)
        });

        Ok(sorted_results)
    }

    /// Search a single query against all targets in the index
    ///
    /// WHY: This enables streaming search where queries are processed one at a time
    /// without accumulating all query signatures in memory. This is critical for
    /// large-scale searches (e.g., 159k human proteins) where loading all queries
    /// into memory would require 50+ GB.
    ///
    /// Uses the inverted k-mer index to find candidate targets (those sharing ≥1 k-mer
    /// with the query) before doing full pairwise comparison. This skips the vast majority
    /// of targets with no k-mer overlap, reducing work from O(all_targets) to O(candidates).
    ///
    /// `filters` is applied inline as each candidate is compared: results that fail are never
    /// pushed into the returned `Vec`, rather than being built up and filtered out afterward.
    pub fn search_one(
        &self,
        query: &ProteinSketch,
        filters: &SearchFilters,
        total_queries: usize,
    ) -> Vec<SearchResult> {
        let prepared = self.prepare_query(query);

        // Find candidate targets via inverted index: only compare against targets
        // that share at least one k-mer with the query
        let mut candidate_set: HashSet<u32> = HashSet::new();
        for &kmer in &prepared.mins {
            if let Some(indices) = self.inverted_index.get(&kmer) {
                for &idx in indices {
                    candidate_set.insert(idx);
                }
            }
        }

        // Compare only against candidates, in parallel.
        // Three paths for signature access (checked in order of cost):
        //   1. In-memory DashMap (populated in slow path for old DBs): zero-copy O(1) lookup
        //   2. sig_cache (DashMap populated on demand): avoids repeated RocksDB reads for hot targets
        //   3. On-demand RocksDB loading: first access per target; result stored in sig_cache
        let sigs = self.index.get_signatures();
        let mut results: Vec<SearchResult> = candidate_set
            .par_iter()
            .filter_map(|&idx| {
                let md5 = &self.target_list[idx as usize];

                // Path 1: in-memory signatures (slow path for old DBs without search cache)
                if let Some(entry) = sigs.get(md5.as_str()) {
                    return self.compare(&prepared, entry.value(), filters, total_queries);
                }

                // Path 2: sig_cache hit (target was loaded by an earlier query)
                if let Some(entry) = self.sig_cache.get(md5.as_str()) {
                    return self.compare(&prepared, entry.value(), filters, total_queries);
                }

                // Path 3: first-time load from RocksDB; store in cache for future queries
                let target = self.index.get_signature_by_md5(md5).ok()??;
                let result = self.compare(&prepared, &target, filters, total_queries);
                self.sig_cache.insert(md5.clone(), target);
                result
            })
            .collect();
        if self.aliases.is_empty() {
            return results;
        }
        let extra: Vec<SearchResult> = results.iter().flat_map(|r| self.alias_results(r)).collect();
        results.extend(extra);
        results
    }

    /// Perform all-vs-all search without cloning signatures
    ///
    /// WHY: This method is optimized for all-vs-all searches where query and target are the same
    /// database. Instead of cloning all signatures (which is expensive for large databases), this
    /// method works directly with references from the index. We collect only the MD5 keys (cheap
    /// String clones) and parallelize over those, looking up the actual signatures in the DashMap.
    /// This avoids cloning the large ProteinSketch objects while still enabling parallel processing.
    /// The method also automatically skips self-matches by comparing MD5 sums in compare().
    ///
    /// `filters` is forwarded to `search_one`, so failing results are never allocated into the
    /// returned `Vec` (see `SearchFilters`).
    ///
    /// # Returns
    /// Vector of SearchResult containing all similarity metrics, sorted by containment score
    #[must_use = "search results should be used to process query matches"]
    pub fn search_all_vs_all(&self, filters: &SearchFilters) -> Result<Vec<SearchResult>> {
        // Every query is a target here, so its hit on itself is dropped.
        let filters = &SearchFilters { skip_self_matches: true, ..*filters };
        // Collect all signatures as owned ProteinSketch values to use as queries.
        // Two paths: in-memory DashMap (slow path / old DBs) or on-demand RocksDB (fast path).
        let all_queries: Vec<ProteinSketch> = {
            let sigs = self.index.get_signatures();
            if !sigs.is_empty() {
                sigs.iter().map(|entry| entry.value().clone()).collect()
            } else {
                // Fast path: load all signatures from RocksDB via target_list
                self.target_list
                    .iter()
                    .filter_map(|md5| {
                        // Check sig_cache first, then RocksDB
                        if let Some(entry) = self.sig_cache.get(md5.as_str()) {
                            return Some(entry.value().clone());
                        }
                        let sig = self.index.get_signature_by_md5(md5).ok()??;
                        self.sig_cache.insert(md5.clone(), sig.clone());
                        Some(sig)
                    })
                    .collect()
            }
        };

        // Create progress bar for tracking query processing
        let progress = ProgressBar::new(all_queries.len() as u64);
        progress.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} queries ({percent}%) | ETA: {eta}",
            )
            .unwrap()
            .progress_chars("=>-"),
        );
        progress.set_message("Searching all-vs-all...");

        let total_queries = all_queries.len();

        // Use search_one (inverted index) for each query instead of exhaustive O(N²) iteration
        let all_results: Vec<SearchResult> = all_queries
            .par_iter()
            .flat_map(|query| {
                let results = self.search_one(query, filters, total_queries);
                progress.inc(1);
                results
            })
            .collect();

        progress.finish_with_message("Search complete");

        // Sort by containment score (descending) - this is the primary ranking metric
        let mut sorted_results = all_results;
        sorted_results.sort_by(|a, b| {
            b.containment.partial_cmp(&a.containment).unwrap_or(std::cmp::Ordering::Equal)
        });

        Ok(sorted_results)
    }

    /// Calculate comprehensive similarity between query and target signatures including TF-IDF and overlap probability
    ///
    /// This method uses the standalone `calculate_similarity` function for the core calculation,
    /// then adds database-specific metrics (TF-IDF and overlap probability) that require the
    /// database context.
    ///
    /// # Why PreparedQuery?
    ///
    /// The `PreparedQuery` struct encapsulates all query-specific data that needs to be computed
    /// once per query and reused across many target comparisons. This avoids redundant allocations
    /// and calculations when comparing the same query against many targets. The standalone
    /// `calculate_similarity` function extracts this data itself, which is fine for 1v1 comparisons
    /// but less efficient for batch operations.
    ///
    /// # Performance Considerations
    ///
    /// For batch searches (1 query vs many targets), this method is more efficient than calling
    /// `calculate_similarity` directly because:
    /// 1. `query.mins` is extracted once and reused
    /// 2. `query.tfidf` is calculated once per query and reused
    /// 3. Database-specific metrics (overlap probability) are calculated using pre-computed stats
    ///
    /// For 1v1 comparisons or testing, use `calculate_similarity` directly.
    pub(crate) fn compare(
        &self,
        query: &PreparedQuery<'_>,
        target: &ProteinSketch,
        filters: &SearchFilters,
        total_queries: usize,
    ) -> Option<SearchResult> {
        // In an all-vs-all search every query is also a target, and its hit on itself says
        // nothing. When a FASTA is searched against an index, the query's identical entry is
        // a real hit (BCL-2 against a database holding BCL2_HUMAN), so the skip is opt-in.
        if filters.skip_self_matches && query.sketch.signature().md5sum == target.signature().md5sum
        {
            return None;
        }

        // Use pre-computed query.mins from PreparedQuery and compute target mins once.
        // WHY: query.mins is already a HashSet built in prepare_query(). Using it directly
        // avoids recomputing the query HashSet on every comparison. target_mins is computed
        // once here and reused for intersection and calculate_similarity_from_precomputed,
        // avoiding a second target HashSet allocation inside calculate_similarity.
        let target_mins = target.mins_as_set();
        let intersection: HashSet<u64> = query.mins.intersection(&target_mins).cloned().collect();

        // Skip if no intersection
        if intersection.is_empty() {
            return None;
        }

        // Containment and n_intersecting_hashes are cheap - just the intersection size - so
        // check them before the expensive per-pair work below (find_matched_regions walks both
        // sequences; abundance_stats sorts the intersection). Candidates failing these never
        // pay for that work.
        //
        // The query p-value and region score checks cannot be hoisted up here: a pair is kept
        // when either scope clears, and the region score isn't known until the regions exist.
        // So this filtering happens after the result is built, and pairs that fail the query
        // scope now pay for region-finding before being rejected.
        let n_intersecting_hashes = intersection.len();
        let containment = n_intersecting_hashes as f64 / query.mins.len() as f64;
        if containment < filters.threshold || n_intersecting_hashes < filters.min_shared_kmers {
            return None;
        }

        let query_expected_shared_kmers =
            self.calculate_expected_shared_kmers(query.sketch, target);
        let query_poisson_pvalue =
            poisson_survival(n_intersecting_hashes as u32, query_expected_shared_kmers);

        // Calculate database-specific overlap metrics
        let mean_matched_kmer_freq = self.calculate_mean_matched_kmer_freq(&intersection);
        let sum_matched_kmer_freq = self.calculate_sum_matched_kmer_freq(&intersection);

        // Build the similarity result using the pre-computed intersection and mins sets
        let mut result = calculate_similarity_from_precomputed(
            query.sketch,
            &query.mins,
            target,
            &target_mins,
            &intersection,
        )?;

        let query_enrichment =
            fold_enrichment(n_intersecting_hashes as u32, query_expected_shared_kmers);

        if let Some(params) = self.extension {
            result.matched_regions = extend_regions(
                std::mem::take(&mut result.matched_regions),
                query.sketch,
                target,
                params,
            );
        }

        // Rescope the same Poisson test to each matched region individually, so a tight local
        // match doesn't get diluted by the whole protein's k-mer count.
        let ksize = query.sketch.protein_ksize() as usize;
        for region in result.matched_regions.iter_mut() {
            // lambda: for each query k-mer positioned inside this region, how many target
            // signatures contain that k-mer's hash, divided by the total number of target
            // signatures, summed. See region_expectation for the exact formula.
            let lambda =
                region_expectation(&query.position_prefix, region.start, region.end, ksize);
            // find_matched_regions never emits a region shorter than ksize or with no
            // shared k-mer in it, and extension only grows a region. debug_assert catches
            // either loudly if that invariant is ever broken.
            debug_assert!(
                region.length >= ksize as u32 && region.n_shared >= 1,
                "bad region: length={}, n_shared={}, ksize={ksize}",
                region.length,
                region.n_shared
            );
            let n_shared = region.n_shared;
            let tail_probability = poisson_survival(n_shared, lambda);

            region.expected_shared_kmers = lambda;
            region.poisson_score = neg_log10_score(tail_probability);
            region.tail_probability = tail_probability;
            region.enrichment = fold_enrichment(n_shared, lambda);
            // Same window as lambda, summing IDF instead of frequency: how rare the k-mers
            // that make up this region are, in total and on average.
            region.tfidf = region_expectation(&query.idf_prefix, region.start, region.end, ksize);
            region.mean_idf = region.tfidf / n_shared as f64;
        }

        // Karlin-Altschul bits and E-value per region, each on the class composition of the
        // two spans that region covers rather than of the two whole proteins. Only
        // meaningful with a mismatch penalty, which is the score's mismatch term.
        if let Some(params) = self.extension {
            if let (Some(q_enc), Some(t_enc)) =
                (query.sketch.get_class_sequence(), target.get_class_sequence())
            {
                let (q_bytes, t_bytes) = (q_enc.as_bytes(), t_enc.as_bytes());
                let m = q_enc.len() as f64;
                let n = self.db_n_kmers as f64;
                for region in result.matched_regions.iter_mut() {
                    region.ka_u = region_match_probability(
                        region_span(q_bytes, region.start, region.end),
                        region_span(t_bytes, region.target_start, region.target_end),
                    );
                    region.ka_lambda = karlin_altschul_lambda(region.ka_u, params.mismatch_penalty)
                        * params.ka_lambda_scale;
                    let matches = region.length as f64 - region.n_mismatches as f64;
                    let raw = matches - params.mismatch_penalty * region.n_mismatches as f64;
                    // A lambda of 0 is "not assessable", not "scored badly": the region's own
                    // composition matches at or past penalty / (1 + penalty), so there is no
                    // null to be surprised against. Leave the no-evidence values in place.
                    if region.ka_lambda > 0.0 && params.ka_k > 0.0 {
                        region.ka_bits = ((region.ka_lambda * raw - params.ka_k.ln())
                            / std::f64::consts::LN_2)
                            .max(0.0);
                        region.evalue = params.ka_k * m * n * (-region.ka_lambda * raw).exp();
                    }
                }
                if params.chain_max_gap > 0 && params.ka_k > 0.0 {
                    result.matched_regions = chain_regions(
                        std::mem::take(&mut result.matched_regions),
                        q_bytes,
                        t_bytes,
                        query.sketch.get_raw_sequence(),
                        target.get_raw_sequence(),
                        params,
                        m,
                        t_enc.len() as f64,
                        self.stats.total_signatures as f64,
                    );
                }
            }
        }

        // Either scope clearing its cap keeps the pair - see SearchFilters::scopes_pass.
        // Bigger poisson_score is more surprising, so the best region is the highest-scoring
        // one.
        let best_region_score = result
            .matched_regions
            .iter()
            .map(|region| region.poisson_score)
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if !filters.scopes_pass(query_poisson_pvalue, best_region_score) {
            return None;
        }

        // Fill in database-specific metrics
        // joint_kmer_freq = Σ freq_query[h] * freq_target[h] for h in intersection.
        // Only computed in two-pass mode (when set_query_frequencies() has been called).
        let joint_kmer_freq = if let Some(qfreqs) = &self.query_kmer_frequencies {
            let total_q = self.total_queries as f64;
            let total_t = self.stats.total_signatures as f64;
            intersection
                .iter()
                .map(|&h| {
                    let fq = qfreqs.get(&h).copied().unwrap_or(0) as f64 / total_q;
                    let ft =
                        self.stats.kmer_frequencies.get(&h).copied().unwrap_or(0) as f64 / total_t;
                    fq * ft
                })
                .sum()
        } else {
            0.0
        };

        result.query_tfidf = query.tfidf;
        result.mean_matched_kmer_freq = mean_matched_kmer_freq;
        result.sum_matched_kmer_freq = sum_matched_kmer_freq;
        result.query_expected_shared_kmers = query_expected_shared_kmers;
        result.query_enrichment = query_enrichment;
        result.joint_kmer_freq = joint_kmer_freq;
        result.query_poisson_pvalue = query_poisson_pvalue;
        result.db_n_targets = self.stats.total_signatures;
        result.db_n_kmers = self.db_n_kmers;
        result.run_n_queries = total_queries;

        Some(result)
    }

    /// Set query-proteome k-mer frequencies for two-pass joint_kmer_freq computation.
    ///
    /// Call this after loading the searcher and before the search loop.
    /// `freqs` maps each k-mer hash to the number of query sequences containing it.
    /// `total` is the total number of query sequences processed.
    pub fn set_query_frequencies(&mut self, freqs: HashMap<u64, usize>, total: usize) {
        self.query_kmer_frequencies = Some(freqs);
        self.total_queries = total;
    }

    /// Calculate TF-IDF score for a query signature
    #[must_use = "TF-IDF score should be used to rank search results"]
    pub fn calculate_tfidf(&self, query: &ProteinSketch) -> f64 {
        let query_mins = query.signature().minhash.mins();
        let mut tfidf_sum = 0.0;

        for min in query_mins {
            if let Some(&idf) = self.stats.idf.get(&min) {
                // TF is the abundance of the k-mer in the query (simplified to 1 if no abundance tracking)
                let tf = 1.0; // Could be enhanced to use actual abundance if available
                tfidf_sum += tf * idf;
            }
        }

        tfidf_sum
    }

    /// Mean frequency of matched k-mers in the target database: mean(freq_target[h]/N) over intersection.
    /// Higher values = matched k-mers are more common in the target DB (less discriminative).
    /// Range: [0, 1].
    pub fn calculate_mean_matched_kmer_freq(&self, intersection: &HashSet<u64>) -> f64 {
        if intersection.is_empty() {
            return 0.0;
        }

        // Calculate average frequency of intersecting k-mers in the database
        // Sequential iter: intersections are typically small (< 200 elements); rayon
        // thread-pool overhead dominates for small collections.
        let sum_freq: f64 = intersection
            .iter()
            .map(|&hashval| {
                // Get frequency of this k-mer in the database (how many signatures contain it)
                let db_frequency =
                    self.stats.kmer_frequencies.get(&hashval).copied().unwrap_or(1) as f64;
                let total_signatures = self.stats.total_signatures as f64;

                // Normalize frequency to [0,1] range
                db_frequency / total_signatures
            })
            .sum();

        // Return average instead of sum
        sum_freq / intersection.len() as f64
    }

    /// Sum of target-DB frequencies for matched k-mers: Σ freq_target[h]/N over intersection.
    /// Takes the pre-computed intersection set directly.
    /// Higher = matched k-mers are collectively more common in the target DB.
    fn calculate_sum_matched_kmer_freq(&self, intersection: &HashSet<u64>) -> f64 {
        let total_signatures = self.stats.total_signatures as f64;
        intersection
            .iter()
            .map(|&hashval| {
                let db_frequency =
                    self.stats.kmer_frequencies.get(&hashval).copied().unwrap_or(1) as f64;
                db_frequency / total_signatures
            })
            .sum()
    }

    /// Sum of target-DB frequencies for matched k-mers, computed from two sketches directly.
    /// Equivalent to `calculate_sum_matched_kmer_freq` but takes sketches instead of a pre-computed
    /// intersection (useful when no intersection HashSet is available).
    pub fn calculate_sum_matched_kmer_freq_from_sketches(
        &self,
        query_sketch: &ProteinSketch,
        target_sketch: &ProteinSketch,
    ) -> f64 {
        let query_mins = query_sketch.signature().minhash.mins();
        let target_mins = target_sketch.mins_as_set();
        let total_signatures = self.stats.total_signatures as f64;

        // Sum database frequencies only for k-mers that actually match (sequential iter)
        query_mins
            .iter()
            .filter(|&&hashval| target_mins.contains(&hashval))
            .map(|&hashval| {
                let db_frequency =
                    self.stats.kmer_frequencies.get(&hashval).copied().unwrap_or(1) as f64;
                db_frequency / total_signatures
            })
            .sum()
    }

    /// Expected number of shared k-mers by chance: Σ freq_target[h]/N over ALL query k-mers.
    ///
    /// For each k-mer in this specific query, sums its probability of appearing in a random
    /// target sequence. Query-dependent (changes with different query sequences).
    /// Compare to n_intersecting_hashes: if observed >> expected, match is significant.
    pub fn calculate_expected_shared_kmers(
        &self,
        query_sketch: &ProteinSketch,
        _target_sketch: &ProteinSketch,
    ) -> f64 {
        let query_mins = query_sketch.signature().minhash.mins();
        let total_signatures = self.stats.total_signatures as f64;

        // For each k-mer in the query (regardless of whether it's in target),
        // calculate the expected probability it would appear in the target.
        // Sequential iter: rayon overhead is unjustified for per-query-sized collections.
        query_mins
            .iter()
            .map(|&hashval| {
                let db_frequency =
                    self.stats.kmer_frequencies.get(&hashval).copied().unwrap_or(1) as f64;
                db_frequency / total_signatures
            })
            .sum()
    }

    /// Get the underlying index
    pub fn index(&self) -> &ProteomeIndex {
        &self.index
    }

    /// Get search statistics
    pub fn stats(&self) -> &SearchStats {
        &self.stats
    }
}

/// Inner implementation: calculate similarity given pre-computed mins sets and intersection.
///
/// WHY: Called from both `calculate_similarity` (which computes the HashSets itself) and
/// `ProteinSearcher::compare` (which has already computed them for candidate pre-filtering).
/// Sharing one implementation prevents the intersection and size calculations from being
/// redundantly repeated across the call stack.
fn calculate_similarity_from_precomputed(
    query: &ProteinSketch,
    query_mins: &HashSet<u64>,
    target: &ProteinSketch,
    target_mins: &HashSet<u64>,
    intersection: &HashSet<u64>,
) -> Option<SearchResult> {
    let n_intersecting_hashes = intersection.len();
    if n_intersecting_hashes == 0 {
        return None;
    }

    let query_size = query_mins.len();
    let target_size = target_mins.len();
    let union_size = query_size + target_size - n_intersecting_hashes;

    let containment = n_intersecting_hashes as f64 / query_size as f64;
    let jaccard = n_intersecting_hashes as f64 / union_size as f64;
    let containment_target_in_query = n_intersecting_hashes as f64 / target_size as f64;
    let max_containment = containment.max(containment_target_in_query);

    let query_abunds = query.signature().minhash.abunds();
    let target_abunds = target.signature().minhash.abunds();

    let (average_abund, median_abund, std_abund) = if let (Some(qa), Some(ta)) =
        (query_abunds.as_ref(), target_abunds.as_ref())
    {
        let query_mins_sorted = query.signature().minhash.mins();
        let target_mins_sorted = target.signature().minhash.mins();
        significance::abundance_stats(intersection, &query_mins_sorted, qa, &target_mins_sorted, ta)
    } else {
        (1.0, 1.0, 0.0)
    };

    let f_weighted_target_in_query = significance::weighted_fraction_target_in_query(
        query_abunds.as_deref(),
        target_abunds.as_deref(),
    );

    let matched_regions = find_matched_regions(query, target, intersection);

    // Number of distinct positions a region could have started at in this query. Regions are
    // chosen after the fact (the best-looking gapless run is kept), so a per-region score
    // needs to be corrected by how many candidate positions it was chosen from. This is that
    // count.
    //
    // saturating_sub guards seq.len() - ksize from underflowing (both are usize; a plain `-`
    // would panic in debug and wrap to a huge number in release) if seq is ever shorter than
    // ksize. In practice that cannot happen here: this function already returned above when
    // the intersection is empty, and a non-empty intersection means the query produced at
    // least one k-mer, which means seq.len() >= ksize. The saturating form is a defensive
    // floor for that invariant, not a case this code path exercises.
    let region_search_space = query
        .get_raw_sequence()
        .map(|seq| seq.len().saturating_sub(query.protein_ksize() as usize) + 1)
        .unwrap_or(0);

    Some(SearchResult {
        query_name: query.signature().name.clone(),
        query_md5: query.signature().md5sum.clone(),
        target_name: target.signature().name.clone(),
        target_md5: target.signature().md5sum.clone(),
        containment,
        n_intersecting_hashes,
        ksize: query.protein_ksize(),
        scaled: query.signature().minhash.scaled(),
        moltype: query.moltype().to_string(),
        jaccard,
        max_containment,
        average_abund,
        median_abund,
        std_abund,
        containment_target_in_query,
        f_weighted_target_in_query,
        query_tfidf: 0.0,                 // requires database context
        mean_matched_kmer_freq: 0.0,      // requires database context
        sum_matched_kmer_freq: 0.0,       // requires database context
        query_expected_shared_kmers: 0.0, // requires database context
        query_enrichment: 0.0,            // requires database context
        joint_kmer_freq: 0.0,             // requires two-pass query frequencies
        query_poisson_pvalue: 1.0,        // requires database context
        region_search_space,
        db_n_targets: 0,  // requires database context
        db_n_kmers: 0,    // requires database context
        run_n_queries: 0, // requires a search run
        matched_regions,
    })
}

/// Calculate similarity between two protein sketches without requiring database context.
///
/// This is a standalone function for 1v1 comparisons that calculates all basic similarity
/// metrics (containment, jaccard, abundance statistics, matched regions) without needing
/// a `ProteinSearcher` instance. For database-specific metrics (TF-IDF, overlap probability),
/// use `ProteinSearcher::compare` instead.
///
/// # Why Standalone?
///
/// This function doesn't require any state from `ProteinSearcher`. It only operates on the
/// sketches provided, making it easier to test and more reusable. This follows idiomatic Rust:
/// functions that don't need state should be standalone. The pattern matches `find_matched_regions`.
///
/// # Arguments
/// * `query` - Query protein sketch
/// * `target` - Target protein sketch to compare against
///
/// # Returns
/// `Some(SearchResult)` if there's any intersection between the sketches, `None` otherwise.
/// TF-IDF is set to 0.0 and overlap probability is set to 1.0, as these metrics require
/// database context and are meaningless for 1v1 comparisons.
///
/// # Example
/// ```
/// use kmerseek::sketch::ProteinSketch;
/// use kmerseek::search::calculate_similarity;
///
/// let query = ProteinSketch::from_protein_sequence("query", "ATCGATCG", 10, 1, "hp_lehninger2").unwrap();
/// let target = ProteinSketch::from_protein_sequence("target", "ATCGATCG", 10, 1, "hp_lehninger2").unwrap();
/// let result = calculate_similarity(&query, &target);
/// ```
#[must_use]
pub fn calculate_similarity(query: &ProteinSketch, target: &ProteinSketch) -> Option<SearchResult> {
    let query_mins = query.mins_as_set();
    let target_mins = target.mins_as_set();
    let intersection: HashSet<u64> = query_mins.intersection(&target_mins).cloned().collect();
    calculate_similarity_from_precomputed(query, &query_mins, target, &target_mins, &intersection)
}

/// Every (query position, target position, hash) triple where a shared k-mer starts, one per
/// pairing of a hash's query starts with its target starts. Unsorted.
fn shared_position_pairs(
    query_sketch: &ProteinSketch,
    target_sketch: &ProteinSketch,
    intersection: &HashSet<u64>,
) -> Vec<(usize, usize, u64)> {
    // Build mapping from hashval to positions for both query and target
    // WHY: We need to maintain correspondence between query and target positions for each
    // k-mer hash. This allows us to find the correct target region for each query region.
    let mut hashval_to_query_positions: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut hashval_to_target_positions: HashMap<u64, Vec<usize>> = HashMap::new();

    for &hashval in intersection {
        if let (Some(query_poss_raw), Some(target_poss_raw)) = (
            query_sketch.kmer_positions().get(&hashval),
            target_sketch.kmer_positions().get(&hashval),
        ) {
            let mut query_poss = query_poss_raw.clone();
            query_poss.sort();
            query_poss.dedup();
            hashval_to_query_positions.insert(hashval, query_poss);

            let mut target_poss = target_poss_raw.clone();
            target_poss.sort();
            target_poss.dedup();
            hashval_to_target_positions.insert(hashval, target_poss);
        }
    }

    // Build (query_pos, target_pos) pairs for each hashval using the maps built above.
    let mut query_target_pairs: Vec<(usize, usize, u64)> = Vec::new();
    for (&hashval, query_poss) in &hashval_to_query_positions {
        if let Some(target_poss) = hashval_to_target_positions.get(&hashval) {
            for &qpos in query_poss {
                for &tpos in target_poss {
                    query_target_pairs.push((qpos, tpos, hashval));
                }
            }
        }
    }
    query_target_pairs
}

/// Whether two stored encoded residues agree under `moltype`: the same class, or an
/// ambiguous letter on either side that can encode to the other side's class. A protein20
/// sketch stores no encoded sequence, so its raw residues are compared, and there an
/// ambiguous letter agrees with either residue it stands for.
fn residues_agree(moltype: &str) -> impl Fn(u8, u8) -> bool {
    let encode = residue_encoder(moltype);
    move |query, target| encoded_residues_agree(query, target, &encode)
}

/// Whether two stored encoded regions of equal length agree residue by residue.
fn encoded_regions_agree(query: &str, target: &str, agree: &impl Fn(u8, u8) -> bool) -> bool {
    query.len() == target.len() && query.bytes().zip(target.bytes()).all(|(q, t)| agree(q, t))
}

/// The maximal stretch of agreeing encoded residues on one diagonal through the seed k-mer
/// at (`qpos`, `tpos`), as `[start, end)` in query coordinates. `None` when the seed window
/// itself disagrees, which can only be a hash collision now that an ambiguous residue
/// agrees with either class it stands for.
fn exact_run_around(
    query: &[u8],
    target: &[u8],
    qpos: usize,
    tpos: usize,
    ksize: usize,
    residues_agree: &impl Fn(u8, u8) -> bool,
) -> Option<(usize, usize)> {
    let offset = tpos as isize - qpos as isize;
    let agree = |i: usize| {
        let j = i as isize + offset;
        i < query.len()
            && j >= 0
            && (j as usize) < target.len()
            && residues_agree(query[i], target[j as usize])
    };
    if !(qpos..qpos + ksize).all(agree) {
        return None;
    }
    let mut start = qpos;
    while start > 0 && agree(start - 1) {
        start -= 1;
    }
    let mut end = qpos + ksize;
    while agree(end) {
        end += 1;
    }
    Some((start, end))
}

/// Region detection when the sketches keep only a `1/scaled` sample of their k-mers.
///
/// With every k-mer present, `find_matched_regions` reads a region straight off the shared
/// k-mers: consecutive shared starts on one diagonal are an exact match, and its length is
/// the count plus `ksize - 1`. A sampled sketch keeps k-mers by hash value, so the shared
/// starts on a diagonal sit about `scaled` apart and are rarely adjacent. That rule would
/// turn every match into a scatter of single k-mers.
///
/// Here each shared k-mer is a seed instead. The sequences are stored, so the seed's diagonal
/// is walked outward while the residues agree, which recovers the full exact match the seed
/// sits in, including the k-mers the sample dropped. A seed inside a run already emitted is
/// skipped. Every region returned is therefore a maximal exact match of at least `ksize`
/// residues holding at least one sampled shared k-mer.
///
/// That is the dense path's set minus any match the sample missed entirely, with one
/// difference: the walk compares residues, while the dense path chains sketched k-mers. A
/// window the sketch never held (one rejected by `disambiguate_kmer`, or dropped as low
/// complexity) breaks a dense run but not a walk, so a run through such a window is reported
/// here as one region where the dense path reports two.
///
/// `n_shared` counts only the sampled k-mers inside the run, since those are the observations
/// the Poisson test's expectation is summed over.
///
/// The walk uses the encoded sequence, falling back to the raw one when no encoded copy is
/// stored (the full alphabet encodes to itself, so `protein20` sketches keep only the raw
/// sequence). `moltype_seq` is empty in that case, as on the dense path. Returns nothing
/// when either sketch stores no sequence at all, since the walk needs one.
fn find_sampled_regions(
    query_sketch: &ProteinSketch,
    target_sketch: &ProteinSketch,
    intersection: &HashSet<u64>,
) -> Vec<MatchedRegion> {
    let ksize = query_sketch.protein_ksize() as usize;
    let (Some(query_raw), Some(target_raw)) =
        (query_sketch.get_raw_sequence(), target_sketch.get_raw_sequence())
    else {
        return Vec::new();
    };
    let has_encoded = query_sketch.get_moltype_sequence().is_some();
    let query_encoded = query_sketch.get_moltype_sequence().unwrap_or(query_raw);
    let target_encoded = target_sketch.get_moltype_sequence().unwrap_or(target_raw);
    let query_name = query_sketch.signature().name.clone();
    let target_name = target_sketch.signature().name.clone();
    let moltype = query_sketch.moltype().clone();
    let agree = residues_agree(&moltype.to_string());

    // Seeds grouped by diagonal (target start minus query start), sorted within each.
    let mut seeds_by_diagonal: BTreeMap<isize, Vec<usize>> = BTreeMap::new();
    for (qpos, tpos, _) in shared_position_pairs(query_sketch, target_sketch, intersection) {
        seeds_by_diagonal.entry(tpos as isize - qpos as isize).or_default().push(qpos);
    }

    let mut regions = Vec::new();
    for (diagonal, mut seeds) in seeds_by_diagonal {
        seeds.sort_unstable();
        seeds.dedup();
        // Query positions below this start inside a run already emitted on this diagonal.
        let mut covered_until = 0;
        for (i, &qpos) in seeds.iter().enumerate() {
            if qpos < covered_until {
                continue;
            }
            let tpos = (qpos as isize + diagonal) as usize;
            let Some((start, end)) = exact_run_around(
                query_encoded.as_bytes(),
                target_encoded.as_bytes(),
                qpos,
                tpos,
                ksize,
                &agree,
            ) else {
                continue;
            };
            let last_start = end - ksize;
            let n_shared = seeds[i..].iter().take_while(|&&p| p <= last_start).count();
            covered_until = last_start + 1;

            let target_start = (start as isize + diagonal) as usize;
            let target_end = target_start + (end - start);
            regions.push(MatchedRegion {
                query_name: query_name.clone(),
                start: start as u32,
                end: end as u32,
                subseq: query_raw[start..end].to_string(),
                target_name: target_name.clone(),
                target_start: target_start as u32,
                target_end: target_end as u32,
                target_subseq: target_raw[target_start..target_end].to_string(),
                moltype: moltype.clone(),
                moltype_seq: if has_encoded {
                    target_encoded[target_start..target_end].to_string()
                } else {
                    String::new()
                },
                length: (end - start) as u32,
                n_shared: n_shared as u32,
                n_mismatches: 0,
                expected_shared_kmers: 0.0,
                poisson_score: 0.0,
                tail_probability: 1.0,
                enrichment: 0.0,
                tfidf: 0.0,
                mean_idf: 0.0,
                ka_u: 0.0,
                ka_lambda: 0.0,
                ka_bits: 0.0,
                n_chained: 1,
                evalue: f64::INFINITY,
            });
        }
    }
    regions.sort_by_key(|r| r.start);
    regions
}

/// Find all consecutive matched regions of k-mer overlap between a query and target sequences
///
/// WHY: This is a standalone function because it doesn't require any state from ProteinSearcher.
/// It only operates on the sketches and intersection provided. This makes it easier to test and
/// more reusable. This is idiomatic Rust - functions that don't need state should be standalone.
///
/// # Panics
///
/// If the two sketches differ in k-mer size, alphabet, or scaled factor. Sketches built by
/// `search` share all three with the index; `--query-is-index` checks them before searching.
#[must_use = "matched regions should be used to analyze query-target alignments"]
pub fn find_matched_regions(
    query_sketch: &ProteinSketch,
    target_sketch: &ProteinSketch,
    intersection: &HashSet<u64>,
) -> Vec<MatchedRegion> {
    // Ensure that query and target protein sketches are the same ksize
    assert_eq!(query_sketch.protein_ksize(), target_sketch.protein_ksize());
    let ksize = query_sketch.protein_ksize() as usize;
    let query_name = query_sketch.signature().name.clone();
    let target_name = target_sketch.signature().name.clone();

    // Ensure that both query and target have the same moltypes
    assert_eq!(query_sketch.moltype(), target_sketch.moltype());
    let moltype = query_sketch.moltype().clone();
    let agree = residues_agree(&moltype.to_string());

    // A sampled sketch has too few adjacent shared k-mers for the consecutive-position rule
    // below; see find_sampled_regions.
    assert_eq!(query_sketch.scaled(), target_sketch.scaled());
    if query_sketch.scaled() > 1 {
        return find_sampled_regions(query_sketch, target_sketch, intersection);
    }

    let mut query_target_pairs = shared_position_pairs(query_sketch, target_sketch, intersection);

    if query_target_pairs.is_empty() {
        return Vec::new();
    }

    // Sort by diagonal (target position minus query position), then by query position, so
    // the walk below only ever compares neighbours on one diagonal. Sorting by query
    // position alone interleaved diagonals whenever a k-mer occurs more than once in the
    // target: pairs (q, t1), (q, t2), (q+1, t1+1) put (q, t2) between the two that continue
    // a run, the consecutive check failed, and one exact match came out as several regions
    // (BCL-2 against its own sequence at hp k=12 reported 144 residues, not 239).
    let diagonal = |p: &(usize, usize, u64)| p.1 as isize - p.0 as isize;
    query_target_pairs.sort_by(|a, b| diagonal(a).cmp(&diagonal(b)).then_with(|| a.0.cmp(&b.0)));

    // Find all consecutive regions where both query and target positions are consecutive
    let mut consecutive_regions = Vec::new();
    let mut i: usize = 0;

    while i < query_target_pairs.len() {
        let (query_start_pos, target_start_pos, _) = query_target_pairs[i];
        let mut consecutive_count: usize = 1;
        let mut j: usize = i + 1;

        // Find consecutive pairs where both query and target positions increment by 1
        // WHY: We need both query and target to be consecutive to form a valid matched region.
        // This ensures we match the correct target region to each query region.
        while j < query_target_pairs.len() {
            let (prev_qpos, prev_tpos, _) = query_target_pairs[j - 1];
            let (curr_qpos, curr_tpos, _) = query_target_pairs[j];

            // Check if both query and target positions are consecutive
            if curr_qpos == prev_qpos + 1 && curr_tpos == prev_tpos + 1 {
                consecutive_count += 1;
                j += 1;
            } else {
                break;
            }
        }

        // Calculate end positions
        let query_end_pos = query_start_pos + consecutive_count + ksize - 1;
        let target_end_pos = target_start_pos + consecutive_count + ksize - 1;

        // Get sequences for extraction
        // WHY: We handle missing sequences gracefully instead of panicking. If sequences aren't
        // stored (because store_raw_sequences was false), we can't extract matched regions, so
        // we return an empty vector. This ensures consistency with the index configuration.
        let (query_raw_sequence, target_raw_sequence) =
            match (query_sketch.get_raw_sequence(), target_sketch.get_raw_sequence()) {
                (Some(q), Some(t)) => (q, t),
                _ => {
                    // Sequences not available - return empty regions
                    // WHY: This can happen when store_raw_sequences is false. Instead of panicking,
                    // we return empty regions, which is the correct behavior when sequences aren't stored.
                    return Vec::new();
                }
            };

        // Extract subsequences using correct positions
        // WHY: Query subsequence uses query positions, target subsequence uses target positions.
        // This is the fix for the bug where both were using query positions.
        // We add bounds checking to prevent panics from out-of-bounds slicing.
        if query_end_pos > query_raw_sequence.len() || target_end_pos > target_raw_sequence.len() {
            // Bounds check failed - skip this region
            i = j;
            continue;
        }
        let query_subseq = &query_raw_sequence[query_start_pos..query_end_pos];
        let target_subseq = &target_raw_sequence[target_start_pos..target_end_pos];

        // Get moltype sequences for validation
        // WHY: We handle missing encoded sequences gracefully. If they're not available,
        // we skip validation but still extract the regions. This ensures consistency with
        // the index configuration where sequences might not be stored.
        let (target_moltype_sequence, query_moltype_sequence) =
            match (target_sketch.get_moltype_sequence(), query_sketch.get_moltype_sequence()) {
                (Some(t), Some(q)) => (t, q),
                _ => {
                    // Encoded sequences not available - skip validation but still extract regions
                    // WHY: This can happen when store_raw_sequences is false. We can still extract
                    // regions from raw sequences, but we skip the moltype validation step.
                    // We reuse the already-extracted subsequences to avoid duplicate bounds checking.
                    // Note: query_subseq and target_subseq are already defined above, so we use them directly.

                    consecutive_regions.push(MatchedRegion {
                        query_name: query_name.clone(),
                        start: query_start_pos as u32,
                        end: query_end_pos as u32,
                        subseq: query_subseq.to_string(),
                        target_name: target_name.clone(),
                        target_start: target_start_pos as u32,
                        target_end: target_end_pos as u32,
                        target_subseq: target_subseq.to_string(),
                        moltype: moltype.clone(),
                        moltype_seq: String::new(), // Empty since we don't have encoded sequence
                        length: (query_end_pos - query_start_pos) as u32,
                        n_shared: consecutive_count as u32,
                        n_mismatches: 0,
                        expected_shared_kmers: 0.0,
                        poisson_score: 0.0,
                        tail_probability: 1.0,
                        enrichment: 0.0,
                        tfidf: 0.0,
                        mean_idf: 0.0,
                        ka_u: 0.0,
                        ka_lambda: 0.0,
                        ka_bits: 0.0,
                        n_chained: 1,
                        evalue: f64::INFINITY,
                    });

                    i = j;
                    continue;
                }
            };

        // Extract moltype subsequences using correct positions
        // WHY: We add bounds checking to prevent panics from out-of-bounds slicing.
        if query_end_pos > query_moltype_sequence.len()
            || target_end_pos > target_moltype_sequence.len()
        {
            // Bounds check failed - skip this region
            i = j;
            continue;
        }
        let query_moltype_seq = &query_moltype_sequence[query_start_pos..query_end_pos];
        let target_moltype_seq = &target_moltype_sequence[target_start_pos..target_end_pos];

        // The two encoded regions must agree residue by residue, with an ambiguous letter
        // agreeing with either class it stands for. They share the same k-mers, so a
        // disagreement can only be a hash collision.
        if !encoded_regions_agree(query_moltype_seq, target_moltype_seq, &agree) {
            i = j;
            continue;
        }

        consecutive_regions.push(MatchedRegion {
            query_name: query_name.clone(),
            start: query_start_pos as u32,
            end: query_end_pos as u32,
            subseq: query_subseq.to_string(),
            target_name: target_name.clone(),
            target_start: target_start_pos as u32,
            target_end: target_end_pos as u32,
            target_subseq: target_subseq.to_string(),
            moltype: moltype.clone(),
            moltype_seq: target_moltype_seq.to_string(),
            length: (query_end_pos - query_start_pos) as u32,
            n_shared: consecutive_count as u32,
            n_mismatches: 0,
            expected_shared_kmers: 0.0,
            poisson_score: 0.0,
            tail_probability: 1.0,
            enrichment: 0.0,
            tfidf: 0.0,
            mean_idf: 0.0,
            ka_u: 0.0,
            ka_lambda: 0.0,
            ka_bits: 0.0,
            n_chained: 1,
            evalue: f64::INFINITY,
        });

        i = j;
    }

    // Sort regions by query position (earliest first)
    // WHY: Regions are already in position order from the loop, but we sort explicitly to ensure
    // correctness and make the ordering clear. Sorting by position makes it easier to understand
    // the sequence of matches along the query sequence.
    //
    // NOTE: We do NOT filter out overlapping regions because query positions may overlap while
    // target positions differ. For example, a 12-mer match at query 170:182 and target 81:93 is
    // distinct from a 19-mer match at query 162:181 and target 138:157, even though the query
    // positions overlap. All matches are reported because they represent different alignments.
    consecutive_regions.sort_by_key(|a| a.start);

    consecutive_regions
}

/// Walk one direction from a seed edge along the encoded sequences with X-drop.
///
/// `positions` yields (query index, target index) pairs stepping away from the seed. Returns
/// how many positions the best-scoring extension covers.
fn xdrop_walk(
    q: &[u8],
    t: &[u8],
    positions: impl Iterator<Item = (usize, usize)>,
    params: ExtensionParams,
) -> usize {
    let mut score = 0.0;
    let mut best = 0.0;
    let mut best_len = 0;
    for (n, (qi, ti)) in positions.enumerate() {
        score += if q[qi] == t[ti] { 1.0 } else { -params.mismatch_penalty };
        if score > best {
            best = score;
            best_len = n + 1;
        } else if best - score > params.xdrop {
            break;
        }
    }
    best_len
}

/// Grow each exact-run region outward with mismatches allowed, then merge regions on one
/// diagonal whose extended spans touch. See `ExtensionParams` for the why.
///
/// `n_shared` is carried from the seeds (summed on merge) and never recomputed from the
/// extended length, so the Poisson test keeps counting shared k-mers rather than residues.
/// `n_mismatches` is recounted over the final span. Regions come back sorted by query start.
///
/// Both sketches must carry encoded and raw sequences; without them the regions are returned
/// unchanged, which is also what `find_matched_regions` does in that case.
pub fn extend_regions(
    regions: Vec<MatchedRegion>,
    query_sketch: &ProteinSketch,
    target_sketch: &ProteinSketch,
    params: ExtensionParams,
) -> Vec<MatchedRegion> {
    if regions.is_empty() || params.mismatch_penalty <= 0.0 {
        return regions;
    }
    let (Some(q_enc), Some(t_enc), Some(q_raw), Some(t_raw)) = (
        query_sketch.get_class_sequence(),
        target_sketch.get_class_sequence(),
        query_sketch.get_raw_sequence(),
        target_sketch.get_raw_sequence(),
    ) else {
        return regions;
    };
    let (q, t) = (q_enc.as_bytes(), t_enc.as_bytes());

    // Extend every seed on its own diagonal.
    let mut extended: Vec<MatchedRegion> = regions
        .into_iter()
        .map(|mut r| {
            let (qs, qe) = (r.start as usize, r.end as usize);
            let (ts, te) = (r.target_start as usize, r.target_end as usize);
            let right = xdrop_walk(q, t, (qe..q.len()).zip(te..t.len()), params);
            let left = xdrop_walk(q, t, (0..qs).rev().zip((0..ts).rev()), params);
            r.start = (qs - left) as u32;
            r.end = (qe + right) as u32;
            r.target_start = (ts - left) as u32;
            r.target_end = (te + right) as u32;
            r
        })
        .collect();

    // Merge on (diagonal, start): two seeds whose extensions meet are one match.
    let diagonal = |r: &MatchedRegion| r.target_start as i64 - r.start as i64;
    extended.sort_by_key(|r| (diagonal(r), r.start));
    let mut merged: Vec<MatchedRegion> = Vec::with_capacity(extended.len());
    for r in extended {
        if let Some(last) = merged.last_mut() {
            if diagonal(last) == diagonal(&r) && r.start <= last.end {
                if r.end > last.end {
                    last.end = r.end;
                    last.target_end = r.target_end;
                }
                last.n_shared += r.n_shared;
                continue;
            }
        }
        merged.push(r);
    }

    // Rebuild the derived fields over the final spans.
    for r in merged.iter_mut() {
        let (qs, qe) = (r.start as usize, r.end as usize);
        let (ts, te) = (r.target_start as usize, r.target_end as usize);
        r.length = (qe - qs) as u32;
        r.n_mismatches = q[qs..qe].iter().zip(&t[ts..te]).filter(|(a, b)| a != b).count() as u32;
        r.subseq = q_raw[qs..qe].to_string();
        r.target_subseq = t_raw[ts..te].to_string();
        r.moltype_seq = t_enc[ts..te].to_string();
    }
    merged.sort_by_key(|r| r.start);
    merged
}

/// Chain extended regions that sit on one diagonal within `chain_max_gap` residues of each
/// other into one region, scored with Karlin-Altschul sum statistics.
///
/// Why. A region is one gapless run and the benchmark's transfer rule labels it only when
/// it covers half the target domain, so a 200-residue domain needs a single clean
/// 100-residue run: the exact-match ceiling in a different coat. Two runs on one diagonal
/// separated by a stretch the X-drop would not cross are one alignment with a bad patch in
/// it, and Karlin & Altschul 1993 give the statistic for the sum of their scores.
///
/// The chain spans from the first region's start to the last region's end on both
/// sequences; `n_shared` is summed, `n_mismatches` is recounted over the whole span (the
/// gap's disagreements included, so the span is honest about what it contains),
/// `n_chained` is the number of members. Its E-value is the sum P-value times the number
/// of targets: each member's normalised score is lambda S_i - ln(K m n_t) with n_t the
/// target length, so for a single region this reduces to the per-region E within the
/// approximation n = N n_t. Members must be colinear (each starts after the previous one
/// ends, on both sequences) and within `chain_max_shift` diagonals of each other, so a
/// chain tolerates a net indel up to that size between members. With a shift the query
/// and target spans differ in length; `length` is the query span, `n_mismatches` is the
/// members' sum (the gaps are not scored either way) and the sum statistic is what
/// carries the evidence.
#[allow(clippy::too_many_arguments)]
pub fn chain_regions(
    regions: Vec<MatchedRegion>,
    q: &[u8],
    t: &[u8],
    q_raw: Option<&str>,
    t_raw: Option<&str>,
    params: ExtensionParams,
    m: f64,
    n_t: f64,
    n_targets: f64,
) -> Vec<MatchedRegion> {
    if regions.len() < 2 {
        return regions;
    }
    let (Some(q_raw), Some(t_raw)) = (q_raw, t_raw) else {
        return regions;
    };
    let diagonal = |r: &MatchedRegion| r.target_start as i64 - r.start as i64;
    let raw_score = |r: &MatchedRegion| {
        (r.length - r.n_mismatches) as f64 - params.mismatch_penalty * r.n_mismatches as f64
    };
    let mut sorted = regions;
    sorted.sort_by_key(|r| (r.start, r.target_start));

    let mut out: Vec<MatchedRegion> = Vec::with_capacity(sorted.len());
    let mut chain: Vec<MatchedRegion> = Vec::new();
    let ln_kmn = (params.ka_k * m * n_t).ln();
    let flush = |chain: &mut Vec<MatchedRegion>, out: &mut Vec<MatchedRegion>| {
        if chain.is_empty() {
            return;
        }
        if chain.len() == 1 {
            out.push(chain.pop().unwrap());
            return;
        }
        let r = chain.len() as u32;
        let first = &chain[0];
        let last = &chain[chain.len() - 1];
        let (qs, qe) = (first.start as usize, last.end as usize);
        let (ts, te) = (first.target_start as usize, last.target_end as usize);
        // One lambda for the chain, from the composition of the span it covers, so the sum
        // statistic adds up members measured on one scale. Taken over the chained span
        // rather than per member because that span is what the row reports, gaps included.
        let ka_u = region_match_probability(&q[qs..qe], &t[ts..te]);
        let ka_lambda =
            karlin_altschul_lambda(ka_u, params.mismatch_penalty) * params.ka_lambda_scale;
        let t_sum: f64 = chain.iter().map(|x| ka_lambda * raw_score(x) - ln_kmn).sum();
        let p = karlin_altschul_sum_p(t_sum, r);
        let mut merged = first.clone();
        merged.end = qe as u32;
        merged.target_end = te as u32;
        merged.length = (qe - qs) as u32;
        merged.n_shared = chain.iter().map(|x| x.n_shared).sum();
        merged.n_mismatches = if qe - qs == te - ts {
            q[qs..qe].iter().zip(&t[ts..te]).filter(|(a, b)| a != b).count() as u32
        } else {
            chain.iter().map(|x| x.n_mismatches).sum()
        };
        merged.n_chained = r;
        merged.subseq = q_raw[qs..qe].to_string();
        merged.target_subseq = t_raw[ts..te].to_string();
        merged.moltype_seq = String::from_utf8_lossy(&t[ts..te]).into_owned();
        merged.ka_u = ka_u;
        merged.ka_lambda = ka_lambda;
        if ka_lambda > 0.0 {
            merged.evalue = p * n_targets;
            // Bits on the same scale as a single region: the sum statistic's -ln P, in bits.
            merged.ka_bits =
                if p > 0.0 { (-p.ln() / std::f64::consts::LN_2).max(0.0) } else { f64::INFINITY };
        } else {
            // The chained span's own composition is past the boundary: not assessable, the
            // same verdict a single region in that position gets.
            merged.evalue = f64::INFINITY;
            merged.ka_bits = 0.0;
        }
        // The chain's Poisson fields describe the first member only and would mislead; the
        // expectation is re-summed over the span by the caller if it needs it. Keep the
        // strongest member's tail so the pair-level filter sees the evidence it saw before.
        merged.poisson_score = chain.iter().map(|x| x.poisson_score).fold(0.0, f64::max);
        merged.tail_probability = chain.iter().map(|x| x.tail_probability).fold(1.0, f64::min);
        merged.expected_shared_kmers = chain.iter().map(|x| x.expected_shared_kmers).sum();
        merged.enrichment = fold_enrichment(merged.n_shared, merged.expected_shared_kmers);
        out.push(merged);
        chain.clear();
    };
    // Greedy colinear chaining in query order: a region joins the open chain when it
    // starts after the chain's last member on both sequences, within the gap cap on the
    // query, and within the diagonal band. Anything else closes the chain. Greedy, not
    // optimal: two interleaved chains on far-apart diagonals would be split at each
    // alternation, which is the conservative outcome.
    for r in sorted {
        if let Some(last) = chain.last() {
            let colinear = r.start >= last.end && r.target_start >= last.target_end;
            let gap_ok = colinear && r.start - last.end <= params.chain_max_gap;
            let shift_ok =
                (diagonal(last) - diagonal(&r)).unsigned_abs() <= params.chain_max_shift as u64;
            if !(gap_ok && shift_ok) {
                flush(&mut chain, &mut out);
            }
        }
        chain.push(r);
    }
    flush(&mut chain, &mut out);
    out.sort_by_key(|r| r.start);
    out
}

impl ProteinSearcher {
    // Note: find_signature_by_name and get_stored_encoded_sequence were removed as they were unused.
    // If needed in the future, they can be re-added.
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::ProteinSketch;
    use crate::tests::test_fixtures::{TEST_BLC2_FASTA, TEST_CED9_FASTA, TEST_FASTA_GZ};
    use approx::assert_relative_eq;
    use needletail::parse_fastx_file;
    use rstest::{fixture, rstest};
    use std::path::Path;
    use tempfile::TempDir;

    #[test]
    fn test_search_result_csv_from_result_and_region() {
        use crate::types::MolType;

        let region = MatchedRegion {
            query_name: "q".to_string(),
            start: 3,
            end: 9,
            subseq: "QSUBSEQ".to_string(),
            target_name: "t".to_string(),
            target_start: 11,
            target_end: 17,
            target_subseq: "TSUBSEQ".to_string(),
            moltype_seq: "hphph".to_string(),
            moltype: MolType::new("hp_lehninger2").unwrap(),
            length: 6,
            n_shared: 2,
            n_mismatches: 0,
            expected_shared_kmers: 2.0,
            poisson_score: 0.05,
            tail_probability: 0.89,
            enrichment: 1.5,
            tfidf: 7.0,
            mean_idf: 3.5,
            ka_u: 0.0,
            ka_lambda: 0.0,
            ka_bits: 0.0,
            n_chained: 1,
            evalue: f64::INFINITY,
        };

        let result = SearchResult {
            query_name: "query1".to_string(),
            query_md5: "qmd5".to_string(),
            target_name: "target1".to_string(),
            target_md5: "tmd5".to_string(),
            containment: 0.5,
            n_intersecting_hashes: 7,
            ksize: 5,
            scaled: 1,
            moltype: "hp_lehninger2".to_string(),
            jaccard: 0.25,
            max_containment: 0.6,
            average_abund: 2.0,
            median_abund: 2.0,
            std_abund: 0.5,
            containment_target_in_query: 0.4,
            f_weighted_target_in_query: 0.3,
            query_tfidf: 1.5,
            mean_matched_kmer_freq: 0.1,
            sum_matched_kmer_freq: 0.7,
            query_expected_shared_kmers: 3.0,
            query_enrichment: 2.33,
            joint_kmer_freq: 0.05,
            query_poisson_pvalue: 0.01,
            region_search_space: 300,
            db_n_targets: 25,
            db_n_kmers: 7629,
            run_n_queries: 4,
            matched_regions: vec![],
        };

        let row = SearchResultCsv::from_result_and_region(&result, &region, false);

        // Fields carried from the SearchResult.
        assert_eq!(row.query_name, "query1");
        assert_eq!(row.target_md5, "tmd5");
        assert_eq!(row.n_intersecting_hashes, 7);
        assert_eq!(row.ksize, 5);
        assert_eq!(row.containment, 0.5);
        assert_eq!(row.query_poisson_pvalue, 0.01);
        // Fields carried from the MatchedRegion.
        assert_eq!(row.region_start, 3);
        assert_eq!(row.region_end, 9);
        assert_eq!(row.region_subseq, "QSUBSEQ");
        assert_eq!(row.target_start, 11);
        assert_eq!(row.target_subseq, "TSUBSEQ");
        assert_eq!(row.moltype_seq, "hphph");
        assert_eq!(row.region_length, 6);
        // Region-scoped stat columns are copied straight from the region.
        assert_eq!(row.region_n_shared_kmers, 2);
        assert_eq!(row.region_expected_shared_kmers, 2.0);
        assert_eq!(row.region_poisson_score, 0.05);
        assert_eq!(row.region_tail_probability, 0.89);
        assert_eq!(row.region_enrichment, 1.5);
        assert_eq!(row.region_tfidf, 7.0);
        assert_eq!(row.region_mean_idf, 3.5);
        // region_search_space, db_n_targets, db_n_kmers, and run_n_queries travel as
        // separate columns, never folded into a p-value.
        assert_eq!(row.region_search_space, 300);
        assert_eq!(row.db_n_targets, 25);
        assert_eq!(row.db_n_kmers, 7629);
        assert_eq!(row.run_n_queries, 4);
    }

    #[allow(dead_code)] // Test data structure - fields may be used for comparison
    struct ExpectedSimilarity {
        ksize: usize,
        n_intersecting_hashes: usize,
        containment: f64,
        jaccard: f64,
        max_containment: f64,
        containment_target_in_query: f64,
        average_abund: f64,
        median_abund: f64,
        std_abund: f64,
        matched_regions_count: usize,
    }

    /// Expected matched region structure for testing
    struct ExpectedMatchedRegion {
        subseq: &'static str,
        moltype_seq: &'static str,
        target_subseq: &'static str,
        start: u32,
        end: u32,
        target_start: u32,
        target_end: u32,
    }

    /// Expected matched regions structure for testing
    #[allow(dead_code)] // Test data structure - fields may be used for comparison
    struct ExpectedMatchedRegions {
        ksize: usize,
        total_regions: usize,
        // Just check first, last, and any specific "landmark" regions
        first: ExpectedMatchedRegion,
        last: ExpectedMatchedRegion,
        // Optional: a specific region to find by subseq
        landmark: Option<ExpectedMatchedRegion>,
    }

    const MATCHED_REGIONS_K12: ExpectedMatchedRegions = ExpectedMatchedRegions {
        ksize: 12,
        total_regions: 13,
        first: ExpectedMatchedRegion {
            subseq: "FTHRIRQNGMEW",
            moltype_seq: "hppphppphhph",
            target_subseq: "FSRRYRRDFAEM",
            start: 87,
            end: 99,
            target_start: 103,
            target_end: 115,
        },
        last: ExpectedMatchedRegion {
            subseq: "GVVVCGRMMFSLK",
            moltype_seq: "hhhhphphhhphp",
            target_subseq: "LYGPSMRPLFDFS",
            start: 267,
            end: 280,
            target_start: 200,
            target_end: 213,
        },
        landmark: Some(ExpectedMatchedRegion {
            subseq: "QCPMSYGRLIGLISFGGFV",
            moltype_seq: "pphhphhphhhhhphhhhh",
            target_subseq: "RDGVNWGRIVAFFEFGGVM",
            start: 162,
            end: 181,
            target_start: 138,
            target_end: 157,
        }),
    };

    const MATCHED_REGIONS_K15: ExpectedMatchedRegions = ExpectedMatchedRegions {
        ksize: 15,
        total_regions: 1,
        first: ExpectedMatchedRegion {
            subseq: "QCPMSYGRLIGLISFGGFV",
            moltype_seq: "pphhphhphhhhhphhhhh",
            target_subseq: "RDGVNWGRIVAFFEFGGVM",
            start: 162,
            end: 181,
            target_start: 138,
            target_end: 157,
        },
        last: ExpectedMatchedRegion {
            subseq: "QCPMSYGRLIGLISFGGFV",
            moltype_seq: "pphhphhphhhhhphhhhh",
            target_subseq: "RDGVNWGRIVAFFEFGGVM",
            start: 162,
            end: 181,
            target_start: 138,
            target_end: 157,
        },
        landmark: None,
    };

    #[fixture]
    fn temp_dir() -> TempDir {
        TempDir::new().unwrap()
    }

    #[fixture]
    fn ced9_record() -> (String, String) {
        read_first_fasta_record(TEST_CED9_FASTA).unwrap()
    }

    #[fixture]
    fn bcl2_record() -> (String, String) {
        read_first_fasta_record(TEST_BLC2_FASTA).unwrap()
    }

    #[fixture]
    fn ced9_sketch_k12() -> ProteinSketch {
        // WHY: This fixture reads the FASTA file directly. While it could depend on ced9_record(),
        // rstest's #[case] attributes don't support fixtures with dependencies. However, rstest
        // still caches this fixture, so if multiple test cases use it, the FASTA is only read once.
        // This gives us the caching benefit while maintaining compatibility with #[case] attributes.
        let (name, seq) = read_first_fasta_record(TEST_CED9_FASTA).unwrap();
        ProteinSketch::from_protein_sequence(&name, &seq, 12, 1, "hp_lehninger2").unwrap()
    }

    #[fixture]
    fn ced9_sketch_k15() -> ProteinSketch {
        let (name, seq) = read_first_fasta_record(TEST_CED9_FASTA).unwrap();
        ProteinSketch::from_protein_sequence(&name, &seq, 15, 1, "hp_lehninger2").unwrap()
    }

    #[fixture]
    fn bcl2_sketch_k12() -> ProteinSketch {
        let (name, seq) = read_first_fasta_record(TEST_BLC2_FASTA).unwrap();
        ProteinSketch::from_protein_sequence(&name, &seq, 12, 1, "hp_lehninger2").unwrap()
    }

    #[fixture]
    fn bcl2_sketch_k15() -> ProteinSketch {
        let (name, seq) = read_first_fasta_record(TEST_BLC2_FASTA).unwrap();
        ProteinSketch::from_protein_sequence(&name, &seq, 15, 1, "hp_lehninger2").unwrap()
    }

    /// A sequence against itself is one exact match over its whole length, whatever k-mers
    /// repeat inside it. BCL-2's Ala/Pro/Gly loop repeats several hp 12-mers, which used to
    /// split this into a 45-residue and a 144-residue region.
    #[rstest]
    fn self_hit_is_one_full_length_region(bcl2_sketch_k12: ProteinSketch) {
        let intersection = bcl2_sketch_k12.intersect(&bcl2_sketch_k12);
        let regions = find_matched_regions(&bcl2_sketch_k12, &bcl2_sketch_k12, &intersection);
        let full: Vec<&MatchedRegion> = regions.iter().filter(|r| r.length == 239).collect();
        assert_eq!(
            full.len(),
            1,
            "regions: {:?}",
            regions.iter().map(|r| (r.start, r.end, r.target_start)).collect::<Vec<_>>()
        );
        assert_eq!(
            (full[0].start, full[0].end, full[0].target_start, full[0].target_end),
            (0, 239, 0, 239)
        );
        // The other 16 are repeats of the Ala/Pro/Gly loop (and the tail's copy of it), all
        // off the main diagonal.
        assert_eq!(regions.len(), 17);
        assert!(regions.iter().all(|r| r.length == 239 || r.start != r.target_start));
    }

    /// Two target positions for one query k-mer must not break the run the first one
    /// continues: (80, 253) and (81, 254) are consecutive on their diagonal even though
    /// (81, 170) sorts between them by query position.
    #[rstest]
    fn a_repeated_kmer_does_not_split_the_run_on_the_other_diagonal(
        bcl2_sketch_k12: ProteinSketch,
        ced9_sketch_k12: ProteinSketch,
    ) {
        let intersection = bcl2_sketch_k12.intersect(&ced9_sketch_k12);
        let regions = find_matched_regions(&bcl2_sketch_k12, &ced9_sketch_k12, &intersection);
        let run =
            regions.iter().find(|r| r.start == 80 && r.target_start == 253).expect("run at 80/253");
        assert_eq!((run.end, run.target_end, run.length), (93, 266, 13));
        assert!(regions.iter().any(|r| r.start == 81 && r.target_start == 170 && r.length == 12));
    }

    /// Two database entries with the same sequence are one sketch in the index; a query that
    /// hits it is reported under both names, and a query identical to it is not dropped as a
    /// self-match unless the search is all-vs-all.
    #[test]
    fn identical_entries_are_each_reported_and_the_query_itself_is_a_hit() -> Result<()> {
        let (bcl2_name, bcl2) = read_first_fasta_record(TEST_BLC2_FASTA)?;
        let (ced9_name, ced9) = read_first_fasta_record(TEST_CED9_FASTA)?;
        let temp_dir = TempDir::new()?;
        let fasta = temp_dir.path().join("db.fasta");
        std::fs::write(
            &fasta,
            format!(">{bcl2_name}\n{bcl2}\n>copy of BCL2\n{bcl2}\n>{ced9_name}\n{ced9}\n"),
        )?;
        let index = ProteomeIndex::new(temp_dir.path().join("db"), 12, 1, "hp_lehninger2", true)?;
        index.process_fasta(&fasta, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(index)?;

        let query =
            ProteinSketch::from_protein_sequence("bcl2 query", &bcl2, 12, 1, "hp_lehninger2")?;
        let mut names: Vec<String> = searcher
            .search_one(&query, &SearchFilters::default(), 1)
            .iter()
            .map(|r| r.target_name.clone())
            .collect();
        names.sort();
        assert_eq!(names, vec!["copy of BCL2".to_string(), bcl2_name.clone(), ced9_name.clone()]);
        let results = searcher.search_one(&query, &SearchFilters::default(), 1);
        let copy = results.iter().find(|r| r.target_name == "copy of BCL2").unwrap();
        let original = results.iter().find(|r| r.target_name == bcl2_name).unwrap();
        assert_eq!(copy.n_intersecting_hashes, original.n_intersecting_hashes);
        assert_eq!(copy.target_md5, original.target_md5);
        assert!(copy.matched_regions.iter().all(|m| m.target_name == "copy of BCL2"));

        // All-vs-all drops every query's hit on its own sketch, under either name.
        let all = searcher.search_all_vs_all(&SearchFilters::default())?;
        assert!(all.iter().all(|r| r.query_md5 != r.target_md5), "self-hits should be skipped");
        Ok(())
    }

    /// Read the first record from a FASTA file and return name and sequence.
    ///
    /// WHY: This helper function eliminates code duplication in tests. It provides a simple
    /// way to read FASTA records for testing purposes without the complexity of the full
    /// fasta module API.
    fn read_first_fasta_record<P: AsRef<Path>>(path: P) -> Result<(String, String)> {
        let mut reader = parse_fastx_file(path)
            .map_err(|e| anyhow::anyhow!("Failed to parse FASTA file: {}", e))?;

        let record = match reader.next() {
            Some(Ok(record)) => record,
            Some(Err(e)) => return Err(anyhow::anyhow!("Failed to read FASTA record: {}", e)),
            None => return Err(anyhow::anyhow!("No sequence found in FASTA file")),
        };

        let sequence = String::from_utf8(record.seq().to_vec())
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence: {}", e))?;
        let name = String::from_utf8(record.id().to_vec())
            .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in name: {}", e))?;

        Ok((name, sequence))
    }

    /// The first record whose id contains `id`, as (name, sequence).
    fn read_fasta_record<P: AsRef<Path>>(path: P, id: &str) -> Result<(String, String)> {
        let mut reader = parse_fastx_file(path)?;
        while let Some(record) = reader.next() {
            let record = record?;
            let name = String::from_utf8(record.id().to_vec())?;
            if name.contains(id) {
                return Ok((name, String::from_utf8(record.seq().to_vec())?));
            }
        }
        anyhow::bail!("no record with id containing {id}")
    }

    /// Test search functionality similar to the Python tests
    #[test]
    fn test_search_basic() -> Result<()> {
        // Create temporary directory for test data
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path();

        // Use BCL2 as query and CED9 as target - they should match via HP encoding
        // WHY: The compare() function skips self-matches by checking MD5 sums, so we
        // need to use different proteins. BCL2 and CED9 are known to match via HP encoding.
        let query_fasta = TEST_BLC2_FASTA;
        let target_fasta = TEST_CED9_FASTA;

        // Create target index (CED9)
        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(
            &target_index_path,
            15,              // ksize - k=15 is where BCL2/CED9 have good HP overlap
            1,               // scaled
            "hp_lehninger2", // moltype
            false,           // store_raw_sequences
        )?;

        target_index.process_fasta(target_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        // Create searcher
        let searcher = ProteinSearcher::new(target_index)?;

        // Create query index (BCL2)
        let query_index_path = temp_path.join("query_index");
        let query_index = ProteomeIndex::new(
            &query_index_path,
            15,              // ksize
            1,               // scaled
            "hp_lehninger2", // moltype
            false,           // store_raw_sequences
        )?;

        query_index.process_fasta(query_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        // Get query signatures
        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        assert!(!query_signatures.is_empty(), "Should have at least one query signature");

        // Perform search
        let results = searcher.search(&query_signatures, &SearchFilters::default())?;

        // Should find at least one match (BCL2 vs CED9 via HP encoding)
        assert!(!results.is_empty(), "Should find at least one match between BCL2 and CED9");

        // Check that the first result has reasonable values
        let first_result = &results[0];
        assert!(
            first_result.query_name.contains("BCL2_HUMAN")
                || first_result.query_name.contains("Q07817"),
            "Query should be BCL2, got: {}",
            first_result.query_name
        );
        assert!(
            first_result.target_name.contains("CED9_CAEEL")
                || first_result.target_name.contains("P41958"),
            "Target should be CED9, got: {}",
            first_result.target_name
        );
        // BCL2 and CED9 share HP k-mers, so containment should be > 0
        assert!(first_result.containment > 0.0, "Should have positive containment");
        assert!(first_result.jaccard > 0.0, "Should have positive jaccard");
        assert!(first_result.n_intersecting_hashes > 0, "Should have intersecting hashes");

        Ok(())
    }

    /// Each of `SearchFilters`'s rejection checks, exercised on its own, rejects an
    /// otherwise-real BCL2/CED9 match: threshold, min_shared_kmers, and (together, since
    /// either alone passing keeps the pair) the query p-value and region score. Every other
    /// field is left at its permissive default in each case, so the field under test is what
    /// causes the rejection, not some other, stricter field.
    #[test]
    fn test_search_filters_reject_candidates() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path();

        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 15, 1, "hp_lehninger2", false)?;
        target_index.process_fasta(
            TEST_CED9_FASTA,
            DEFAULT_PROGRESS_INTERVAL,
            DEFAULT_BATCH_SIZE,
        )?;
        let searcher = ProteinSearcher::new(target_index)?;

        let query_index_path = temp_path.join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, 15, 1, "hp_lehninger2", false)?;
        query_index.process_fasta(
            TEST_BLC2_FASTA,
            DEFAULT_PROGRESS_INTERVAL,
            DEFAULT_BATCH_SIZE,
        )?;
        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        // Sanity check: with permissive filters, this pair does match.
        let unfiltered = searcher.search(&query_signatures, &SearchFilters::default())?;
        assert!(!unfiltered.is_empty(), "BCL2 vs CED9 should match with no filtering");

        // containment is always <= 1.0, so this threshold rejects every candidate.
        let threshold_filtered = searcher.search(
            &query_signatures,
            &SearchFilters { threshold: 1.1, ..SearchFilters::default() },
        )?;
        assert!(threshold_filtered.is_empty(), "threshold: 1.1 should reject every candidate");

        // No real match shares more k-mers than usize::MAX.
        let min_shared_kmers_filtered = searcher.search(
            &query_signatures,
            &SearchFilters { min_shared_kmers: usize::MAX, ..SearchFilters::default() },
        )?;
        assert!(
            min_shared_kmers_filtered.is_empty(),
            "min_shared_kmers: usize::MAX should reject every candidate"
        );

        // max_query_pvalue: 0.0 can never be cleared (query_pvalue < 0.0 is impossible), and
        // min_region_score: INFINITY can never be cleared either (no finite score is >
        // infinity, and poisson_score is always finite). Either one alone would already
        // reject everything; setting both makes that explicit.
        let pvalue_filtered = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.0,
                min_region_score: f64::INFINITY,
                ..SearchFilters::default()
            },
        )?;
        assert!(pvalue_filtered.is_empty(), "capping both scopes at their extremes rejects all");

        Ok(())
    }

    /// The query p-value and region score combine with OR: a match is kept if either one
    /// passes. BCL2/CED9 at k=15 is a weak whole-query match (p ~ 0.99) carrying one strong
    /// region (score ~ 3.16, i.e. p ~ 0.0007). Capping only the query scope must not discard
    /// it, and capping only the region scope must keep it.
    #[test]
    fn test_pvalue_scopes_combine_with_or() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 15, 1, "hp_lehninger2", true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(target_index)?;

        let query_index_path = temp_dir.path().join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, 15, 1, "hp_lehninger2", true)?;
        query_index.process_fasta(TEST_CED9_FASTA, 0, DEFAULT_BATCH_SIZE)?;
        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        let find_bcl2 =
            |results: &[SearchResult]| results.iter().any(|r| r.target_name.contains("BCL2_HUMAN"));

        // -log10(0.05): the score-scale equivalent of the CLI's default 0.05 cutoff.
        let default_min_region_score = -0.05_f64.log10();

        // Whole-query scope alone rejects it: its query p-value is ~0.99, far above 0.05.
        // min_region_score: INFINITY disables the region scope entirely (no finite score
        // clears it).
        let query_only = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.05,
                min_region_score: f64::INFINITY,
                ..SearchFilters::default()
            },
        )?;
        assert!(!find_bcl2(&query_only), "query scope alone should reject the diluted BCL2 match");

        // Region scope alone keeps it, on the strength of the one significant region.
        // max_query_pvalue: 0.0 disables the query scope entirely (query_pvalue < 0.0 is
        // impossible).
        let region_only = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.0,
                min_region_score: default_min_region_score,
                ..SearchFilters::default()
            },
        )?;
        assert!(find_bcl2(&region_only), "region scope should keep the sub-protein domain hit");

        // Both at the CLI default: OR means the region rescues it.
        let both = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.05,
                min_region_score: default_min_region_score,
                ..SearchFilters::default()
            },
        )?;
        assert!(find_bcl2(&both), "default OR semantics should surface the hit");

        Ok(())
    }

    /// `ProteinSearcher::new()` (used by most tests) reads the cache the index just finalized,
    /// so `search_one()` always takes its "Path 1" branch. Only `ProteinSearcher::load()` (the
    /// fast path used by the real CLI, backed by an on-demand `sig_cache`) exercises Path 3
    /// (first RocksDB load of a target) and Path 2 (sig_cache hit on a later query that shares
    /// that target). Two queries against the same small target database, both known to match
    /// entries in it, should touch at least one common target twice.
    #[rstest]
    fn test_search_one_sig_cache_paths(
        bcl2_sketch_k12: ProteinSketch,
        ced9_sketch_k12: ProteinSketch,
    ) -> Result<()> {
        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 12, 1, "hp_lehninger2", false)?;
        target_index.process_fasta(TEST_FASTA_GZ, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;
        target_index.save_state()?;
        drop(target_index);

        // ProteinSearcher::load() takes the fast path: no signatures loaded up front, so
        // search_one() must go through sig_cache (Path 2/3), never Path 1.
        let searcher = ProteinSearcher::load(&target_index_path)?;
        let query_signatures = vec![bcl2_sketch_k12, ced9_sketch_k12];
        let results = searcher.search(&query_signatures, &SearchFilters::default())?;
        assert!(!results.is_empty(), "BCL2/CED9 queries should match this BCL2-family database");

        // At least one target must be shared by both queries' candidate sets - otherwise Path 2
        // (sig_cache hit) is never exercised, and this test would pass without covering it.
        let bcl2_targets: std::collections::HashSet<_> = results
            .iter()
            .filter(|r| r.query_name.contains("BCL2"))
            .map(|r| r.target_md5.clone())
            .collect();
        let ced9_targets: std::collections::HashSet<_> = results
            .iter()
            .filter(|r| r.query_name.contains("CED9"))
            .map(|r| r.target_md5.clone())
            .collect();
        assert!(
            bcl2_targets.intersection(&ced9_targets).next().is_some(),
            "BCL2 and CED9 queries should share at least one matched target in this database"
        );

        Ok(())
    }

    /// Tests if the correct k-mer overlap region for a query-target pair is found
    // # BCL2 & Ced9 `hp` k-mer match via **`pphhphhphhhhhphhhhh`, yay!**
    //
    // This is awesome because it's evidence that the hp k-mer method can work! The k-mer sizes
    // that this works for are k ≤ 19, since this region is of length 19
    //
    // → Need to figure out how we can auto-detect the k-mer size necessary for a query protein.
    //
    // ```
    // Ced9 pr: …RTVGNAQTD**QCPMSYGRLIGLISFGGFV**AAKMMESVE…
    // Ced9 hp: …pphhphppp**pphhphhphhhhhphhhhh**hhphhpphp…
    //                    |||||||||||||||||||
    // BCL2 hp: …hhphhpphh**pphhphhphhhhhphhhhh**phpphppph…
    // BCL2 pr: …FATVVEELF**RDGVNWGRIVAFFEFGGVM**CVESVNREM…
    //
    // - RDG starts at 138-157
    // ```
    /// Helper function to assert that a matched region matches expected values
    ///
    /// WHY: This helper function eliminates code duplication in tests and makes assertions
    /// more readable. It's idiomatic Rust to extract common assertion logic into helper functions.
    fn assert_region_matches(region: &MatchedRegion, expected: &ExpectedMatchedRegion) {
        assert_eq!(region.subseq, expected.subseq);
        assert_eq!(region.moltype_seq, expected.moltype_seq);
        assert_eq!(region.target_subseq, expected.target_subseq);
        assert_eq!(region.start, expected.start);
        assert_eq!(region.end, expected.end);
        assert_eq!(region.target_start, expected.target_start);
        assert_eq!(region.target_end, expected.target_end);
    }

    #[rstest]
    #[case::k12(ced9_sketch_k12(), bcl2_sketch_k12(), MATCHED_REGIONS_K12)]
    #[case::k15(ced9_sketch_k15(), bcl2_sketch_k15(), MATCHED_REGIONS_K15)]
    fn test_find_matched_regions(
        #[case] query_sketch: ProteinSketch,
        #[case] target_sketch: ProteinSketch,
        #[case] expected: ExpectedMatchedRegions,
    ) -> Result<()> {
        // Calculate intersection for find_matched_regions using the new intersect() method
        // WHY: We use a standalone function that doesn't require a searcher/index, making tests
        // simpler and more focused. This is idiomatic Rust - functions that don't need state
        // should be standalone.
        let intersection = query_sketch.intersect(&target_sketch);

        // Find matched regions using the standalone function
        let matched_regions = find_matched_regions(&query_sketch, &target_sketch, &intersection);

        // Verify total number of regions
        assert_eq!(
            matched_regions.len(),
            expected.total_regions,
            "Should find {} matched regions",
            expected.total_regions
        );

        // Verify first region
        assert_region_matches(&matched_regions[0], &expected.first);

        // Verify last region
        assert_region_matches(&matched_regions[matched_regions.len() - 1], &expected.last);

        // Verify landmark region if specified
        if let Some(landmark) = &expected.landmark {
            let found = matched_regions
                .iter()
                .find(|r| r.subseq == landmark.subseq)
                .expect("Should find landmark region");
            assert_region_matches(found, landmark);
        }

        Ok(())
    }

    /// Two sequences whose HP encodings agree on 25 positions except one flip in the
    /// middle. The pattern has no repeated 8-mer, so at k=8 the only shared k-mers sit on
    /// the main diagonal; residues cycle so the raw sequences differ too.
    fn one_flip_pair() -> (ProteinSketch, ProteinSketch) {
        // h-class residues cycle through AFILMV, p-class through DEKNQR (Lehninger).
        let pattern = "hpphhpphphphhppphhhppphhp";
        let flip_at = 12;
        let build = |name: &str, flip: bool| {
            let (h, p) = ("AFILMV".as_bytes(), "DEKNQR".as_bytes());
            let mut hi = 0usize;
            let mut pi = 0usize;
            let seq: String = pattern
                .bytes()
                .enumerate()
                .map(|(i, c)| {
                    let is_h = (c == b'h') != (flip && i == flip_at);
                    let r = if is_h {
                        hi += 1;
                        h[hi % h.len()]
                    } else {
                        pi += 1;
                        p[pi % p.len()]
                    };
                    r as char
                })
                .collect();
            ProteinSketch::from_protein_sequence(name, &seq, 8, 1, "hp").unwrap()
        };
        (build("q", false), build("t", true))
    }

    #[test]
    fn test_karlin_altschul_lambda() {
        // u e^x + (1-u) e^(-2x) = 1 at u = 0.5: e^x = 1.618..., x = ln(golden ratio).
        let lam = karlin_altschul_lambda(0.5, 2.0);
        assert_relative_eq!(lam, ((1.0 + 5f64.sqrt()) / 2.0).ln(), epsilon = 1e-9);
        // No positive root once agreement is expected: u >= penalty / (1 + penalty).
        assert_eq!(karlin_altschul_lambda(2.0 / 3.0, 2.0), 0.0);
        assert_eq!(karlin_altschul_lambda(0.9, 2.0), 0.0);
        // Rarer agreement, larger lambda.
        assert!(karlin_altschul_lambda(0.3, 2.0) > lam);
        // The root satisfies the equation for a 20-letter-like composition too.
        let u = 0.06;
        let l = karlin_altschul_lambda(u, 1.0);
        assert_relative_eq!(u * l.exp() + (1.0 - u) * (-l).exp(), 1.0, epsilon = 1e-9);
    }

    #[test]
    fn test_karlin_altschul_sum_p() {
        // r = 1 is the plain exponential tail.
        assert_relative_eq!(karlin_altschul_sum_p(3.0, 1), (-3.0f64).exp(), epsilon = 1e-12);
        // r = 2: e^-t t / (2! 1!).
        assert_relative_eq!(
            karlin_altschul_sum_p(4.0, 2),
            (-4.0f64).exp() * 4.0 / 2.0,
            epsilon = 1e-12
        );
        assert_eq!(karlin_altschul_sum_p(0.0, 2), 1.0);
        assert_eq!(karlin_altschul_sum_p(-1.0, 1), 1.0);
        assert!(karlin_altschul_sum_p(50.0, 3) < karlin_altschul_sum_p(20.0, 3));
        assert_relative_eq!(libm_lgamma(5.0), (24.0f64).ln(), epsilon = 1e-9);
    }

    #[test]
    fn test_chain_regions_joins_one_diagonal_within_gap() {
        // Two exact seeds on one diagonal with a bad patch between them that a strict
        // extension will not cross; chaining joins them, a gap cap below the patch does not.
        let (q, t) = one_flip_pair();
        let intersection = q.intersect(&t);
        let exact = find_matched_regions(&q, &t, &intersection);
        assert_eq!(exact.len(), 2);
        let strict = ExtensionParams {
            mismatch_penalty: 9.0,
            xdrop: 8.0,
            ka_k: 0.03,
            ka_lambda_scale: 1.0,
            chain_max_gap: 5,
            chain_max_shift: 0,
        };
        let ext = extend_regions(exact.clone(), &q, &t, strict);
        assert_eq!(ext.len(), 2);
        let (qe, te) = (q.get_moltype_sequence().unwrap(), t.get_moltype_sequence().unwrap());
        let chained = chain_regions(
            ext.clone(),
            qe.as_bytes(),
            te.as_bytes(),
            q.get_raw_sequence(),
            t.get_raw_sequence(),
            strict,
            25.0,
            25.0,
            100.0,
        );
        assert_eq!(chained.len(), 1, "{chained:?}");
        let c = &chained[0];
        assert_eq!((c.start, c.end, c.target_start, c.target_end), (0, 25, 0, 25));
        assert_eq!(c.n_chained, 2);
        assert_eq!(c.n_shared, 10);
        assert_eq!(c.n_mismatches, 1);
        assert_eq!(c.length, 25);
        // The chain's lambda comes from the composition of the span it covers, here the
        // whole 25 residues of both sequences: 12 h of 25 against 11 h of 25.
        assert_relative_eq!(c.ka_u, 0.48 * 0.44 + 0.52 * 0.56, epsilon = 1e-12);
        assert_relative_eq!(
            c.ka_lambda,
            karlin_altschul_lambda(c.ka_u, strict.mismatch_penalty),
            epsilon = 1e-12
        );
        assert!(c.ka_lambda > 0.0, "assessable, so the E-value is a number");
        assert!(c.evalue.is_finite());
        let tight = ExtensionParams { chain_max_gap: 0, ..strict };
        let kept = chain_regions(
            ext.clone(),
            qe.as_bytes(),
            te.as_bytes(),
            q.get_raw_sequence(),
            t.get_raw_sequence(),
            tight,
            25.0,
            25.0,
            100.0,
        );
        assert_eq!(kept.len(), 2);
        assert!(kept.iter().all(|r| r.n_chained == 1));
    }

    /// Human BNIP3 residues 69-108 and human BNIP3L residues 84-123 (UniProt Q12983 and
    /// O60238, both in `tests/testdata/index/bcl2_first25_...fasta`). Two ordinary
    /// apoptosis proteins, each carrying one polar-rich stretch.
    const BNIP3_POLAR_STRETCH: &[u8] = b"RSQTPQDTNRASETDTHSIGEKNSSQSEEDDIERRKEVES";
    const BNIP3L_POLAR_STRETCH: &[u8] = b"QSSSRGSSHCDSPSPQEDGQIMFDVEMHTSRDHSSQSEEE";

    #[test]
    fn test_region_lambda_is_zero_for_a_polar_stretch_in_an_ordinary_protein() {
        // Taken from the whole proteins, the compositions are ordinary and every region
        // between them is scored as if its matches were evidence.
        let (q_whole, t_whole) = (b"hpphhpphphphhppphhhppphhp", b"hphhpphhpphphhpphhpphhphp");
        let u_whole = region_match_probability(q_whole, t_whole);
        assert!(u_whole < 2.0 / 3.0, "the two proteins look ordinary: u = {u_whole}");
        assert!(karlin_altschul_lambda(u_whole, 2.0) > 0.0);

        // Taken from the two stretches themselves, 85% and 78% polar, they match at more
        // than two positions in three by composition alone. At penalty 2 the boundary is
        // 2/3, there is no positive root, and the region is not assessable.
        let q = ProteinSketch::from_protein_sequence(
            "bnip3",
            std::str::from_utf8(BNIP3_POLAR_STRETCH).unwrap(),
            12,
            1,
            "hp",
        )
        .unwrap();
        let t = ProteinSketch::from_protein_sequence(
            "bnip3l",
            std::str::from_utf8(BNIP3L_POLAR_STRETCH).unwrap(),
            12,
            1,
            "hp",
        )
        .unwrap();
        let (q_enc, t_enc) = (q.get_moltype_sequence().unwrap(), t.get_moltype_sequence().unwrap());
        let u_region = region_match_probability(q_enc.as_bytes(), t_enc.as_bytes());
        assert_relative_eq!(u_region, 0.6925, epsilon = 1e-4);
        assert!(u_region > 2.0 / 3.0, "the two stretches match for free: u = {u_region}");
        assert_eq!(karlin_altschul_lambda(u_region, 2.0), 0.0);
    }

    #[test]
    fn test_region_span_clamps_to_the_sequence() {
        let seq = b"hpphhpph";
        assert_eq!(region_span(seq, 2, 5), b"phh");
        assert_eq!(region_span(seq, 0, 8), seq);
        // Past the end, and inverted: clamped, not a panic.
        assert_eq!(region_span(seq, 6, 99), b"ph");
        assert_eq!(region_span(seq, 99, 120), b"");
        assert_eq!(region_span(seq, 5, 2), b"");
    }

    #[test]
    fn test_class_composition_match_probability() {
        let p = class_composition(b"hhpp");
        let q = class_composition(b"hhhp");
        assert_relative_eq!(match_probability(&p, &q), 0.5 * 0.75 + 0.5 * 0.25, epsilon = 1e-12);
    }

    #[test]
    fn test_extend_regions_bridges_one_flip() {
        let (q, t) = one_flip_pair();
        let intersection = q.intersect(&t);
        let exact = find_matched_regions(&q, &t, &intersection);
        // Exact runs stop at the flip: positions 0..12 and 13..25, five 8-mers each.
        assert_eq!(exact.len(), 2, "{exact:?}");
        assert_eq!((exact[0].start, exact[0].end, exact[0].n_shared), (0, 12, 5));
        assert_eq!((exact[1].start, exact[1].end, exact[1].n_shared), (13, 25, 5));
        assert!(exact.iter().all(|r| r.n_mismatches == 0));

        let params = ExtensionParams {
            mismatch_penalty: 2.0,
            xdrop: 8.0,
            ka_k: 0.1,
            ka_lambda_scale: 1.0,
            chain_max_gap: 0,
            chain_max_shift: 0,
        };
        let extended = extend_regions(exact.clone(), &q, &t, params);
        assert_eq!(extended.len(), 1, "{extended:?}");
        let r = &extended[0];
        assert_eq!((r.start, r.end, r.target_start, r.target_end), (0, 25, 0, 25));
        assert_eq!(r.length, 25);
        assert_eq!(r.n_shared, 10, "shared k-mers are summed from the seeds, not recomputed");
        assert_eq!(r.n_mismatches, 1);
        assert_eq!(r.subseq.len(), 25);
        assert_eq!(r.target_subseq.len(), 25);
        assert_eq!(r.moltype_seq, t.get_moltype_sequence().unwrap());

        // A penalty larger than the X-drop cannot cross the flip: the two seeds stay apart,
        // and nothing else changes about them.
        let strict = ExtensionParams {
            mismatch_penalty: 9.0,
            xdrop: 8.0,
            ka_k: 0.1,
            ka_lambda_scale: 1.0,
            chain_max_gap: 0,
            chain_max_shift: 0,
        };
        let kept = extend_regions(exact.clone(), &q, &t, strict);
        assert_eq!(kept.len(), 2);
        for (a, b) in kept.iter().zip(&exact) {
            assert_eq!(
                (a.start, a.end, a.n_shared, a.n_mismatches),
                (b.start, b.end, b.n_shared, 0)
            );
        }

        // Penalty 0 is "off" and returns the regions untouched.
        let off = ExtensionParams {
            mismatch_penalty: 0.0,
            xdrop: 8.0,
            ka_k: 0.1,
            ka_lambda_scale: 1.0,
            chain_max_gap: 0,
            chain_max_shift: 0,
        };
        let same = extend_regions(exact.clone(), &q, &t, off);
        assert_eq!(same.len(), exact.len());
    }

    #[test]
    fn test_extend_regions_contains_seeds_on_real_pair() -> Result<()> {
        // CED9 vs BCL2 at hp k=12: every extended region must contain the seed it grew from,
        // stay on its diagonal, and count its mismatches correctly.
        let (ced9_name, ced9_sequence) = read_first_fasta_record(TEST_CED9_FASTA)?;
        let (bcl2_name, bcl2_sequence) = read_first_fasta_record(TEST_BLC2_FASTA)?;
        let q = ProteinSketch::from_protein_sequence(&ced9_name, &ced9_sequence, 12, 1, "hp")?;
        let t = ProteinSketch::from_protein_sequence(&bcl2_name, &bcl2_sequence, 12, 1, "hp")?;
        let intersection = q.intersect(&t);
        let exact = find_matched_regions(&q, &t, &intersection);
        assert!(!exact.is_empty());
        let params = ExtensionParams {
            mismatch_penalty: 2.0,
            xdrop: 8.0,
            ka_k: 0.1,
            ka_lambda_scale: 1.0,
            chain_max_gap: 0,
            chain_max_shift: 0,
        };
        let extended = extend_regions(exact.clone(), &q, &t, params);
        assert!(!extended.is_empty());
        assert!(extended.len() <= exact.len());
        let (qe, te) = (
            q.get_moltype_sequence().unwrap().as_bytes(),
            t.get_moltype_sequence().unwrap().as_bytes(),
        );
        let mut seeds_covered = 0;
        for r in &extended {
            assert_eq!(
                r.target_start as i64 - r.start as i64,
                r.target_start as i64 - r.start as i64
            );
            assert_eq!(r.length, r.end - r.start);
            assert_eq!(r.target_end - r.target_start, r.length);
            let recount = qe[r.start as usize..r.end as usize]
                .iter()
                .zip(&te[r.target_start as usize..r.target_end as usize])
                .filter(|(a, b)| a != b)
                .count() as u32;
            assert_eq!(r.n_mismatches, recount);
            let inside: Vec<_> = exact
                .iter()
                .filter(|s| {
                    s.target_start as i64 - s.start as i64 == r.target_start as i64 - r.start as i64
                        && s.start >= r.start
                        && s.end <= r.end
                })
                .collect();
            assert!(!inside.is_empty(), "extended region {r:?} contains no seed");
            assert_eq!(r.n_shared, inside.iter().map(|s| s.n_shared).sum::<u32>());
            seeds_covered += inside.len();
        }
        assert_eq!(seeds_covered, exact.len(), "every seed lands in exactly one extended region");
        let total_growth: u32 = extended.iter().map(|r| r.length).sum::<u32>();
        let total_seed: u32 = exact.iter().map(|r| r.length).sum::<u32>();
        assert!(total_growth >= total_seed);
        Ok(())
    }

    #[test]
    fn test_find_matched_regions_multiple() -> Result<()> {
        // 14 is the minimum k-mersize that finds multiple match regions from Delilah's analyses
        let ksize = 12;
        let scaled = 1;
        let moltype = "hp_lehninger2";

        // Read CED9 sequence from FASTA file
        let (ced9_name, ced9_sequence) = read_first_fasta_record(TEST_CED9_FASTA)?;

        // Read BCL2 sequence from FASTA file
        let (bcl2_name, bcl2_sequence) = read_first_fasta_record(TEST_BLC2_FASTA)?;

        // Create sketches using from_protein_sequence - this now handles everything:
        // minhash, kmer_infos, raw sequence, and encoded sequence storage
        // WHY: The enhanced from_protein_sequence method does all the heavy lifting,
        // making tests simple and avoiding boilerplate. This is idiomatic Rust - make
        // the common case easy by having the method do what users almost always need.
        let query_sketch = ProteinSketch::from_protein_sequence(
            &ced9_name,
            &ced9_sequence,
            ksize,
            scaled,
            moltype,
        )?;

        let target_sketch = ProteinSketch::from_protein_sequence(
            &bcl2_name,
            &bcl2_sequence,
            ksize,
            scaled,
            moltype,
        )?;

        // Calculate intersection for find_matched_regions using the new intersect() method
        // WHY: We use a standalone function that doesn't require a searcher/index, making tests
        // simpler and more focused. This is idiomatic Rust - functions that don't need state
        // should be standalone.
        let intersection = query_sketch.intersect(&target_sketch);

        // Find matched regions using the standalone function
        let matched_regions = find_matched_regions(&query_sketch, &target_sketch, &intersection);

        // Verify we found exactly 3 matches
        // assert_eq!(matched_regions.len(), 3, "Should find exactly 3 matches");
        assert_eq!(matched_regions.len(), 13, "Should find exactly 13 matches");

        let first_match = &matched_regions[0];
        // Positions 87-99 in CED9 (query) and 103-115 in BCL2 (target)
        assert_eq!(first_match.subseq, "FTHRIRQNGMEW");
        assert_eq!(first_match.moltype_seq, "hppphppphhph");
        assert_eq!(first_match.target_subseq, "FSRRYRRDFAEM");
        assert_eq!(first_match.start, 87, "First match query start should be 87");
        assert_eq!(first_match.end, 99, "First match query end should be 99");
        assert_eq!(first_match.target_start, 103, "First match target start should be 103");
        assert_eq!(first_match.target_end, 115, "First match target end should be 115");

        let last_match = &matched_regions[matched_regions.len() - 1];
        // Positions 267-280 in CED9 (query) and 200-213 in BCL2 (target)
        assert_eq!(last_match.subseq, "GVVVCGRMMFSLK");
        assert_eq!(last_match.moltype_seq, "hhhhphphhhphp");
        assert_eq!(last_match.target_subseq, "LYGPSMRPLFDFS");
        assert_eq!(last_match.start, 267, "Last match query start should be 267");
        assert_eq!(last_match.end, 280, "Last match query end should be 280");
        assert_eq!(last_match.target_start, 200, "Last match target start should be 200");
        assert_eq!(last_match.target_end, 213, "Last match target end should be 213");

        // Verify the expected match region at positions 162-181 in CED9 (query) and 138-157 in BCL2 (target)
        // Query subsequence: "QCPMSYGRLIGLISFGGFV"
        // Target subsequence: "RDGVNWGRIVAFFEFGGVM"
        // Moltype sequence: "pphhphhphhhhhphhhhh"
        // WHY: Since regions are now sorted by position (not length), we need to find
        // the specific region by its subsequence rather than assuming it's at index 0.
        let largest_match = matched_regions
            .iter()
            .find(|r| r.subseq == "QCPMSYGRLIGLISFGGFV")
            .expect("Should find the expected match region");
        assert_eq!(largest_match.subseq, "QCPMSYGRLIGLISFGGFV");
        assert_eq!(largest_match.moltype_seq, "pphhphhphhhhhphhhhh");
        assert_eq!(largest_match.target_subseq, "RDGVNWGRIVAFFEFGGVM");
        assert_eq!(largest_match.start, 162, "Largest match query start should be 162");
        assert_eq!(largest_match.end, 181, "Largest match query end should be 181");
        assert_eq!(largest_match.target_start, 138, "Largest match target start should be 138");
        assert_eq!(largest_match.target_end, 157, "Largest match target end should be 157");

        Ok(())
    }

    /// Test TF-IDF calculation
    #[test]
    fn test_tfidf_calculation() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path();

        // Create target FASTA with multiple sequences
        let target_fasta = temp_path.join("target.fasta");
        std::fs::write(&target_fasta, ">seq1\nATCGATCGATCGATCG\n>seq2\nGCTAGCTAGCTAGCTA")?;

        // Create target index
        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp_lehninger2", false)?;

        target_index.process_fasta(&target_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        // Create searcher
        let searcher = ProteinSearcher::new(target_index)?;

        // Create query FASTA
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">query\nATCGATCGATCGATCG")?;

        let query_index =
            ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp_lehninger2", false)?;

        query_index.process_fasta(&query_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        // Calculate TF-IDF
        let tfidf = searcher.calculate_tfidf(&query_signatures[0]);
        assert!(tfidf >= 0.0, "TF-IDF should be non-negative");

        Ok(())
    }

    /// Test search result structure matches expected format
    #[test]
    fn test_search_result_structure() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path();

        // Create test data
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">test_query\nATCGATCGATCGATCG")?;

        let target_fasta = temp_path.join("target.fasta");
        std::fs::write(&target_fasta, ">test_target\nATCGATCGATCGATCG")?;

        // Create indices
        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp_lehninger2", false)?;
        target_index.process_fasta(&target_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        let searcher = ProteinSearcher::new(target_index)?;

        let query_index =
            ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp_lehninger2", false)?;
        query_index.process_fasta(&query_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        let results = searcher.search(&query_signatures, &SearchFilters::default())?;

        if !results.is_empty() {
            let result = &results[0];

            // Test that all required fields are present and have reasonable values
            assert!(!result.query_name.is_empty());
            assert!(!result.query_md5.is_empty());
            assert!(!result.target_name.is_empty());
            assert!(!result.target_md5.is_empty());
            assert!(!result.moltype.is_empty());

            assert!(result.containment >= 0.0 && result.containment <= 1.0);
            assert!(result.jaccard >= 0.0 && result.jaccard <= 1.0);
            assert!(result.max_containment >= 0.0 && result.max_containment <= 1.0);
            assert!(result.n_intersecting_hashes > 0);
            assert!(result.ksize > 0);
            assert!(result.scaled > 0);

            // Test abundance statistics
            assert!(result.average_abund >= 0.0);
            assert!(result.median_abund >= 0.0);
            assert!(result.std_abund >= 0.0);

            assert!(
                result.containment_target_in_query >= 0.0
                    && result.containment_target_in_query <= 1.0
            );
            assert!(result.f_weighted_target_in_query >= 0.0);
        }

        Ok(())
    }

    /// Test that search results are sorted by containment score
    #[test]
    fn test_search_results_sorted() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path();

        // Create target with multiple sequences of different similarity
        let target_fasta = temp_path.join("target.fasta");
        std::fs::write(&target_fasta, ">exact_match\nATCGATCGATCGATCG\n>partial_match\nATCGATCGATCGATCA\n>no_match\nGGGGGGGGGGGGGGGG")?;

        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp_lehninger2", false)?;
        target_index.process_fasta(&target_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        let searcher = ProteinSearcher::new(target_index)?;

        // Create query
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">query\nATCGATCGATCGATCG")?;

        let query_index =
            ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp_lehninger2", false)?;
        query_index.process_fasta(&query_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        let results = searcher.search(&query_signatures, &SearchFilters::default())?;

        // Results should be sorted by containment (descending)
        for i in 1..results.len() {
            assert!(
                results[i - 1].containment >= results[i].containment,
                "Results should be sorted by containment score"
            );
        }

        Ok(())
    }

    /// Test search stats calculation
    #[test]
    fn test_search_stats_calculation() {
        // This would need a proper test with actual signatures
        // For now, just test the structure
        let stats = SearchStats {
            total_signatures: 100,
            idf: HashMap::new(),
            kmer_frequencies: HashMap::new(),
        };

        assert_eq!(stats.total_signatures, 100);
    }

    /// Test basic TF-IDF calculation structure
    #[test]
    fn test_tfidf_calculation_structure() -> Result<()> {
        // Create a temporary directory for the test
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path().join("test.db");

        // Create a proper index with a database path
        let index = ProteomeIndex::new(&temp_path, 10, 5, "hp_lehninger2", false)?;

        let query = ProteinSketch::new("test", 10, 5, "hp_lehninger2")?;
        let stats = SearchStats {
            total_signatures: 100,
            idf: HashMap::new(),
            kmer_frequencies: HashMap::new(),
        };

        let searcher = ProteinSearcher {
            index,
            stats,
            target_list: Vec::new(),
            inverted_index: HashMap::new(),
            sig_cache: DashMap::new(),
            aliases: HashMap::new(),
            query_kmer_frequencies: None,
            total_queries: 0,
            db_n_kmers: 0,
            extension: None,
        };

        let tfidf = searcher.calculate_tfidf(&query);
        assert!(tfidf >= 0.0);

        Ok(())
    }

    // From Sourmash values:
    // $ sourmash sig overlap -k 12 ced9.fasta.hp.k12-15.scaled1.sig.zip bcl2.fasta.hp.k12-15.scaled1.sig.zip

    // == This is sourmash version 4.9.4. ==
    // == Please cite Irber et. al (2024), doi:10.21105/joss.06830. ==

    // loaded one signature each from ced9.fasta.hp.k12-15.scaled1.sig.zip and bcl2.fasta.hp.k12-15.scaled1.sig.zip
    // size_estimate_inaccurate: False
    // first signature:
    //   signature filename: ced9.fasta.hp.k12-15.scaled1.sig.zip
    //   signature name: sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
    //   source filename: ced9.fasta
    //   md5: 5baa6059c3306c2b6a500abd929d539d
    //   k=12 molecule=hp num=0 scaled=1 track_abundance=False
    //   size: 264
    //   sum hashes: 264
    //   signature license: CC0

    // second signature:
    //   signature filename: bcl2.fasta.hp.k12-15.scaled1.sig.zip
    //   signature name: sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=BCL2 PE=1 SV=2
    //   source filename: bcl2.fasta
    //   md5: b6da406482d741a9d49f406684c17e17
    //   k=12 molecule=hp num=0 scaled=1 track_abundance=False
    //   size: 220
    //   sum hashes: 220
    //   signature license: CC0

    // --- Similarity measures ---
    // jaccard similarity:          0.05217
    // first contained in second:   0.09091 (cANI: 0.81887)
    // second contained in first:   0.10909 (cANI: 0.83141)
    // average containment ANI:     0.82514

    // --- Hash overlap summary ---
    // number of hashes in first:   264
    // number of hashes in second:  220

    // number of hashes in common:  24
    // only in first:               240
    // only in second:              196
    // total (union):               460
    const BCL2_CED9_K12: ExpectedSimilarity = ExpectedSimilarity {
        ksize: 12,
        n_intersecting_hashes: 24,
        containment: 0.09091,
        jaccard: 0.05217,
        max_containment: 0.10909,
        containment_target_in_query: 0.10909,
        average_abund: 1.0625,
        median_abund: 1.0,
        std_abund: 0.099_913_156_735_681_66,
        matched_regions_count: 13,
    };

    // From Sourmash values:
    // $ sourmash sig overlap -k 15 ced9.fasta.hp.k12-15.scaled1.sig.zip bcl2.fasta.hp.k12-15.scaled1.sig.zip
    // == This is sourmash version 4.9.4. ==
    // == Please cite Irber et. al (2024), doi:10.21105/joss.06830. ==

    // loaded one signature each from ced9.fasta.hp.k12-15.scaled1.sig.zip and bcl2.fasta.hp.k12-15.scaled1.sig.zip
    // size_estimate_inaccurate: False
    // first signature:
    // signature filename: ced9.fasta.hp.k12-15.scaled1.sig.zip
    // signature name: sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
    // source filename: ced9.fasta
    // md5: 61094124a51b6d4802c37cd7bb43fad5
    // k=15 molecule=hp num=0 scaled=1 track_abundance=False
    // size: 266
    // sum hashes: 266
    // signature license: CC0

    // second signature:
    // signature filename: bcl2.fasta.hp.k12-15.scaled1.sig.zip
    // signature name: sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=BCL2 PE=1 SV=2
    // source filename: bcl2.fasta
    // md5: 26546d1eea2143522424bc59a5ded0f5
    // k=15 molecule=hp num=0 scaled=1 track_abundance=False
    // size: 225
    // sum hashes: 225
    // signature license: CC0

    // --- Similarity measures ---
    // jaccard similarity:          0.01029
    // first contained in second:   0.01880 (cANI: 0.76725)
    // second contained in first:   0.02222 (cANI: 0.77586)
    // average containment ANI:     0.77156

    // --- Hash overlap summary ---
    // number of hashes in first:   266
    // number of hashes in second:  225

    // number of hashes in common:  5
    // only in first:               261
    // only in second:              220
    // total (union):               486
    const BCL2_CED9_K15: ExpectedSimilarity = ExpectedSimilarity {
        ksize: 15,
        n_intersecting_hashes: 5,
        containment: 0.018796992481203006,
        jaccard: 0.0102880658436214,
        max_containment: 0.022222222222222223,
        containment_target_in_query: 0.022222222222222223,
        average_abund: 1.0,
        median_abund: 1.0,
        std_abund: 0.0,
        matched_regions_count: 1,
    };

    #[rstest]
    #[case::k12(ced9_sketch_k12(), bcl2_sketch_k12(), BCL2_CED9_K12)]
    #[case::k15(ced9_sketch_k15(), bcl2_sketch_k15(), BCL2_CED9_K15)]
    fn test_calculate_similarity_bcl2_ced9(
        #[case] query_sketch: ProteinSketch,
        #[case] target_sketch: ProteinSketch,
        #[case] expected: ExpectedSimilarity,
    ) -> Result<()> {
        // Use the standalone calculate_similarity function for simple 1v1 comparisons
        // WHY: This is much simpler than creating a searcher and index just to test similarity.
        // The standalone function is designed exactly for this use case - 1v1 comparisons without
        // database context. TF-IDF and overlap probability will be defaults (0.0 and 1.0), which
        // is correct for 1v1 comparisons where these metrics are meaningless.
        let result = calculate_similarity(&query_sketch, &target_sketch)
            .expect("Should find similarity between CED9 and BCL2");

        // Assertions are now readable
        assert_eq!(result.n_intersecting_hashes, expected.n_intersecting_hashes);
        assert_relative_eq!(result.containment, expected.containment, epsilon = 1e-5);
        assert_relative_eq!(result.jaccard, expected.jaccard, epsilon = 1e-5);
        assert_relative_eq!(result.max_containment, expected.max_containment, epsilon = 1e-5);
        assert_relative_eq!(result.average_abund, expected.average_abund, epsilon = 1e-5);
        assert_eq!(result.matched_regions.len(), expected.matched_regions_count);

        // Verify database-specific metrics are defaults for 1v1 comparisons
        assert_eq!(result.query_tfidf, 0.0, "query_tfidf should be 0.0 for 1v1 comparisons");
        assert_eq!(
            result.mean_matched_kmer_freq, 0.0,
            "mean_matched_kmer_freq should be 0.0 for 1v1 comparisons"
        );
        assert_eq!(
            result.sum_matched_kmer_freq, 0.0,
            "sum_matched_kmer_freq should be 0.0 for 1v1 comparisons"
        );
        assert_eq!(
            result.query_expected_shared_kmers, 0.0,
            "query_expected_shared_kmers should be 0.0 for 1v1 comparisons"
        );
        assert_eq!(
            result.query_enrichment, 0.0,
            "query_enrichment should be 0.0 for 1v1 comparisons"
        );
        assert_eq!(
            result.joint_kmer_freq, 0.0,
            "joint_kmer_freq should be 0.0 without query frequencies"
        );

        Ok(())
    }

    /// Test the public search() method with a real database (multiple signatures)
    ///
    /// This test uses:
    /// - Target database: bcl2_first25 (multiple BCL2 family proteins) + bcl2.fasta (single BCL2)
    /// - Query: ced9.fasta (single CED9 sequence)
    ///
    /// This properly tests TF-IDF and overlap probability metrics since we have multiple
    /// signatures in the database where some k-mers are common and others are rare.
    #[test]
    fn test_search_database_bcl2_ced9() -> Result<()> {
        let ksize = 12;
        let scaled = 1;
        let moltype = "hp_lehninger2";

        // Create temporary directory for the index
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path();

        // Create target index with both bcl2_first25 and bcl2.fasta
        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;

        // Process the bcl2_first25 FASTA (contains multiple BCL2 family proteins)
        // WHY: This gives us a database with multiple signatures, enabling meaningful
        // TF-IDF and overlap probability calculations. Some k-mers will appear in many
        // signatures (common) while others will appear in few (rare).
        // This fasta already includes the BCL2 sequence, so we don't need to add it again.
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;

        // Verify we have multiple signatures in the index
        let signature_count = target_index.signature_count();
        assert!(
            signature_count == 25,
            "Should have 26 signatures in the database, got {}",
            signature_count
        );

        // Create searcher from the index
        let searcher = ProteinSearcher::new(target_index)?;

        // Create query index from CED9 in a temporary directory
        // WHY: new_with_auto_filename creates the database next to the input file, which causes
        // RocksDB lock conflicts when tests run in parallel. Using a temporary directory ensures
        // each test run has its own isolated database path.
        let query_index_path = temp_path.join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, ksize, scaled, moltype, true)?;

        query_index.process_fasta(TEST_CED9_FASTA, 0, DEFAULT_BATCH_SIZE)?;

        // Get query signatures
        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        assert_eq!(query_signatures.len(), 1, "Should have exactly one query signature (CED9)");

        // Perform search using the public search() method
        let results = searcher.search(&query_signatures, &SearchFilters::default())?;

        // Should find at least one match (CED9 should match BCL2 and potentially other BCL2 family members)
        assert_eq!(
            results.len(),
            25,
            "Should find exactly 25 matches between CED9 and BCL2 database"
        );

        // Find the result for the canonical BCL2 sequence
        let bcl2_result = results
            .iter()
            .find(|r| r.target_name.contains("BCL2_HUMAN"))
            .expect("Should find a match with BCL2_HUMAN");

        // Verify basic metadata
        assert!(
            bcl2_result.query_name.contains("CED9_CAEEL"),
            "Query name should contain CED9_CAEEL"
        );
        assert!(
            bcl2_result.target_name.contains("BCL2_HUMAN"),
            "Target name should contain BCL2_HUMAN"
        );
        assert_eq!(bcl2_result.ksize, ksize);
        assert_eq!(bcl2_result.scaled, scaled);
        // Results carry the normalized moltype, not the spelling the index was created with.
        assert_eq!(bcl2_result.moltype, "hp_lehninger2");

        // Verify we have intersecting k-mers
        assert!(
            bcl2_result.n_intersecting_hashes > 0,
            "Should have intersecting k-mers between CED9 and BCL2"
        );

        // Verify similarity metrics are within valid ranges
        assert!(
            bcl2_result.containment > 0.0 && bcl2_result.containment <= 1.0,
            "Containment should be in [0, 1], got {}",
            bcl2_result.containment
        );
        assert!(
            bcl2_result.jaccard > 0.0 && bcl2_result.jaccard <= 1.0,
            "Jaccard should be in [0, 1], got {}",
            bcl2_result.jaccard
        );

        // Verify TF-IDF is meaningful (should not be 0 with multiple signatures)
        // WHY: With multiple signatures in the database, some k-mers will be rare and have
        // higher IDF values, making TF-IDF > 0. In a real database search, TF-IDF helps
        // identify matches based on rare, significant k-mers.
        // Compared approximately, not with ==: this is a sum over a HashMap keyed by k-mer
        // hash, so the summation order (and therefore the last bits) depends on the hash
        // values themselves.
        approx::assert_relative_eq!(bcl2_result.query_tfidf, 565.119680433367, epsilon = 1e-9);

        assert!(
            bcl2_result.mean_matched_kmer_freq > 0.0,
            "mean_matched_kmer_freq should be > 0.0, got {}",
            bcl2_result.mean_matched_kmer_freq
        );

        // joint_kmer_freq is 0.0 because set_query_frequencies() was not called
        assert_eq!(
            bcl2_result.joint_kmer_freq, 0.0,
            "joint_kmer_freq should be 0.0 without query frequencies, got {}",
            bcl2_result.joint_kmer_freq
        );

        // Verify matched regions are present
        assert!(
            !bcl2_result.matched_regions.is_empty(),
            "Should find matched regions between CED9 and BCL2"
        );

        // Verify results are sorted by containment (descending)
        // WHY: The search() method sorts results by containment score, with the best matches first.
        // This is a key feature of the search API - users expect results in order of relevance.
        for i in 1..results.len() {
            assert!(
                results[i - 1].containment >= results[i].containment,
                "Results should be sorted by containment (descending), but result {} has containment {} and result {} has containment {}",
                i - 1,
                results[i - 1].containment,
                i,
                results[i].containment
            );
        }

        Ok(())
    }

    /// Test joint_kmer_freq is computed correctly when set_query_frequencies() is called.
    ///
    /// Key invariant: when total_queries=1 and every query k-mer has frequency 1
    /// (i.e. the "query proteome" is a single sequence),
    ///   joint_kmer_freq = Σ (freq_q[h]/1) * (freq_t[h]/N_targets)
    ///               = Σ freq_t[h]/N_targets
    ///               = sum_matched_kmer_freq
    /// so the two metrics must be equal in this degenerate case.
    #[test]
    fn test_joint_kmer_freq_two_pass() -> Result<()> {
        let ksize = 12;
        let scaled = 1;
        let moltype = "hp_lehninger2";

        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;

        let mut searcher = ProteinSearcher::new(target_index)?;

        // Build a query sketch for CED9
        let mut query_sig = ProteinSketch::new("ced9_query", ksize, scaled, moltype)?;
        let ced9_seq = {
            let mut reader = needletail::parse_fastx_file(TEST_CED9_FASTA)
                .map_err(|e| anyhow::anyhow!("{}", e))?;
            let record = reader.next().unwrap().map_err(|e| anyhow::anyhow!("{}", e))?;
            std::str::from_utf8(&record.seq()).unwrap().to_uppercase()
        };
        query_sig.add_protein(&ced9_seq, true)?;

        // Simulate a "query proteome" of 1 sequence: every k-mer in ced9 has freq=1
        let qfreqs: HashMap<u64, usize> =
            query_sig.signature().minhash.mins().iter().map(|&h| (h, 1usize)).collect();
        let total_queries = 1;
        searcher.set_query_frequencies(qfreqs, total_queries);

        let results = searcher.search_one(&query_sig, &SearchFilters::default(), total_queries);

        // Find bcl2 result
        let bcl2_result = results
            .iter()
            .find(|r| r.target_name.contains("BCL2_HUMAN"))
            .expect("Should find BCL2_HUMAN in results");

        // With total_queries=1 and all query k-mer frequencies=1:
        //   joint_kmer_freq = sum_matched_kmer_freq (exact equality)
        // joint_kmer_freq must equal sum_matched_kmer_freq when total_queries=1, all query freqs=1
        assert_relative_eq!(
            bcl2_result.joint_kmer_freq,
            bcl2_result.sum_matched_kmer_freq,
            epsilon = 1e-12
        );

        // Also verify it's non-zero (there's an actual intersection)
        assert!(
            bcl2_result.joint_kmer_freq > 0.0,
            "joint_kmer_freq should be > 0.0, got {}",
            bcl2_result.joint_kmer_freq
        );

        Ok(())
    }

    fn sketch_pair(ksize: u32, scaled: u32) -> (ProteinSketch, ProteinSketch) {
        let (qn, qs) = read_first_fasta_record(TEST_CED9_FASTA).unwrap();
        let (tn, ts) = read_first_fasta_record(TEST_BLC2_FASTA).unwrap();
        (
            ProteinSketch::from_protein_sequence(&qn, &qs, ksize, scaled, "hp_lehninger2").unwrap(),
            ProteinSketch::from_protein_sequence(&tn, &ts, ksize, scaled, "hp_lehninger2").unwrap(),
        )
    }

    fn span(r: &MatchedRegion) -> (u32, u32, u32, u32) {
        (r.start, r.end, r.target_start, r.target_end)
    }

    /// Recounts `n_shared` without the region code: shared k-mer starts on the region's own
    /// diagonal whose whole window lies inside the span.
    fn recount_shared_on_diagonal(
        query: &ProteinSketch,
        target: &ProteinSketch,
        region: &MatchedRegion,
    ) -> u32 {
        let ksize = query.protein_ksize();
        let diagonal = region.target_start as i64 - region.start as i64;
        let mut n = 0;
        for hash in query.intersect(target) {
            for &qpos in &query.kmer_positions()[&hash] {
                let qpos = qpos as u32;
                let inside = region.start <= qpos && qpos + ksize <= region.end;
                let on_diagonal = target.kmer_positions()[&hash]
                    .iter()
                    .any(|&tpos| tpos as i64 - qpos as i64 == diagonal);
                if inside && on_diagonal {
                    n += 1;
                }
            }
        }
        n
    }

    /// A sampled sketch keeps a k-mer by hash value, so at scaled=s roughly 1/s of the shared
    /// k-mers survive and they are rarely adjacent. The sampled path grows each survivor
    /// back out to the full exact match it sits in, so every region it reports must be one of
    /// the dense path's regions, span for span, and never a fragment of one. What it cannot
    /// do is report a match none of whose k-mers survived, so the count only falls.
    ///
    /// CED9 vs BCL2 at hp k=12 has 13 dense regions. The survivors at each scaled are fixed
    /// by the hash cutoff, so they are asserted exactly.
    #[rstest]
    #[case::scaled_2(2, 7)]
    #[case::scaled_5(5, 4)]
    #[case::scaled_10(10, 1)]
    fn test_sampled_regions_are_whole_dense_regions(
        #[case] scaled: u32,
        #[case] expected_regions: usize,
    ) {
        let (dq, dt) = sketch_pair(12, 1);
        let dense = find_matched_regions(&dq, &dt, &dq.intersect(&dt));
        assert_eq!(dense.len(), 13);

        let (sq, st) = sketch_pair(12, scaled);
        let sampled = find_matched_regions(&sq, &st, &sq.intersect(&st));
        assert_eq!(sampled.len(), expected_regions);

        let dense_spans: Vec<_> = dense.iter().map(span).collect();
        for r in &sampled {
            let i = dense_spans
                .iter()
                .position(|&d| d == span(r))
                .unwrap_or_else(|| panic!("sampled region {:?} is not a dense region", span(r)));
            assert_eq!(r.subseq, dense[i].subseq);
            assert_eq!(r.target_subseq, dense[i].target_subseq);
            assert_eq!(r.moltype_seq, dense[i].moltype_seq);
            assert_eq!(r.length, dense[i].length);
            assert_eq!(r.n_shared, recount_shared_on_diagonal(&sq, &st, r), "{:?}", span(r));
            assert!(r.n_shared <= dense[i].n_shared);
        }
    }

    /// The four survivors at scaled=5, with the sampled k-mer count each one rests on.
    /// GVVVCGRMMFSLK kept 2 of its 2 k-mers; the other three kept 1 each, and the 19-residue
    /// QCPMSYGRLIGLISFGGFV match was recovered in full from that single k-mer.
    #[test]
    fn test_sampled_regions_scaled_5_exact() {
        let (q, t) = sketch_pair(12, 5);
        let regions = find_matched_regions(&q, &t, &q.intersect(&t));
        let got: Vec<_> = regions
            .iter()
            .map(|r| (r.start, r.end, r.target_start, r.target_end, r.n_shared, r.subseq.as_str()))
            .collect();
        assert_eq!(
            got,
            [
                (145, 159, 130, 144, 1, "FSLYQDVVRTVGNA"),
                (162, 181, 138, 157, 1, "QCPMSYGRLIGLISFGGFV"),
                (253, 266, 80, 93, 1, "MIGAGVTAGAIGI"),
                (267, 280, 200, 213, 2, "GVVVCGRMMFSLK"),
            ]
        );
    }

    /// At k=15 the only dense region is the 19-residue landmark, held up by 5 shared k-mers.
    /// The sampled path reports it at every scaled here, with the span intact and `n_shared`
    /// falling to however many of the 5 the cutoff kept: all 5 at scaled=2 (the hash values
    /// happen to land low), 3 at scaled=5, 1 at scaled=10.
    #[rstest]
    #[case::scaled_1(1, 5)]
    #[case::scaled_2(2, 5)]
    #[case::scaled_5(5, 3)]
    #[case::scaled_10(10, 1)]
    fn test_sampled_landmark_k15(#[case] scaled: u32, #[case] n_shared: u32) {
        let (q, t) = sketch_pair(15, scaled);
        let regions = find_matched_regions(&q, &t, &q.intersect(&t));
        assert_eq!(regions.len(), 1);
        let r = &regions[0];
        assert_eq!(span(r), (162, 181, 138, 157));
        assert_eq!(r.subseq, "QCPMSYGRLIGLISFGGFV");
        assert_eq!(r.length, 19);
        assert_eq!(r.n_shared, n_shared);
        assert_eq!(r.n_shared, recount_shared_on_diagonal(&q, &t, r));
    }

    /// The seed window's own residues must agree; a run is never grown from a k-mer pair that
    /// only shares a hash. Two identical 20-mers on a diagonal with a substitution between
    /// them are two runs, not one bridged run: exact_run_around stops at the mismatch.
    #[test]
    fn test_exact_run_stops_at_mismatch() {
        // BCL2 positions 138..157 and the same stretch with one residue changed in the middle.
        let q = b"RDGVNWGRIVAFFEFGGVM";
        let t = b"RDGVNWGRIVKFFEFGGVM"; // A -> K at index 10
        let agree = residues_agree("protein20");
        assert_eq!(exact_run_around(q, t, 0, 0, 5, &agree), Some((0, 10)));
        assert_eq!(exact_run_around(q, t, 12, 12, 5, &agree), Some((11, 19)));
        assert_eq!(
            exact_run_around(q, t, 8, 8, 5, &agree),
            None,
            "window 8..13 crosses the mismatch"
        );
        // A seed near the end grows left to the mismatch and right to the sequence end.
        assert_eq!(exact_run_around(q, t, 14, 14, 5, &agree), Some((11, 19)));
    }

    /// An ambiguous residue agrees with either residue it stands for, so a run grows through
    /// it. Under protein20 the raw sequences are compared: BCL2 residues 1-19 with Asp10
    /// written as B run against the real fragment end to end, and a query Z against a target
    /// Asp is a mismatch, since Z stands for Glu or Gln. Under an HP alphabet the stored
    /// sequence already holds the class, so the same run is exact.
    #[test]
    fn test_exact_run_grows_through_an_ambiguous_residue() {
        let q = b"MAHAGRTGYBNREIVMKYI";
        let t = b"MAHAGRTGYDNREIVMKYI";
        let agree = residues_agree("protein20");
        assert_eq!(exact_run_around(q, t, 0, 0, 5, &agree), Some((0, 19)));
        assert_eq!(exact_run_around(q, t, 7, 7, 5, &agree), Some((0, 19)));

        let z = b"MAHAGRTGYZNREIVMKYI";
        assert_eq!(exact_run_around(z, t, 0, 0, 5, &agree), Some((0, 9)));
        assert_eq!(exact_run_around(z, t, 7, 7, 5, &agree), None);

        let hp = residues_agree("hp_lehninger2");
        let q_hp = b"hhphhpphpppphhhhpph";
        assert_eq!(exact_run_around(q_hp, q_hp, 0, 0, 5, &hp), Some((0, 19)));
        // sdm12 keeps Asp and Asn apart, so the stored query keeps its B, and B agrees with
        // the class of D and of N but not with the class of A.
        let sdm12 = residues_agree("sdm12");
        let encode = residue_encoder("sdm12");
        assert!(sdm12(b'B', encode(b'D')));
        assert!(sdm12(b'B', encode(b'N')));
        assert!(!sdm12(b'B', encode(b'A')));
    }

    /// `protein20` sketches store no encoded copy of the sequence (the full alphabet encodes to
    /// itself), so the walk has to fall back to the raw sequence. BCL2 and BCL-xL share the
    /// 16-residue BH1 stretch ELFRDGVNWGRIVAFF, 7 k-mers at k=10. Every sampled sketch that
    /// keeps any of the 7 reports the whole stretch, with `moltype_seq` empty as on the dense
    /// path; at scaled=10 none of the 7 survive and the match is missed outright.
    #[rstest]
    #[case::scaled_1(1, 7)]
    #[case::scaled_2(2, 4)]
    #[case::scaled_5(5, 3)]
    #[case::scaled_10(10, 0)]
    fn test_sampled_regions_protein20_fall_back_to_raw_sequence(
        #[case] scaled: u32,
        #[case] n_shared: u32,
    ) {
        let (qn, qs) = read_fasta_record(TEST_BLC2_FASTA, "BCL2_HUMAN").unwrap();
        let (tn, ts) = read_fasta_record(TEST_FASTA_GZ, "B2CL1_HUMAN").unwrap();
        let q = ProteinSketch::from_protein_sequence(&qn, &qs, 10, scaled, "protein20").unwrap();
        let t = ProteinSketch::from_protein_sequence(&tn, &ts, 10, scaled, "protein20").unwrap();
        assert!(q.get_moltype_sequence().is_none());
        let intersection = q.intersect(&t);
        assert_eq!(intersection.len() as u32, n_shared);

        let regions = find_matched_regions(&q, &t, &intersection);
        if n_shared == 0 {
            assert!(regions.is_empty());
            return;
        }
        assert_eq!(regions.len(), 1);
        let r = &regions[0];
        assert_eq!(span(r), (135, 151, 128, 144));
        assert_eq!(r.subseq, "ELFRDGVNWGRIVAFF");
        assert_eq!(r.target_subseq, "ELFRDGVNWGRIVAFF");
        assert_eq!(r.moltype_seq, "");
        assert_eq!(r.length, 16);
        assert_eq!(r.n_shared, n_shared);
        assert_eq!(r.n_shared, recount_shared_on_diagonal(&q, &t, r));
    }

    /// The sampled path needs the stored sequences to grow a seed; without them it reports
    /// nothing rather than a scatter of single k-mers.
    #[test]
    fn test_sampled_regions_need_stored_sequences() {
        let (q, t) = sketch_pair(12, 5);
        let mut bare = ProteinSketch::new("bare", 12, 5, "hp_lehninger2").unwrap();
        for (h, positions) in q.kmer_positions() {
            bare.kmer_positions_mut().insert(*h, positions.clone());
        }
        assert!(bare.get_moltype_sequence().is_none());
        assert_eq!(find_matched_regions(&bare, &t, &q.intersect(&t)).len(), 0);
    }

    /// Checks `MatchedRegion::n_shared` (surfaced as the CSV's `region_n_shared_kmers`)
    /// against an independent count of real k-mer positions, for every named HP alphabet.
    ///
    /// At scaled=1 FracMinHash keeps every k-mer, so the dense region path sets `n_shared` to
    /// its consecutive-k-mer count, which must equal `length - ksize + 1`. Rather than
    /// trusting that reasoning, this test recomputes the count a different way: for each
    /// region, it walks `kmer_positions` (the sketch's own record of where each retained k-mer
    /// starts) and counts how many positions fall inside the region's span, then asserts
    /// that matches both the field and the formula. Repeated for every named HP
    /// alphabet (including hp_thomas_dill_no_c) since each partitions residues into H/P
    /// differently, and the position bookkeeping has to hold for all of them, not just one.
    #[test]
    fn test_region_shared_kmer_count_exact_at_scaled_one_all_alphabets() {
        use crate::alphabets::Alphabet;

        // Real sequence (already used elsewhere in this codebase for HP-alphabet regression
        // tests, see test_hp_encoding.rs) - self-hit so every alphabet reliably finds a region.
        let seq = "MKTAYIAKQRFLVSNSQLAGKRILVTQADTFMGPTLCEVFAEMG";
        let ksize = 8;

        for alpha in Alphabet::hp_family() {
            let moltype = alpha.to_moltype();
            let sketch = ProteinSketch::from_protein_sequence("self", seq, ksize, 1, moltype)
                .unwrap_or_else(|e| panic!("{moltype}: from_protein_sequence failed: {e}"));

            let shared = sketch.mins_as_set();
            let regions = find_matched_regions(&sketch, &sketch, &shared);
            assert!(!regions.is_empty(), "{moltype}: self-hit should find at least one region");

            let kmer_positions = sketch.kmer_positions();
            for region in &regions {
                let window_end = (region.end as usize).saturating_sub(ksize as usize) + 1;
                let n_in_region = kmer_positions
                    .values()
                    .flatten()
                    .filter(|&&p| p >= region.start as usize && p < window_end)
                    .count();
                assert_eq!(
                    n_in_region as u32, region.n_shared,
                    "{moltype}: region {:?} n_shared mismatch",
                    region.subseq
                );
                assert_eq!(region.n_shared, region.length - ksize + 1);
            }
        }
    }

    /// The same 5 shared k-mers that read as a weak whole-protein match (containment ~0.019
    /// against CED9's 266 k-mers) are the entire signal inside their own 19aa region. Runs a
    /// real database search (needed for region.poisson_score's DB context) and checks the
    /// region-scoped Poisson score independently: recomputes lambda by hand from the searcher's
    /// own background frequencies, restricted to the region's span, and checks it against
    /// region.expected_shared_kmers/region.poisson_score rather than trusting the same code
    /// path that produced them.
    #[test]
    fn test_region_poisson_score_independently_recomputed() -> Result<()> {
        let ksize = 15;
        let scaled = 1;
        let moltype = "hp_lehninger2";

        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(target_index)?;

        let query_index_path = temp_dir.path().join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, ksize, scaled, moltype, true)?;
        query_index.process_fasta(TEST_CED9_FASTA, 0, DEFAULT_BATCH_SIZE)?;
        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        let results = searcher.search(&query_signatures, &SearchFilters::default())?;
        let bcl2_result = results
            .iter()
            .find(|r| r.target_name.contains("BCL2_HUMAN"))
            .expect("Should find BCL2_HUMAN in results");

        // Matches the sourmash-verified BCL2_CED9_K15 constants (see BCL2_CED9_K15 above):
        // 5 shared k-mers, whole-protein containment ~0.019 - a real, weak-looking match.
        assert_eq!(bcl2_result.n_intersecting_hashes, 5);
        assert_relative_eq!(bcl2_result.containment, 0.018796992481203006, epsilon = 1e-9);

        let region = bcl2_result
            .matched_regions
            .iter()
            .find(|r| r.subseq == "QCPMSYGRLIGLISFGGFV")
            .expect("Should find the landmark region");
        assert_eq!(region.length - ksize + 1, 5, "region should account for all 5 shared k-mers");

        // Recompute expected_shared_kmers by hand: sum db_frequency[h]/total_signatures over
        // every retained query k-mer position inside the region's own span - independent of
        // compare()'s implementation, using only the searcher's public stats().
        let stats = searcher.stats();
        let query_sketch = &query_signatures[0];
        let window_end = (region.end as usize).saturating_sub(ksize as usize) + 1;
        let expected_by_hand: f64 = query_sketch
            .kmer_positions()
            .iter()
            .map(|(hashval, positions)| {
                let freq = stats.kmer_frequencies.get(hashval).copied().unwrap_or(1) as f64
                    / stats.total_signatures as f64;
                let n_in_region = positions
                    .iter()
                    .filter(|&&p| p >= region.start as usize && p < window_end)
                    .count();
                freq * n_in_region as f64
            })
            .sum();
        assert_relative_eq!(region.expected_shared_kmers, expected_by_hand, epsilon = 1e-12);

        let pvalue_by_hand = if expected_by_hand > 0.0 {
            Poisson::new(expected_by_hand)
                .map(|dist| (1.0 - dist.cdf(4)).max(0.0)) // k - 1 = 5 - 1 = 4
                .unwrap_or(1.0)
        } else {
            1.0
        };
        let score_by_hand = -pvalue_by_hand.max(f64::MIN_POSITIVE).log10();
        assert_relative_eq!(region.poisson_score, score_by_hand, epsilon = 1e-12);
        assert_relative_eq!(region.tail_probability, pvalue_by_hand, epsilon = 1e-12);

        // Region TF-IDF by hand over the same window: sum ln(N / freq[h]) per retained
        // k-mer position inside the region, and its per-k-mer mean over the 5 shared k-mers.
        let tfidf_by_hand: f64 = query_sketch
            .kmer_positions()
            .iter()
            .map(|(hashval, positions)| {
                let idf = (stats.total_signatures as f64
                    / stats.kmer_frequencies.get(hashval).copied().unwrap_or(1) as f64)
                    .ln();
                let n_in_region = positions
                    .iter()
                    .filter(|&&p| p >= region.start as usize && p < window_end)
                    .count();
                idf * n_in_region as f64
            })
            .sum();
        assert_relative_eq!(region.tfidf, tfidf_by_hand, epsilon = 1e-12);
        assert_relative_eq!(region.mean_idf, tfidf_by_hand / 5.0, epsilon = 1e-12);
        // Pinned values for this fixture (N = 25 signatures): 5 shared k-mers whose database
        // frequencies are a mix of 1 and more than 1, so the mean IDF sits below ln(25) = 3.22.
        assert_relative_eq!(region.tfidf, 10.326058128547245, epsilon = 1e-12);
        assert_relative_eq!(region.mean_idf, 2.0652116257094493, epsilon = 1e-12);

        // The region-scoped null (over ~5 background-frequency k-mers) is a much smaller number
        // than the whole-protein null (over all 266 of CED9's k-mers), so the two numbers are
        // computed from different lambdas and shouldn't coincide.
        assert_ne!(region.poisson_score, bcl2_result.query_poisson_pvalue);

        Ok(())
    }

    /// region_search_space, db_n_targets, db_n_kmers, and run_n_queries are reported as
    /// separate columns rather than multiplied into the p-value, so a hit's reported
    /// significance never shifts depending on what else was in the same search run. Checks
    /// that each column counts what it claims, and that running the same query alone leaves
    /// its p-values and region_search_space unchanged: only run_n_queries moves.
    #[test]
    fn test_multiplicity_components_are_reported_separately() -> Result<()> {
        let ksize = 12;
        let scaled = 1;
        let moltype = "hp_lehninger2";

        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(target_index)?;

        let (ced9_name, ced9_sequence) = read_first_fasta_record(TEST_CED9_FASTA)?;
        let (bcl2_name, bcl2_sequence) = read_first_fasta_record(TEST_BLC2_FASTA)?;
        let ced9 = ProteinSketch::from_protein_sequence(
            &ced9_name,
            &ced9_sequence,
            ksize,
            scaled,
            moltype,
        )?;
        let bcl2 = ProteinSketch::from_protein_sequence(
            &bcl2_name,
            &bcl2_sequence,
            ksize,
            scaled,
            moltype,
        )?;

        let two_query_run = searcher.search(&[ced9.clone(), bcl2], &SearchFilters::default())?;
        let paired = two_query_run
            .iter()
            .find(|r| r.query_name.contains("CED9") && r.target_name.contains("BCL2_HUMAN"))
            .expect("CED9 query should match BCL2 target");

        // TEST_FASTA_GZ holds 25 signatures, and every k-mer occurrence across them is counted.
        assert_eq!(paired.db_n_targets, 25);
        assert_eq!(paired.db_n_kmers, searcher.stats().kmer_frequencies.values().sum::<usize>());
        assert_eq!(paired.run_n_queries, 2);
        // Candidate region start positions in CED9: len - k + 1.
        assert_eq!(paired.region_search_space, ced9_sequence.len() - ksize as usize + 1);

        // Rerunning the same query on its own must not change what the hit itself means.
        let solo_run = searcher.search(&[ced9], &SearchFilters::default())?;
        let solo = solo_run
            .iter()
            .find(|r| r.query_name.contains("CED9") && r.target_name.contains("BCL2_HUMAN"))
            .expect("CED9 query should match BCL2 target when run alone");

        assert_eq!(solo.run_n_queries, 1, "only the run-level count should move");
        assert_eq!(solo.region_search_space, paired.region_search_space);
        assert_eq!(solo.db_n_targets, paired.db_n_targets);
        assert_relative_eq!(
            solo.query_poisson_pvalue,
            paired.query_poisson_pvalue,
            epsilon = 1e-12
        );
        assert_eq!(solo.matched_regions.len(), paired.matched_regions.len());
        for (solo_region, paired_region) in
            solo.matched_regions.iter().zip(paired.matched_regions.iter())
        {
            assert_relative_eq!(
                solo_region.poisson_score,
                paired_region.poisson_score,
                epsilon = 1e-12
            );
        }

        Ok(())
    }

    /// The query p-value and region score are combined with OR, so the full truth table
    /// matters: a single scope clearing is enough, and only both failing rejects. Region
    /// scores below use 0.9 for "weak" (score 0.9 is p ~ 0.126, unremarkable) and 3.1549 for
    /// "strong" (the real BCL2/CED9 landmark: -log10(0.0007)). Exercised directly here because
    /// through `search()` the degenerate caps are hard to reach.
    #[rstest]
    // query passes, region fails -> kept on the query scope
    #[case(0.001, Some(0.9), 0.05, 1.301, true)]
    // query fails, region passes -> kept on the region scope (the BCL2/CED9 shape)
    #[case(0.99, Some(3.1549), 0.05, 1.301, true)]
    // both pass
    #[case(0.001, Some(3.1549), 0.05, 1.301, true)]
    // both fail -> rejected
    #[case(0.99, Some(0.9), 0.05, 1.301, false)]
    // no regions at all: the region disjunct is vacuously false, query alone decides
    #[case(0.001, None, 0.05, 1.301, true)]
    #[case(0.99, None, 0.05, 1.301, false)]
    // min_region_score: INFINITY can never be cleared (no finite score is > infinity), which
    // is how the deprecated --max-pvalue alias reduces to whole-query filtering
    #[case(0.001, Some(3.1549), 0.05, f64::INFINITY, true)]
    #[case(0.99, Some(3.1549), 0.05, f64::INFINITY, false)]
    // NEG_INFINITY accepts anything, including the score = 0.0 (p = 1.0) of a no-DB-context
    // result; max_query_pvalue: INFINITY is the query-scope mirror of that same case
    #[case(1.0, Some(0.0), f64::INFINITY, f64::NEG_INFINITY, true)]
    fn test_scopes_pass_truth_table(
        #[case] query_pvalue: f64,
        #[case] best_region_score: Option<f64>,
        #[case] max_query_pvalue: f64,
        #[case] min_region_score: f64,
        #[case] expected: bool,
    ) {
        let filters =
            SearchFilters { max_query_pvalue, min_region_score, ..SearchFilters::default() };
        assert_eq!(filters.scopes_pass(query_pvalue, best_region_score), expected);
    }

    /// The Poisson survival function is the shared engine behind both the query p-value and
    /// the region score. The degenerate inputs return 1.0 (no evidence) instead of erroring or
    /// producing NaN, which a real search cannot reach but a caller can.
    #[test]
    fn test_poisson_survival_degenerate_inputs_yield_no_evidence() {
        assert_eq!(poisson_survival(0, 2.0), 1.0, "nothing observed is not surprising");
        assert_eq!(poisson_survival(5, 0.0), 1.0, "no null to be surprised against");
        assert_eq!(poisson_survival(5, -1.0), 1.0, "negative rate is not a distribution");
        assert_eq!(poisson_survival(5, f64::NAN), 1.0, "NaN rate is rejected by Poisson::new");
    }

    /// Observing what is expected is unsurprising; observing far more is not. Values are
    /// checked against the closed form of the survival function rather than restated constants.
    #[test]
    fn test_poisson_survival_matches_closed_form() {
        // P(X >= 1 | lambda) = 1 - e^-lambda
        assert_relative_eq!(poisson_survival(1, 0.5), 1.0 - (-0.5f64).exp(), epsilon = 1e-12);
        // P(X >= 2 | lambda) = 1 - e^-lambda(1 + lambda)
        assert_relative_eq!(poisson_survival(2, 0.5), 1.0 - (-0.5f64).exp() * 1.5, epsilon = 1e-12);
        // Strongly enriched observations get vanishing p-values, and p is monotonically
        // decreasing in the observed count for a fixed null.
        assert!(poisson_survival(20, 0.5) < poisson_survival(10, 0.5));
        assert!(poisson_survival(10, 0.5) < poisson_survival(2, 0.5));
        // Every value stays a probability.
        for observed in [1u32, 3, 10] {
            for lambda in [0.1f64, 1.0, 7.5] {
                let pvalue = poisson_survival(observed, lambda);
                assert!((0.0..=1.0).contains(&pvalue), "p={pvalue} out of range");
            }
        }
    }

    /// Enrichment reports 0.0 rather than +inf when there is no expectation to divide by, so
    /// downstream sorting and serialization never see an infinity.
    #[test]
    fn test_fold_enrichment_guards_zero_expectation() {
        assert_eq!(fold_enrichment(5, 0.0), 0.0);
        assert_eq!(fold_enrichment(0, 0.0), 0.0);
        assert_eq!(fold_enrichment(5, 2.0), 2.5);
        assert_eq!(fold_enrichment(0, 2.0), 0.0);
    }

    /// region_expectation only counts k-mers that fit entirely inside the region, which is what
    /// makes the shared-k-mer count exact at scaled=1. Checked against a hand-computed sum: with
    /// one signature in the database every k-mer has frequency 1/1, so lambda is just the number
    /// of query k-mer positions inside the window.
    #[test]
    fn test_region_expectation_counts_only_fully_contained_kmers() -> Result<()> {
        let ksize = 12;
        let temp_dir = TempDir::new()?;
        let index_path = temp_dir.path().join("index");
        let index = ProteomeIndex::new(&index_path, ksize, 1, "hp_lehninger2", true)?;
        index.process_fasta(TEST_CED9_FASTA, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(index)?;

        let (name, sequence) = read_first_fasta_record(TEST_CED9_FASTA)?;
        let sketch =
            ProteinSketch::from_protein_sequence(&name, &sequence, ksize, 1, "hp_lehninger2")?;
        let prefix = searcher.build_position_prefix(&sketch);

        // A window holding one k-mer: [0, 0 + 1) after the ksize adjustment.
        let single = region_expectation(&prefix, 0, ksize, ksize as usize);
        assert_relative_eq!(single, 1.0, epsilon = 1e-12);

        // Widening the region by one residue admits one more k-mer start.
        let double = region_expectation(&prefix, 0, ksize + 1, ksize as usize);
        assert_relative_eq!(double, 2.0, epsilon = 1e-12);

        // A span shorter than k contains no whole k-mer, so there is nothing to expect.
        let too_short = region_expectation(&prefix, 0, ksize - 1, ksize as usize);
        assert_eq!(too_short, 0.0);

        Ok(())
    }

    /// A k-mer window holding an ambiguous residue is sketched under every reading, so one
    /// query position carries several hashes. The prefix arrays must add those up: then the
    /// whole-query prefix total is the same sum `calculate_expected_shared_kmers` and
    /// `calculate_tfidf` take over every query hash, and a region covering the ambiguous
    /// residue counts both readings instead of whichever one HashMap iteration visited last.
    ///
    /// protein20 keeps Asp and Asn distinct, so the B really does yield two different hashes
    /// per window covering it.
    #[test]
    fn test_prefix_sums_add_every_reading_of_an_ambiguous_residue() -> Result<()> {
        let ksize = 5;
        let temp_dir = TempDir::new()?;
        let index_path = temp_dir.path().join("index");
        let index = ProteomeIndex::new(&index_path, ksize, 1, "protein20", true)?;
        index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(index)?;

        // BCL2_HUMAN (P10415) residues 1-32 with the Asp at index 9 written as B (Asp or Asn).
        let with_b = "MAHAGRTGYBNREIVMKYIHYKLSQRGYEWDA";
        let sketch = ProteinSketch::from_protein_sequence("bcl2_b", with_b, ksize, 1, "protein20")?;
        let prepared = searcher.prepare_query(&sketch);

        // 28 windows, plus one extra hash for each of the 5 windows covering the B.
        assert_eq!(sketch.mins_as_set().len(), 33);
        assert_eq!(prepared.position_prefix.len(), 29);

        // Whole-query totals agree with the per-hash sums.
        let expected_total = searcher.calculate_expected_shared_kmers(&sketch, &sketch);
        let total_prefix = *prepared.position_prefix.last().unwrap();
        let total_idf_prefix = *prepared.idf_prefix.last().unwrap();
        assert_relative_eq!(total_prefix, expected_total, epsilon = 1e-12);
        assert_relative_eq!(total_idf_prefix, prepared.tfidf, epsilon = 1e-12);

        // The window [5, 14) holds k-mer starts 5..=9, all covering the B. The Asp readings
        // are BCL2's own k-mers; the Asn readings are in no target, so they count as
        // frequency 1 for lambda and contribute nothing to IDF.
        let k = ksize as usize;
        let lambda = region_expectation(&prepared.position_prefix, 5, 14, k);
        let tfidf = region_expectation(&prepared.idf_prefix, 5, 14, k);
        // 10 hashes at frequency 1 out of 25 targets; 5 Asp-reading hashes at ln(25 / 1).
        assert_relative_eq!(lambda, 10.0 / 25.0, epsilon = 1e-12);
        assert_relative_eq!(tfidf, 5.0 * 25f64.ln(), epsilon = 1e-12);
        assert_relative_eq!(expected_total, 1.4, epsilon = 1e-12);
        assert_relative_eq!(prepared.tfidf, 88.74222873518976, epsilon = 1e-12);

        Ok(())
    }

    /// A region whose expectation is zero must not produce NaN or infinity downstream. This is
    /// the pairing the guards in poisson_survival/fold_enrichment exist for.
    #[test]
    fn test_zero_expectation_region_stays_finite() {
        let pvalue = poisson_survival(5, 0.0);
        let enrichment = fold_enrichment(5, 0.0);
        assert_eq!(pvalue, 1.0);
        assert_eq!(enrichment, 0.0);
        assert!(pvalue.is_finite() && enrichment.is_finite());
    }

    /// all-vs-all searches the database against itself. It had no test at all, and it is one of
    /// the three entry points that has to supply its own query count, so it is the place a
    /// wrong `total_queries` would go unnoticed.
    #[test]
    fn test_search_all_vs_all_skips_self_and_counts_its_own_queries() -> Result<()> {
        let ksize = 12;
        let temp_dir = TempDir::new()?;
        let index_path = temp_dir.path().join("index");
        let index = ProteomeIndex::new(&index_path, ksize, 1, "hp_lehninger2", true)?;
        index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        assert_eq!(index.signature_count(), 25);

        let searcher = ProteinSearcher::new(index)?;
        let results = searcher.search_all_vs_all(&SearchFilters::default())?;
        assert!(!results.is_empty(), "a family database should match itself across members");

        for result in &results {
            // compare() drops self-matches by md5, so no protein may match itself.
            assert_ne!(
                result.query_md5, result.target_md5,
                "self-match leaked into all-vs-all: {}",
                result.query_name
            );
            // Every query is also a target here, so the run count is the database size.
            assert_eq!(result.run_n_queries, 25);
            assert_eq!(result.db_n_targets, 25);
        }

        // Sorted by containment descending, same contract as search().
        for pair in results.windows(2) {
            assert!(pair[0].containment >= pair[1].containment, "results must stay sorted");
        }

        Ok(())
    }

    /// Complements `test_pvalue_scopes_combine_with_or`'s BCL2/CED9 example, where the region
    /// scope rescues a hit the query scope rejects, with a real case running the other
    /// direction. BCL2 vs RTN3 (reticulon-3) at k=9, in the same 25-sequence fixture database,
    /// is an overwhelming whole-protein match (152 shared k-mers scattered across 217 short
    /// regions) with no single region concentrated enough to pass on its own. Its strongest
    /// region only reaches p=0.0533 (score ~1.27), below the ~1.301 default cap (p=0.05).
    /// This shows the OR only needs one scope to hold, in either direction.
    ///
    /// The pair used to be BCL2A1 vs ASPP2 (115 k-mers, 333 regions, best p=0.0956). Chaining
    /// seeds per diagonal merges the pieces a repeated k-mer used to split, and that pair's
    /// best region became 13 residues at p=0.0491, just over the cap.
    #[test]
    fn test_query_scope_alone_keeps_a_diffuse_match_with_no_standout_region() -> Result<()> {
        let ksize = 9;
        let scaled = 1;
        let moltype = "hp_lehninger2";

        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(target_index)?;

        let query_index_path = temp_dir.path().join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, ksize, scaled, moltype, true)?;
        query_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        query_index.load_state()?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        fn find_hit(results: &[SearchResult]) -> Option<&SearchResult> {
            results
                .iter()
                .find(|r| r.query_name.contains("BCL2_HUMAN") && r.target_name.contains("RTN3"))
        }

        // Unfiltered, to inspect the pair's raw numbers.
        let all_results = searcher.search(&query_signatures, &SearchFilters::default())?;
        let hit = find_hit(&all_results).expect("BCL2 vs RTN3 should be found");

        assert_eq!(hit.n_intersecting_hashes, 152);
        assert_eq!(hit.matched_regions.len(), 217);
        assert_relative_eq!(hit.query_poisson_pvalue, 9.961_523_828e-10, epsilon = 1e-18);
        assert!(hit.query_poisson_pvalue < 0.05, "whole-query scope should clearly pass");

        // Bigger poisson_score is more surprising, so the best region is the highest-scoring
        // one (the lowest underlying p-value).
        let best_region_score = hit
            .matched_regions
            .iter()
            .map(|region| region.poisson_score)
            .fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(best_region_score, 1.272_904_943_544_298, epsilon = 1e-12);
        let default_min_region_score = -0.05_f64.log10();
        assert!(
            best_region_score < default_min_region_score,
            "no single region should clear the default cap"
        );

        // Region scope alone: nothing to rescue it, since no region is significant on its own.
        let region_only = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.0,
                min_region_score: default_min_region_score,
                ..SearchFilters::default()
            },
        )?;
        assert!(
            find_hit(&region_only).is_none(),
            "region scope alone should reject a match with no standout region"
        );

        // Query scope alone: the diffuse whole-protein signal is sufficient by itself.
        // min_region_score: INFINITY disables the region scope entirely.
        let query_only = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.05,
                min_region_score: f64::INFINITY,
                ..SearchFilters::default()
            },
        )?;
        assert!(
            find_hit(&query_only).is_some(),
            "query scope alone should keep this diffuse whole-protein match"
        );

        Ok(())
    }
}

#[cfg(test)]
mod ka_calibration_tests {
    use super::*;
    use crate::tests::test_fixtures::TEST_FASTA_GZ;
    use tempfile::TempDir;

    fn shuffled_settings(mismatch_penalty: f64, n_queries: usize) -> KaCalibrationSettings {
        KaCalibrationSettings {
            mismatch_penalty,
            xdrop: 8.0,
            null: DecoyNull::Shuffled,
            reference: DecoyNull::Shuffled,
            reference_shuffles: 1,
            n_queries,
            seed: 1,
        }
    }

    fn searcher_on_first25() -> Result<(TempDir, ProteinSearcher)> {
        let temp_dir = TempDir::new()?;
        let index_path = temp_dir.path().join("index");
        let index = ProteomeIndex::new(&index_path, 12, 1, "hp_lehninger2", true)?;
        index.process_fasta(TEST_FASTA_GZ, 0, 1000)?;
        Ok((temp_dir, ProteinSearcher::new(index)?))
    }

    /// The fit is reproducible from the seed, is stored under its (penalty, X-drop) and
    /// found again by `resolve_ka`, and `--ka-k` wins over it.
    #[test]
    fn test_calibrate_ka_on_first25_is_stored_and_reused() -> Result<()> {
        let (_dir, mut searcher) = searcher_on_first25()?;
        let report = searcher.calibrate_ka(shuffled_settings(2.0, 25))?;
        let fit =
            report.fitted.clone().expect("25 shuffled BCL2 queries give thousands of regions");
        assert_eq!(searcher.extension, None, "calibration restores the extension setting");
        // 25 shuffled queries against the 25 BCL2-family proteins at hp k=12, penalty 2:
        // 9,288 query residues against 8,340 database k-mers, no homolog excess, the line
        // read off the top 8 bins of x = lambda_region S (half a nat each) with at least 30
        // regions, x 7.0 to 11.0 nats.
        //
        // The same 9,561 regions as before lambda went per-region, on a different x axis:
        // each region is now placed by the composition of its own two spans, so the window
        // and the line moved (slope 0.870 -> 0.819, K 0.0178 -> 0.0065). The residual grew
        // with them, 0.054 -> 0.134 in ln count, which is the size counting noise alone
        // gives at the 30-region floor (1/sqrt(30) = 0.18): the old curve was smooth
        // because one lambda per pair put every region of a pair on one scale.
        assert_eq!((fit.n_queries, fit.n_regions), (25, 9561));
        assert_eq!((fit.query_residues, fit.database_kmers), (9288, 8340));
        assert_eq!((fit.score_lo, fit.score_hi, fit.bend_score), (14, 21, None));
        assert_eq!(fit.x_range(), (7.0, 11.0));
        assert!((fit.slope - 0.819_148_194_939_836).abs() < 1e-12, "{}", fit.slope);
        assert!((fit.k - 0.006_512_282_435_710_007_4).abs() < 1e-12, "{}", fit.k);
        assert!((fit.rms_residual - 0.134_104_578_247_110_86).abs() < 1e-12);
        // The BCL2 family is half hydrophobic in the Lehninger classes, so the database's
        // u is 0.5 and the closed-form lambda is ln of the golden ratio. This one stays a
        // whole-database number: it is reported next to the fit, not used to score.
        assert!((fit.match_probability - 0.500_001_136_008_712_7).abs() < 1e-12);
        assert!((fit.lambda_analytic - 0.481_208_536_966_019_95).abs() < 1e-12);
        assert_eq!(fit.survival[0].1, 9561);

        // Same seed, same fit.
        let again = searcher.calibrate_ka(shuffled_settings(2.0, 25))?;
        assert_eq!(again, report);

        searcher.index().put_ka_calibration(&fit)?;
        let (params, source) = searcher.resolve_ka(None, shuffled_settings(2.0, 0))?;
        assert_eq!(params, KaParams { k: fit.k, lambda_scale: fit.slope });
        assert_eq!(source, KaSource::Index(fit.clone()));

        let (params, source) = searcher.resolve_ka(Some(0.03), shuffled_settings(2.0, 0))?;
        assert_eq!((params, source), (KaParams { k: 0.03, lambda_scale: 1.0 }, KaSource::Flag));

        // No stored fit for penalty 3 and no queries allowed: refused, not guessed.
        let err = searcher.resolve_ka(None, shuffled_settings(3.0, 0)).unwrap_err();
        assert!(err.to_string().contains("no Karlin-Altschul fit for penalty 3"), "{err}");
        Ok(())
    }

    /// A refused fit still reports the histogram it refused: one shuffled query gives a
    /// few hundred regions, nearly all in the peak bin, so fewer than MIN_FIT_POINTS tail
    /// bins reach MIN_BIN_COUNT and `fitted` is None, but `survival` is the full curve
    /// and its first entry counts every region.
    #[test]
    fn test_calibrate_ka_refused_fit_still_carries_survival() -> Result<()> {
        let (_dir, mut searcher) = searcher_on_first25()?;
        let report = searcher.calibrate_ka(shuffled_settings(2.0, 1))?;
        assert_eq!(report.fitted, None, "{} regions fitted", report.n_regions);
        assert_eq!(report.n_queries, 1);
        assert!(report.n_regions > 0);
        assert_eq!(report.survival[0].1 as usize, report.n_regions);
        assert!(report.survival.windows(2).all(|w| w[0].1 >= w[1].1), "survival is monotone");
        assert!(report.reference_survival.is_empty(), "no reference under the shuffled null");
        assert!((report.match_probability - 0.500_001_136_008_712_7).abs() < 1e-12);
        assert!(report.query_residues > 0);
        assert_eq!(report.database_kmers, 8340);
        Ok(())
    }

    /// Shuffling each reference query several times multiplies the chance curve and
    /// leaves the real one alone, which is what lets a starved fit reach the 30 regions a
    /// score bin needs. Under the database null with 25 queries: four shuffles give four
    /// times the reference queries and about four times the chance regions, while the
    /// real curve is identical to the one shuffle run.
    #[test]
    fn test_reference_shuffles_multiply_only_the_chance_curve() -> Result<()> {
        let (_dir, mut searcher) = searcher_on_first25()?;
        let settings = |shuffles| KaCalibrationSettings {
            mismatch_penalty: 2.0,
            xdrop: 8.0,
            null: DecoyNull::Database,
            reference: DecoyNull::ShuffledDipeptide,
            reference_shuffles: shuffles,
            n_queries: 25,
            seed: 1,
        };
        let one = searcher.calibrate_ka(settings(1))?;
        let four = searcher.calibrate_ka(settings(4))?;
        assert_eq!((one.n_queries, four.n_queries), (25, 25), "the real side is untouched");
        assert_eq!(one.n_regions, four.n_regions);
        assert_eq!(one.survival, four.survival);
        assert_eq!((one.n_reference_queries, four.n_reference_queries), (25, 100));
        let ratio = four.n_reference_regions as f64 / one.n_reference_regions as f64;
        assert!((3.0..5.0).contains(&ratio), "{ratio} from {one:?} and {four:?}");
        Ok(())
    }
}
