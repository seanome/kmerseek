//! The E-value of an extended region: Karlin-Altschul K and r_database fitted on the index
//! itself, and the decoys the fit is read against. `docs/evalue.md` walks through it with
//! the figures.
//!
//! E = K m n e^(-lambda S) counts the regions with score >= S expected between an
//! unrelated query of m residues and a database of n residues. lambda is solved per pair
//! from the two class compositions (`search::karlin_altschul_lambda`), so a pair of two
//! hydrophobic sequences, whose agreement is what their compositions do by chance, gets
//! lambda 0 and no significance. K, and whether that per-pair lambda has the right scale,
//! are read off the search itself: a few hundred database sequences are searched against
//! the index and every region's normalised score x = lambda_pair S is binned. Under the
//! model the count at x is K L N (1 - e^-w) e^-x, L the calibration residues, N the
//! database residues and w the bin width, so ln(count) against x is a line of slope -1
//! and intercept ln(K L N (1 - e^-w)) (Altschul & Gish 1996; Pearson 1998). Minus the
//! fitted slope is r_database, the factor every pair's lambda is multiplied by at search
//! time; 1 means the closed form holds. Related pairs lift the counts at high x; the fit
//! stops below that, where the real curve starts to rise above the same queries shuffled.
//! The fit goes to the
//! count at each x, not the count at or above it, because a plateau of relatives far up the
//! axis adds a constant to every survival count below it and would flatten the slope.
//!
//! Fitting x rather than the raw score S matters on a proteome: Swiss-Prot holds pairs of
//! membrane and low-complexity proteins whose raw scores run to 60 and beyond with a slope
//! near 0.1, while ordinary pairs fall at 0.45. One line cannot serve both. In x each pair
//! is already on its own scale, and those pairs sit at x = 0.

use serde::{Deserialize, Serialize};

use crate::errors::IndexResult;
use crate::index::ProteomeIndex;
use crate::search::{
    karlin_altschul_lambda, ExtensionScoring, KaCalibrationSettings, KaParams, KaSource,
    ProteinSearcher,
};

/// Karlin-Altschul K for +1 / -penalty scoring with match probability `a`, counting every
/// high-scoring segment of an ungapped comparison (no seed requirement, no give-up margin).
///
/// When the only positive score is +1 every ascending ladder step of the random walk is
/// exactly 1, and K has the closed form E[X e^(lambda X)] (1 - e^-lambda), with
/// E[X e^(lambda X)] = H / lambda (Karlin & Altschul 1990, PNAS 87:2264; this is the
/// `high == 1` branch of BLAST's BlastKarlinLHtoK). For a = 0.3, penalty 1 it reduces to
/// the textbook (q - p)^2 / q = 0.2286. The lattice argument needs a whole-number penalty;
/// any other penalty returns None. None also when no positive lambda exists.
pub fn karlin_altschul_k_theory(a: f64, penalty: f64) -> Option<f64> {
    if penalty <= 0.0 || penalty.fract() != 0.0 {
        return None;
    }
    let lambda = karlin_altschul_lambda(a, penalty);
    if lambda <= 0.0 {
        return None;
    }
    let b = 1.0 - a;
    let mean_score_tilted = a * lambda.exp() - penalty * b * (-penalty * lambda).exp();
    Some(mean_score_tilted * (1.0 - (-lambda).exp()))
}

/// Which sequences are searched to fit r_database and K.
///
/// What each choice keeps and loses, and what it does to the fit, is measured in
/// `docs/evalue.md` and drawn in `docs/images/ka_fit_grid_nulls_by_database.png` (every
/// null against SCOPe40, a Swiss-Prot sample and a UniRef50 sample) and
/// `docs/images/ka_fit_scope40_four_nulls.png`. In short: on SCOPe40 domains the three
/// scrambled nulls agree (r_database 1.04) and real domains give 0.95; on full-length
/// proteins real sequences give 0.83 to 0.87, a plain shuffle 1.0, and keeping dipeptides
/// already pulls the shuffle to 0.94, so hydrophobic runs alone explain a third of the gap.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
pub enum DecoyNull {
    /// Database sequences as they are, searched against the index. Everything real stays
    /// in; the fit stops where the counts start to rise above the same queries shuffled,
    /// which are searched alongside as the reference (Altschul's tutorial, approach i;
    /// Collins et al. 1988; Pearson 1998). Costs two calibration searches.
    Database,
    /// Each query is a database sequence with its residues shuffled: the independent-letter
    /// model BLAST's tables are fitted on (Altschul & Gish 1996). Too easy a null on real
    /// proteins, whose hydrophobic runs and helix and strand periodicity a shuffle destroys.
    Shuffled,
    /// Each query is a database sequence shuffled so that every dipeptide (each pair of
    /// neighbouring residues) occurs as often as in the original (Altschul & Erickson 1985;
    /// sampled as a random Eulerian path, Kandel et al. 1996, the uShuffle k = 2 method).
    /// Keeps the rate at which a hydrophobic residue follows a hydrophobic one, and with it
    /// the lengths of hydrophobic runs. The default reference for `Database`.
    ShuffledDipeptide,
    /// Each query is a database sequence read back to front. In a hydrophobic/polar
    /// alphabet a helix or a strand reads much the same backwards, so reversed family
    /// members still match the query one element at a time; a check, not a null to fit on.
    Reversed,
}

impl std::fmt::Display for DecoyNull {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DecoyNull::Database => write!(f, "database"),
            DecoyNull::Reversed => write!(f, "reversed"),
            DecoyNull::Shuffled => write!(f, "shuffled"),
            DecoyNull::ShuffledDipeptide => write!(f, "shuffled-dipeptide"),
        }
    }
}

/// One fitted (r_database, K), stored in the index under `ka_calibration` and looked up
/// at search time by its `scoring`. The seed length and alphabet are the index's.
///
/// The fit is a straight line through ln(regions in each bin of x = lambda_pair S).
/// Minus its slope is `r_database`; its height gives `k` (see `line_through`). Bins are
/// `bin_width` nats wide, so every field named "score" below is a bin index: bin b covers
/// x in [b w, (b + 1) w).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KaCalibration {
    /// Mismatch penalty and give-up margin the fit is for. A search with a different pair
    /// cannot use this fit.
    pub scoring: ExtensionScoring,
    /// What the calibration queries were.
    pub null: DecoyNull,
    /// Seed of the query sampler, so the fit can be reproduced.
    pub seed: u64,
    /// Calibration queries searched.
    pub n_queries: usize,
    /// Residues in the calibration queries added up: L in ln(K L N).
    pub query_residues: u64,
    /// The database size the E-value uses for n: the index's k-mer count stands in for its
    /// residue count.
    pub database_kmers: u64,
    /// Regions the calibration queries produced, at any score.
    pub n_regions: usize,
    /// Chance that two positions drawn from the sampled database sequences share a class
    /// (a of the database against itself). For reporting only: the search solves lambda
    /// per pair.
    pub match_probability: f64,
    /// The closed-form lambda at `match_probability`, for reporting next to `r_database`.
    pub lambda_analytic: f64,
    /// Minus the slope of ln(count) against x = lambda_pair S, per nat. 1 means the
    /// closed-form per-pair lambda has the right scale; a search multiplies every pair's
    /// lambda by this.
    pub r_database: f64,
    /// Karlin-Altschul K: the line's height with the query residues, database size and bin
    /// width divided out.
    pub k: f64,
    /// Width of one x bin in nats (`BIN_WIDTH`).
    pub bin_width: f64,
    /// First bin the line was fitted on.
    pub score_lo: i64,
    /// Last bin the line was fitted on, inclusive.
    pub score_hi: i64,
    /// First bin above the fit that sat above the line: where related pairs begin. None
    /// when the fit ran out of counts first.
    pub bend_score: Option<i64>,
    /// Root mean square of the fit residuals in ln count.
    pub rms_residual: f64,
    /// Regions with score >= s, for every bin s from the smallest score seen up to the
    /// largest, for plotting the curve the fit was read from.
    pub survival: Vec<(i64, u64)>,
    /// For the `Database` null: the same queries shuffled, searched the same way. The ratio
    /// of the two curves is what decides where the fit stops
    /// (`fit_scores_with_reference`). Empty for other nulls.
    pub reference_survival: Vec<(i64, u64)>,
    /// Minus the slope of the reference curve over the fit window, per bin. None for nulls
    /// other than `Database`.
    pub reference_lambda: Option<f64>,
    /// How the reference queries were made (`Shuffled` or `ShuffledDipeptide`). None for
    /// nulls other than `Database`.
    pub reference: Option<DecoyNull>,
}

impl KaCalibration {
    /// The two numbers a search takes from the fit.
    pub fn ka_params(&self) -> KaParams {
        KaParams { k: self.k, r_database: self.r_database }
    }

    /// The fit window in nats of x = lambda_pair S.
    pub fn x_range(&self) -> (f64, f64) {
        (self.score_lo as f64 * self.bin_width, (self.score_hi + 1) as f64 * self.bin_width)
    }

    /// Number of score bins the line was fitted on.
    pub fn n_fit_points(&self) -> i64 {
        self.score_hi - self.score_lo + 1
    }

    /// A record with nothing in it but its key and K, for tests of how fits are stored
    /// and looked up. No number in it comes from a fit; a real one is in
    /// `search::tests::test_calibrate_ka_on_first25_is_stored_and_reused`.
    #[cfg(test)]
    pub(crate) fn placeholder(scoring: ExtensionScoring, k: f64) -> Self {
        Self {
            scoring,
            null: DecoyNull::Shuffled,
            seed: 0,
            n_queries: 0,
            query_residues: 0,
            database_kmers: 0,
            n_regions: 0,
            match_probability: 0.0,
            lambda_analytic: 0.0,
            r_database: 0.0,
            k,
            bin_width: BIN_WIDTH,
            score_lo: 0,
            score_hi: 0,
            bend_score: None,
            rms_residual: 0.0,
            survival: Vec::new(),
            reference_survival: Vec::new(),
            reference_lambda: None,
            reference: None,
        }
    }
}

/// Width of one bin of the normalised score x = lambda_pair S, in nats. Half a nat is about
/// one raw score unit at lambda 0.45.
pub const BIN_WIDTH: f64 = 0.5;

/// Fewest regions a score bin needs to enter the fit; below this the Poisson noise in
/// ln count (about 1 / sqrt(count)) is larger than the effects being fitted.
pub const MIN_BIN_COUNT: u64 = 30;

/// The fit uses at most this many score bins, the highest ones below the related pairs,
/// so the slope is read as close to the decision tail as those pairs allow.
pub const FIT_WINDOW: i64 = 8;

/// Fewest bins a fit is accepted on.
pub const MIN_FIT_POINTS: i64 = 4;

/// A bin whose ln count sits more than this many Poisson standard deviations above the
/// line fitted to the bins below it is where related pairs begin.
const BEND_SIGMAS: f64 = 2.0;

/// Slack added to the bend test, in ln count, so a bin one region above the line at large
/// counts does not end the fit.
const BEND_SLACK: f64 = 0.05;

/// The line fitted to ln(regions in the bin at score S) against S, and where it was read.
/// Scores here are bin indices, not nats; `ProteinSearcher::run_calibration` converts.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoreFit {
    /// Minus the fitted slope, per bin. Becomes `KaCalibration::r_database` once divided
    /// by the bin width.
    pub lambda: f64,
    /// The line's value at score 0: ln(K L N (1 - e^-lambda)). Becomes
    /// `KaCalibration::k` once L, N and the bin width are divided out.
    pub ln_intercept: f64,
    /// First bin the line was fitted on.
    pub score_lo: i64,
    /// Last bin the line was fitted on, inclusive.
    pub score_hi: i64,
    /// First bin above the fit that sat above the line; None when the fit ran out of
    /// counts first.
    pub bend_score: Option<i64>,
    /// Root mean square of the fit residuals in ln count.
    pub rms_residual: f64,
    /// Regions with score >= s for every bin s, the curve the fit was read from.
    pub survival: Vec<(i64, u64)>,
    /// The reference curve when the fit stopped against one
    /// (`fit_scores_with_reference`); empty otherwise.
    pub reference_survival: Vec<(i64, u64)>,
    /// Minus the reference curve's slope over the same window, per bin; None without a
    /// reference.
    pub reference_lambda: Option<f64>,
}

/// Count of regions with score >= s for every integer s from the smallest to the largest
/// score in `scores`. Fractional scores (a non-integer penalty) fall into the bin below.
pub fn survival_counts(scores: &[f64]) -> Vec<(i64, u64)> {
    let bins: Vec<i64> = scores.iter().map(|s| s.floor() as i64).collect();
    let (Some(&lo), Some(&hi)) = (bins.iter().min(), bins.iter().max()) else {
        return Vec::new();
    };
    let mut per_bin = vec![0u64; (hi - lo + 1) as usize];
    for b in bins {
        per_bin[(b - lo) as usize] += 1;
    }
    let mut running = 0u64;
    let mut out: Vec<(i64, u64)> = Vec::with_capacity(per_bin.len());
    for (i, c) in per_bin.iter().enumerate().rev() {
        running += c;
        out.push((lo + i as i64, running));
    }
    out.reverse();
    out
}

/// Count of regions with score exactly s (fractional scores in the bin below), for every
/// integer s from the smallest to the largest score in `scores`.
pub fn bin_counts(scores: &[f64]) -> Vec<(i64, u64)> {
    let survival = survival_counts(scores);
    survival
        .iter()
        .enumerate()
        .map(|(i, &(s, at_least))| {
            let above = survival.get(i + 1).map_or(0, |&(_, c)| c);
            (s, at_least - above)
        })
        .collect()
}

/// The bins above the most populated one. Below the peak sit the bare seeds and the pairs
/// whose lambda is small; the Karlin-Altschul tail is what comes after it.
fn tail_bins(bins: &[(i64, u64)]) -> Vec<(i64, u64)> {
    match bins.iter().enumerate().max_by_key(|(_, b)| b.1) {
        Some((peak, _)) => bins[peak + 1..].to_vec(),
        None => Vec::new(),
    }
}

/// Least squares of ln count on score over `points`; returns (slope, intercept).
///
/// Under the model ln(count at score s) = ln(K L N (1 - e^-w)) - r_database s, so the two
/// numbers map onto the two constants: minus the slope is r_database (the factor on every
/// pair's closed-form lambda), and the intercept, once L, N and the bin width w are
/// divided out, is K (Altschul & Gish 1996; Pearson 1998).
fn line_through(points: &[(i64, u64)]) -> (f64, f64) {
    let n = points.len() as f64;
    let (sx, sy) =
        points.iter().fold((0.0, 0.0), |(sx, sy), &(s, c)| (sx + s as f64, sy + (c as f64).ln()));
    let (mx, my) = (sx / n, sy / n);
    let (sxx, sxy) = points.iter().fold((0.0, 0.0), |(sxx, sxy), &(s, c)| {
        let dx = s as f64 - mx;
        (sxx + dx * dx, sxy + dx * ((c as f64).ln() - my))
    });
    let slope = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    (slope, my - slope * mx)
}

fn rms_residual(points: &[(i64, u64)], slope: f64, intercept: f64) -> f64 {
    let ss: f64 = points
        .iter()
        .map(|&(s, c)| {
            let r = (c as f64).ln() - (intercept + slope * s as f64);
            r * r
        })
        .sum();
    (ss / points.len() as f64).sqrt()
}

/// Fit ln(regions with score S) against S, stopping below the related pairs.
///
/// Bins with fewer than `MIN_BIN_COUNT` regions are ignored, and so is everything up to and
/// including the most populated bin: below it sit the bare seeds and the pairs whose
/// lambda is small. Starting from the next `MIN_FIT_POINTS` bins, the line is extended one
/// bin at a time upward and a bin joins
/// while its ln count is within `BEND_SIGMAS / sqrt(count) + BEND_SLACK` above the line
/// fitted so far (below it is fine: the seed requirement makes the true curve concave).
/// The first bin above that is where homologs start to show and where the fit stops. The
/// line is then read off the highest `FIT_WINDOW` bins that made it in, the ones nearest
/// the scores that decide a hit. None when fewer than `MIN_FIT_POINTS` bins qualify.
pub fn fit_scores(scores: &[f64]) -> Option<ScoreFit> {
    let usable: Vec<(i64, u64)> =
        tail_bins(&bin_counts(scores)).into_iter().filter(|&(_, c)| c >= MIN_BIN_COUNT).collect();
    if (usable.len() as i64) < MIN_FIT_POINTS {
        return None;
    }
    let mut accepted = MIN_FIT_POINTS as usize;
    let mut bend_score = None;
    while accepted < usable.len() {
        let (slope, intercept) = line_through(&usable[..accepted]);
        let (s, c) = usable[accepted];
        let excess = (c as f64).ln() - (intercept + slope * s as f64);
        if excess > BEND_SIGMAS / (c as f64).sqrt() + BEND_SLACK {
            bend_score = Some(s);
            break;
        }
        accepted += 1;
    }
    let window = &usable[accepted.saturating_sub(FIT_WINDOW as usize)..accepted];
    let (slope, ln_intercept) = line_through(window);
    if slope >= 0.0 {
        return None;
    }
    Some(ScoreFit {
        lambda: -slope,
        ln_intercept,
        score_lo: window[0].0,
        score_hi: window[window.len() - 1].0,
        bend_score,
        rms_residual: rms_residual(window, slope, ln_intercept),
        survival: survival_counts(scores),
        reference_survival: Vec::new(),
        reference_lambda: None,
    })
}

/// Bins the ratio baseline is fitted on: the lowest usable ones, where real and shuffled
/// sequences agree best (on SCOPe40 their slopes differ by 1% over scores 20 to 27).
const REFERENCE_BASE: usize = 8;

/// Fit ln(regions with score S) against S for `scores`, stopping where the curve starts to
/// rise relative to `reference`, the same queries shuffled.
///
/// The single-curve test in `fit_scores` sees a step. Relatives below 40% identity come
/// in as a ramp: on SCOPe40 the real-sequence slope falls away from the shuffled one from
/// score 32 while the step test fires at 40. Here the ratio ln(count) - ln(reference count)
/// is taken per bin; a line through its lowest `REFERENCE_BASE` usable bins is the
/// baseline (a constant ratio is a K difference, a gentle slope is sequence structure),
/// and a bin joins while its ratio is within `BEND_SIGMAS` Poisson standard deviations
/// plus `BEND_SLACK` of that baseline. Because the reference has the same finite-length
/// concavity as the real curve, the ratio also cancels that, which the step test could not.
/// The line is then read off the highest `FIT_WINDOW` accepted bins, as in `fit_scores`.
pub fn fit_scores_with_reference(scores: &[f64], reference: &[f64]) -> Option<ScoreFit> {
    let reference_bins: std::collections::HashMap<i64, u64> =
        bin_counts(reference).into_iter().collect();
    let usable: Vec<(i64, u64, u64)> = tail_bins(&bin_counts(scores))
        .into_iter()
        .filter_map(|(s, c)| {
            let r = *reference_bins.get(&s)?;
            (c >= MIN_BIN_COUNT && r >= MIN_BIN_COUNT).then_some((s, c, r))
        })
        .collect();
    if (usable.len() as i64) < MIN_FIT_POINTS {
        return None;
    }
    let ratio = |&(s, c, r): &(i64, u64, u64)| (s as f64, (c as f64).ln() - (r as f64).ln());
    let base: Vec<(f64, f64)> = usable.iter().take(REFERENCE_BASE).map(ratio).collect();
    let (slope, intercept) = line_through_f64(&base);
    let mut accepted = base.len();
    let mut bend_score = None;
    while accepted < usable.len() {
        let (s, c, r) = usable[accepted];
        let (x, y) = ratio(&usable[accepted]);
        let sigma = (1.0 / c as f64 + 1.0 / r as f64).sqrt();
        if y - (intercept + slope * x) > BEND_SIGMAS * sigma + BEND_SLACK {
            bend_score = Some(s);
            break;
        }
        accepted += 1;
    }
    let window: Vec<(i64, u64)> = usable[accepted.saturating_sub(FIT_WINDOW as usize)..accepted]
        .iter()
        .map(|&(s, c, _)| (s, c))
        .collect();
    let reference_window: Vec<(i64, u64)> = usable
        [accepted.saturating_sub(FIT_WINDOW as usize)..accepted]
        .iter()
        .map(|&(s, _, r)| (s, r))
        .collect();
    let (slope, ln_intercept) = line_through(&window);
    if slope >= 0.0 {
        return None;
    }
    Some(ScoreFit {
        lambda: -slope,
        ln_intercept,
        score_lo: window[0].0,
        score_hi: window[window.len() - 1].0,
        bend_score,
        rms_residual: rms_residual(&window, slope, ln_intercept),
        survival: survival_counts(scores),
        reference_survival: survival_counts(reference),
        reference_lambda: Some(-line_through(&reference_window).0),
    })
}

/// Least squares of y on x; returns (slope, intercept).
fn line_through_f64(points: &[(f64, f64)]) -> (f64, f64) {
    let n = points.len() as f64;
    let (mx, my) =
        (points.iter().map(|p| p.0).sum::<f64>() / n, points.iter().map(|p| p.1).sum::<f64>() / n);
    let (sxx, sxy) = points.iter().fold((0.0, 0.0), |(sxx, sxy), &(x, y)| {
        (sxx + (x - mx) * (x - mx), sxy + (x - mx) * (y - my))
    });
    let slope = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    (slope, my - slope * mx)
}

/// Deterministic generator for picking calibration queries, so an index built twice from
/// the same FASTA stores the same fit. splitmix64 (Steele, Lea & Flood 2014).
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniform index below `n`.
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// `count` distinct indices below `n` (all of them when `count >= n`), in order.
    pub fn sample_indices(&mut self, n: usize, count: usize) -> Vec<usize> {
        if count >= n {
            return (0..n).collect();
        }
        let mut chosen = std::collections::BTreeSet::new();
        while chosen.len() < count {
            chosen.insert(self.below(n));
        }
        chosen.into_iter().collect()
    }

    /// Fisher-Yates shuffle.
    pub fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.below(i + 1);
            items.swap(i, j);
        }
    }
}

/// The calibration query built from one database sequence under `null`.
pub fn make_decoy(raw: &str, null: DecoyNull, rng: &mut SplitMix64) -> String {
    match null {
        DecoyNull::Database => raw.to_string(),
        DecoyNull::Reversed => raw.chars().rev().collect(),
        DecoyNull::Shuffled => {
            let mut residues: Vec<char> = raw.chars().collect();
            rng.shuffle(&mut residues);
            residues.into_iter().collect()
        }
        DecoyNull::ShuffledDipeptide => dipeptide_shuffle(raw, rng),
    }
}

/// A uniformly random rearrangement of `raw` with the same dipeptide counts, the same first
/// residue and the same last residue (Kandel, Matias, Unger & Winkler 1996).
///
/// The residues are the vertices of a graph and each neighbouring pair an edge, so a
/// rearrangement with the same dipeptide counts is a path that uses every edge once. Such a
/// path exists when the edges that leave each vertex for the last time form a tree pointing
/// at the final residue. That tree is drawn with Wilson's loop-erased random walk, weighted
/// by edge multiplicity, which is the distribution the theorem needs; the other edges out
/// of each vertex are then shuffled and the path is walked from the first residue.
pub fn dipeptide_shuffle(raw: &str, rng: &mut SplitMix64) -> String {
    let bytes = raw.as_bytes();
    if bytes.len() < 3 {
        return raw.to_string();
    }
    let last = bytes[bytes.len() - 1];
    // Outgoing edges per residue, as the residue they lead to, in sequence order.
    let mut out: [Vec<u8>; 256] = std::array::from_fn(|_| Vec::new());
    for pair in bytes.windows(2) {
        out[pair[0] as usize].push(pair[1]);
    }
    let next = last_edges(&out, last, rng);
    for (v, edges) in out.iter_mut().enumerate() {
        if edges.is_empty() {
            continue;
        }
        if v as u8 != last {
            let target = next[v];
            let pos = edges.iter().position(|&t| t == target).expect("last edge is an edge");
            edges.swap_remove(pos);
            rng.shuffle(edges);
            edges.push(target);
        } else {
            rng.shuffle(edges);
        }
    }
    let mut cursor = [0usize; 256];
    let mut path = Vec::with_capacity(bytes.len());
    let mut v = bytes[0];
    path.push(v);
    for _ in 1..bytes.len() {
        let t = out[v as usize][cursor[v as usize]];
        cursor[v as usize] += 1;
        path.push(t);
        v = t;
    }
    String::from_utf8(path).expect("a rearrangement of ASCII residues")
}

/// Wilson's algorithm: for every residue with outgoing edges (other than `root`), the edge
/// it leaves by for the last time, drawn as a random spanning tree pointing at `root`.
fn last_edges(out: &[Vec<u8>; 256], root: u8, rng: &mut SplitMix64) -> [u8; 256] {
    let mut next = [0u8; 256];
    let mut in_tree = [false; 256];
    in_tree[root as usize] = true;
    for start in 0..256usize {
        if out[start].is_empty() || in_tree[start] {
            continue;
        }
        let mut u = start;
        while !in_tree[u] {
            let edges = &out[u];
            next[u] = edges[rng.below(edges.len())];
            u = next[u] as usize;
        }
        let mut u = start;
        while !in_tree[u] {
            in_tree[u] = true;
            u = next[u] as usize;
        }
    }
    next
}

/// Say so when the related pairs left the fit fewer bins than `FIT_WINDOW`: the slope is
/// then read from the seed end of the curve, where the seed requirement still shapes it.
pub fn warn_on_short_fit(fit: &KaCalibration) {
    if fit.n_fit_points() < FIT_WINDOW {
        eprintln!(
            "  WARNING: the fit has only {} bins (x {:.1}..{:.1}) below the relatives at x {}. \
             Related sequences are dense in this database; the slope is read close to the \
             seed. More --ka-queries gives a second opinion.",
            fit.n_fit_points(),
            fit.x_range().0,
            fit.x_range().1,
            fit.bend_score
                .map_or("none".to_string(), |b| format!("{:.1}", b as f64 * fit.bin_width))
        );
    }
}

/// Fit r_database and K on `settings.n_queries` calibration queries of the index just
/// built and store the fit in the index. Also prints the closed-form lambda and K for the
/// database's own composition, so the effect of the seed requirement and of real sequence
/// structure on each is visible. `survival_out` gets the curve the fit was read from.
pub fn calibrate_index(
    index: ProteomeIndex,
    settings: KaCalibrationSettings,
    survival_out: Option<&std::path::Path>,
) -> IndexResult<()> {
    let KaCalibrationSettings { scoring, null, reference, n_queries, .. } = settings;
    let ExtensionScoring { mismatch_penalty, xdrop } = scoring;
    eprintln!(
        "Fitting r_database and K on {n_queries} {null} sequences (mismatch penalty {mismatch_penalty}, give-up margin {xdrop}){}...",
        if null == DecoyNull::Database {
            format!("; the fit stops where their counts rise above the same sequences {reference}")
        } else {
            String::new()
        }
    );
    let mut searcher = ProteinSearcher::new(index)?;
    let report = searcher.calibrate_ka(settings)?;
    let theory_k = karlin_altschul_k_theory(report.match_probability, mismatch_penalty)
        .map_or("none".to_string(), |k| format!("{k:.4}"));
    eprintln!(
        "  Closed form at the database's own match probability {:.3}: K {theory_k} (independent positions, no seed, one lambda for every pair)",
        report.match_probability
    );
    match report.fitted {
        Some(fit) => {
            eprintln!("  {}", KaSource::Fitted(fit.clone()));
            warn_on_short_fit(&fit);
            if let Some(path) = survival_out {
                write_survival_csv(path, &fit)?;
                eprintln!("  Survival curve written to {}", path.display());
            }
            searcher.index().put_ka_calibration(&fit)?;
            eprintln!(
                "  Stored in the index for --extend-mismatch-penalty {mismatch_penalty} --extend-xdrop {xdrop}"
            );
        }
        None => eprintln!(
            "  {} queries gave only {} regions, too few score bins to fit; nothing stored. \
             A search will have to fit its own r_database and K (--ka-queries) or be given --ka-k.",
            report.n_queries, report.n_regions
        ),
    }
    Ok(())
}

/// One row per bin of x = lambda_pair S: the count of regions at or above it, the fitted
/// line's count, the reference count, and whether the bin was inside the fit window.
/// `scripts/plot_ka_survival.py` draws it.
pub fn write_survival_csv(path: &std::path::Path, fit: &KaCalibration) -> IndexResult<()> {
    let mut w = csv::Writer::from_path(path)?;
    w.write_record([
        "x",
        "n_regions_at_least",
        "fitted_n_regions_at_least",
        "reference_n_regions_at_least",
        "in_fit",
        "r_database",
        "k",
        "lambda_analytic",
        "match_probability",
        "null",
        "reference",
        "mismatch_penalty",
        "xdrop",
        "n_queries",
        "query_residues",
        "database_kmers",
        "bin_width",
    ])?;
    // The fit is a line through ln(regions in the bin at x); its survival is the same line
    // divided by (1 - e^(-r_database w)).
    let per_bin = 1.0 - (-fit.r_database * fit.bin_width).exp();
    let ln_intercept =
        (fit.k * fit.query_residues as f64 * fit.database_kmers as f64 * per_bin).ln();
    let reference: std::collections::HashMap<i64, u64> =
        fit.reference_survival.iter().copied().collect();
    for &(bin, count) in &fit.survival {
        let x = bin as f64 * fit.bin_width;
        let fitted = (ln_intercept - fit.r_database * x).exp() / per_bin;
        w.write_record([
            format!("{x:.3}"),
            count.to_string(),
            format!("{fitted:.3}"),
            reference.get(&bin).map_or(String::new(), |r| r.to_string()),
            (fit.score_lo <= bin && bin <= fit.score_hi).to_string(),
            fit.r_database.to_string(),
            fit.k.to_string(),
            fit.lambda_analytic.to_string(),
            fit.match_probability.to_string(),
            fit.null.to_string(),
            fit.reference.map_or(String::new(), |r| r.to_string()),
            fit.scoring.mismatch_penalty.to_string(),
            fit.scoring.xdrop.to_string(),
            fit.n_queries.to_string(),
            fit.query_residues.to_string(),
            fit.database_kmers.to_string(),
            fit.bin_width.to_string(),
        ])?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::test_fixtures::TEST_BLC2_FASTA;

    /// Steps of +1 with probability p and -1 with probability q = 1 - p, p < q. Karlin &
    /// Altschul 1990 (PNAS 87:2264) give lambda as the root of p e^lambda + q e^-lambda = 1,
    /// so e^lambda = q / p, and K = E[X e^(lambda X)] (1 - e^-lambda) when the only positive
    /// step is +1. Written out: E[X e^(lambda X)] = p (q/p) - q (p/q) = q - p and
    /// 1 - e^-lambda = 1 - p/q, so K = (q - p)(q - p) / q. Ewens & Grant, Statistical
    /// Methods in Bioinformatics (2nd ed., 2005), reach the same (q - p)^2 / q for this walk
    /// in their BLAST chapter.
    #[test]
    fn test_k_theory_matches_the_plus_minus_one_random_walk() {
        let (p, q) = (0.3, 0.7);
        let k = karlin_altschul_k_theory(p, 1.0).unwrap();
        assert!((k - (q - p) * (q - p) / q).abs() < 1e-12, "{k}");
    }

    /// Balanced two-class composition (a = 0.5) with penalty 2: the lambda equation
    /// 0.5 e^lambda + 0.5 e^-2lambda = 1 becomes y^3 - 2 y^2 + 1 = 0 in y = e^lambda, whose
    /// root above 1 is the golden ratio phi. Then E[X e^(lambda X)] = a phi - 2 (1 - a) / phi^2
    /// = phi / 2 - 1 / phi^2 and 1 - e^-lambda = 1 - 1 / phi; K is their product, 0.1631.
    /// A simulated 4-million-step walk with these steps gave 458 high-scoring segments where
    /// this K predicts 478.
    #[test]
    fn test_k_theory_balanced_hp_at_penalty_two() {
        let (a, penalty) = (0.5, 2.0);
        let phi = (1.0 + 5f64.sqrt()) / 2.0;
        let expected = (a * phi - penalty * (1.0 - a) / (phi * phi)) * (1.0 - 1.0 / phi);
        let k = karlin_altschul_k_theory(a, penalty).unwrap();
        assert!((k - expected).abs() < 1e-12, "{k} vs {expected}");
        assert!((k - 0.1631).abs() < 5e-4, "{k}");
    }

    #[test]
    fn test_k_theory_none_without_positive_lambda_or_whole_penalty() {
        assert_eq!(karlin_altschul_k_theory(2.0 / 3.0, 2.0), None);
        assert_eq!(karlin_altschul_k_theory(0.9, 2.0), None);
        assert_eq!(karlin_altschul_k_theory(0.5, 1.5), None);
        assert_eq!(karlin_altschul_k_theory(0.5, 0.0), None);
    }

    #[test]
    fn test_survival_counts() {
        assert_eq!(
            survival_counts(&[12.0, 12.0, 13.0, 15.5, 15.0]),
            vec![(12, 5), (13, 3), (14, 2), (15, 2)]
        );
        assert!(survival_counts(&[]).is_empty());
    }

    /// Made-up counts that follow the model exactly: N_SEED regions at score >= SEED_SCORE,
    /// falling as e^(-TRUE_LAMBDA (s - SEED_SCORE)) up to score 30, so each bin holds
    /// N_SEED (1 - e^-TRUE_LAMBDA) e^(-TRUE_LAMBDA (s - SEED_SCORE)) regions. Then a bump
    /// of BUMP extra regions at score 27, standing in for related pairs. The fit reads
    /// TRUE_LAMBDA off the bins below the bump, to within rounding of the counts to whole
    /// regions, and the bump does not reach the bins below it the way it would on the
    /// survival curve.
    #[test]
    fn test_fit_scores_recovers_slope_and_stops_at_homolog_excess() {
        const N_SEED: f64 = 4000.0;
        const SEED_SCORE: i64 = 12;
        const TRUE_LAMBDA: f64 = 0.4;
        const BUMP: usize = 400;
        let at_least =
            |x: i64| (N_SEED * (-TRUE_LAMBDA * (x - SEED_SCORE) as f64).exp()).round() as i64;
        let mut scores = Vec::new();
        for s in SEED_SCORE..=30 {
            let exactly = (at_least(s) - at_least(s + 1)).max(0) as usize;
            scores.extend(std::iter::repeat_n(s as f64, exactly));
        }
        let clean = fit_scores(&scores).unwrap();
        // Bins 13..=21 hold at least MIN_BIN_COUNT regions; the top FIT_WINDOW of them are
        // the window.
        assert_eq!((clean.score_lo, clean.score_hi, clean.bend_score), (14, 21, None));
        assert!((clean.lambda - TRUE_LAMBDA).abs() < 0.01, "{}", clean.lambda);
        let expected_intercept =
            (N_SEED * (1.0 - (-TRUE_LAMBDA).exp())).ln() + TRUE_LAMBDA * SEED_SCORE as f64;
        assert!((clean.ln_intercept - expected_intercept).abs() < 0.1, "{}", clean.ln_intercept);
        assert!(clean.rms_residual < 0.03, "{}", clean.rms_residual);
        assert_eq!(clean.survival[0], (SEED_SCORE, scores.len() as u64));

        scores.extend(std::iter::repeat_n(27.0, BUMP));
        let bent = fit_scores(&scores).unwrap();
        assert_eq!((bent.score_lo, bent.score_hi, bent.bend_score), (14, 21, Some(27)));
        assert!((bent.lambda - clean.lambda).abs() < 1e-12);
        assert_eq!(bent.survival[0], (SEED_SCORE, scores.len() as u64));
    }

    /// Made-up counts again. Reference: N_SEED e^(-TRUE_LAMBDA (s - SEED_SCORE)) regions at
    /// score >= s. Real: K_RATIO times that (a K difference, a constant ratio) plus related
    /// pairs coming in as a ramp from RAMP_START, RAMP_HEIGHT e^(-RAMP_DECAY (s - RAMP_START))
    /// per bin. A step test would not see a ramp; the ratio test stops two bins after it
    /// starts and the slope is read below it, 4% low from the two ramp bins inside the
    /// window.
    #[test]
    fn test_fit_scores_with_reference_stops_where_the_ratio_rises() {
        const N_SEED: f64 = 400_000.0;
        const SEED_SCORE: i64 = 12;
        const TRUE_LAMBDA: f64 = 0.4;
        const K_RATIO: f64 = 1.5;
        const RAMP_START: i64 = 24;
        const RAMP_HEIGHT: f64 = 150.0;
        const RAMP_DECAY: f64 = 0.05;
        let at_least = |x: i64, scale: f64| {
            (scale * N_SEED * (-TRUE_LAMBDA * (x - SEED_SCORE) as f64).exp()).round() as i64
        };
        let (mut reference, mut real, mut plain) = (Vec::new(), Vec::new(), Vec::new());
        for s in SEED_SCORE..=44 {
            let per_bin =
                |scale: f64| (at_least(s, scale) - at_least(s + 1, scale)).max(0) as usize;
            reference.extend(std::iter::repeat_n(s as f64, per_bin(1.0)));
            plain.extend(std::iter::repeat_n(s as f64, per_bin(K_RATIO)));
            let ramp = if s >= RAMP_START {
                (RAMP_HEIGHT * (-RAMP_DECAY * (s - RAMP_START) as f64).exp()).round() as usize
            } else {
                0
            };
            real.extend(std::iter::repeat_n(s as f64, per_bin(K_RATIO) + ramp));
        }
        let fit = fit_scores_with_reference(&real, &reference).unwrap();
        assert_eq!((fit.bend_score, fit.score_lo, fit.score_hi), (Some(26), 18, 25), "{fit:?}");
        assert!((fit.lambda - 0.385).abs() < 0.005, "{}", fit.lambda);
        assert!(
            (fit.reference_lambda.unwrap() - TRUE_LAMBDA).abs() < 0.005,
            "{:?}",
            fit.reference_lambda
        );
        assert_eq!(fit.reference_survival[0], (SEED_SCORE, reference.len() as u64));

        // Without the ramp the fit runs to the count floor; the constant ratio is harmless.
        let fit = fit_scores_with_reference(&plain, &reference).unwrap();
        assert_eq!((fit.bend_score, fit.score_lo, fit.score_hi), (None, 26, 33), "{fit:?}");
        assert!((fit.lambda - TRUE_LAMBDA).abs() < 0.005, "{}", fit.lambda);
    }

    #[test]
    fn test_bin_counts() {
        assert_eq!(
            bin_counts(&[12.0, 12.0, 13.0, 15.5, 15.0]),
            vec![(12, 2), (13, 1), (14, 0), (15, 2)]
        );
    }

    #[test]
    fn test_fit_scores_needs_enough_bins() {
        // Three bins with at least MIN_BIN_COUNT regions after the skipped seed bin: fewer
        // than MIN_FIT_POINTS.
        let mut scores = vec![12.0; 100];
        scores.extend(vec![13.0; 60]);
        scores.extend(vec![14.0; 40]);
        scores.extend(vec![15.0; MIN_BIN_COUNT as usize + 1]);
        assert_eq!(fit_scores(&scores), None);
        assert_eq!(fit_scores(&[]), None);
    }

    /// Counts that grow with the score describe no exponential tail; neither fit will read
    /// a lambda off them.
    #[test]
    fn test_fits_refuse_a_rising_curve() {
        let mut rising = vec![12.0; 100];
        for (score, n) in [(13.0, 30), (14.0, 40), (15.0, 50), (16.0, 60)] {
            rising.extend(std::iter::repeat_n(score, n));
        }
        assert_eq!(fit_scores(&rising), None);
        let mut flat = vec![12.0; 100];
        for score in [13.0, 14.0, 15.0, 16.0] {
            flat.extend(std::iter::repeat_n(score, 45));
        }
        assert_eq!(fit_scores_with_reference(&rising, &flat), None);
        // The same counts falling with the score fit fine. Least squares through ln(count)
        // at scores 13..=16 gives slope -(3 ln(60/30) + ln(50/40)) / 10 = -ln(10) / 10.
        let mut falling = vec![12.0; 100];
        for (score, n) in [(13.0, 60), (14.0, 50), (15.0, 40), (16.0, 30)] {
            falling.extend(std::iter::repeat_n(score, n));
        }
        assert!((fit_scores(&falling).unwrap().lambda - 10f64.ln() / 10.0).abs() < 1e-12);
    }

    /// The window is stored in bins; `x_range` gives it back in nats. Bins 15..=22 of width
    /// 0.5 cover x from 7.5 up to but not including 11.5.
    #[test]
    fn test_calibration_reads_its_fit_window_in_x() {
        let fit = KaCalibration {
            score_lo: 15,
            score_hi: 22,
            bin_width: 0.5,
            r_database: 0.806,
            ..KaCalibration::placeholder(ExtensionScoring::default(), 0.0115)
        };
        assert_eq!(fit.n_fit_points(), 8);
        assert_eq!(fit.x_range(), (7.5, 11.5));
        assert_eq!(fit.ka_params(), KaParams { k: 0.0115, r_database: 0.806 });
    }

    #[test]
    fn test_decoy_null_names_match_the_flag_values() {
        assert_eq!(DecoyNull::Database.to_string(), "database");
        assert_eq!(DecoyNull::Shuffled.to_string(), "shuffled");
        assert_eq!(DecoyNull::ShuffledDipeptide.to_string(), "shuffled-dipeptide");
        assert_eq!(DecoyNull::Reversed.to_string(), "reversed");
    }

    #[test]
    fn test_sampler_is_deterministic_and_distinct() {
        let a = SplitMix64::new(7).sample_indices(1000, 25);
        let b = SplitMix64::new(7).sample_indices(1000, 25);
        assert_eq!(a, b);
        assert_eq!(a.len(), 25);
        assert!(a.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(SplitMix64::new(7).sample_indices(10, 25), (0..10).collect::<Vec<_>>());
    }

    fn dipeptide_counts(s: &str) -> std::collections::BTreeMap<(u8, u8), usize> {
        let mut counts = std::collections::BTreeMap::new();
        for w in s.as_bytes().windows(2) {
            *counts.entry((w[0], w[1])).or_insert(0) += 1;
        }
        counts
    }

    /// The one sequence in a FASTA file, header dropped and lines joined.
    fn read_single_fasta(path: &str) -> String {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|l| !l.starts_with('>'))
            .collect::<Vec<_>>()
            .concat()
    }

    /// Human BCL2 (P10415, `TEST_BLC2_FASTA`). Same length, same first and last residue,
    /// every dipeptide as often as before, a different order, and the same order again from
    /// the same seed.
    #[test]
    fn test_dipeptide_shuffle_keeps_every_dipeptide_count() {
        let bcl2 = read_single_fasta(TEST_BLC2_FASTA);
        let bcl2 = bcl2.as_str();
        assert_eq!(bcl2.len(), 239, "BCL2_HUMAN is 239 residues");
        let shuffled = dipeptide_shuffle(bcl2, &mut SplitMix64::new(3));
        assert_eq!(shuffled.len(), bcl2.len());
        assert_eq!(dipeptide_counts(&shuffled), dipeptide_counts(bcl2));
        assert_eq!(shuffled.as_bytes()[0], b'M');
        assert_eq!(shuffled.as_bytes()[shuffled.len() - 1], b'K');
        assert_ne!(shuffled, bcl2);
        assert_eq!(shuffled, dipeptide_shuffle(bcl2, &mut SplitMix64::new(3)));
        assert_ne!(shuffled, dipeptide_shuffle(bcl2, &mut SplitMix64::new(4)));
        // A plain shuffle of the same sequence does not keep the dipeptides.
        assert_ne!(
            dipeptide_counts(&make_decoy(bcl2, DecoyNull::Shuffled, &mut SplitMix64::new(3))),
            dipeptide_counts(bcl2)
        );
        // Too short to rearrange, or all one residue: returned as is.
        assert_eq!(dipeptide_shuffle("MK", &mut SplitMix64::new(1)), "MK");
        assert_eq!(dipeptide_shuffle("AAAAAA", &mut SplitMix64::new(1)), "AAAAAA");
    }

    #[test]
    fn test_make_decoy() {
        let mut rng = SplitMix64::new(1);
        assert_eq!(make_decoy("MKTAYIAK", DecoyNull::Database, &mut rng), "MKTAYIAK");
        assert_eq!(make_decoy("MKTAYIAK", DecoyNull::Reversed, &mut rng), "KAIYATKM");
        assert_eq!(
            make_decoy("MKTAYIAK", DecoyNull::Shuffled, &mut SplitMix64::new(1)),
            "YATKIAMK"
        );
    }
}
