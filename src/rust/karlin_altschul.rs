//! Karlin-Altschul K, and a check on lambda, for the E-value of an extended region, fitted
//! on the index itself.
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
//! and intercept ln(K L N (1 - e^-w)) (Altschul & Gish 1996; Pearson 1998). The fitted
//! slope is the factor every pair's lambda is multiplied by at search time; 1 means the
//! closed form holds. Relatives lift the counts at high x; the fit stops below that, where
//! the real curve starts to rise relative to the same queries shuffled. The fit goes to the
//! count at each x, not the count at or above it, because a plateau of relatives far up the
//! axis adds a constant to every survival count below it and would flatten the slope.
//!
//! Fitting x rather than the raw score S matters on a proteome: Swiss-Prot holds pairs of
//! membrane and low-complexity proteins whose raw scores run to 60 and beyond with a slope
//! near 0.1, while ordinary pairs fall at 0.45. One line cannot serve both. In x each pair
//! is already on its own scale, and those pairs sit at x = 0.

use serde::{Deserialize, Serialize};

use crate::search::karlin_altschul_lambda;

/// Karlin-Altschul K for +1 / -penalty scoring with match probability `a`, counting every
/// high-scoring segment of an ungapped comparison (no seed requirement, no X-drop).
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

/// Which sequences are searched to fit lambda and K.
///
/// Measured on 300 human BCL2-related proteins at hp k=12, penalty 2, X-drop 8, 200
/// calibration queries each (local slope of ln count vs score, per 4-score window):
///
/// - `Database`: 0.35 to 0.36 up to score 20, then the line bends upward into 1,800
///   regions of family hits reaching score 1,688. The bend is what the fit cuts off.
/// - `Shuffled`: 0.37 rising to 0.44 between scores 12 and 32, no bend. Straight enough,
///   but shuffling removes the hydrophobic runs and the periodicity real unrelated proteins
///   have, so the null is too easy and its E-values too small.
/// - `Reversed`: 0.36 falling to 0.18 between scores 12 and 48. A helix or a strand reads
///   much the same backwards in a hydrophobic/polar alphabet, so reversed family members
///   still hit the query one element at a time; the leak sits at the scores the fit needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, clap::ValueEnum)]
pub enum DecoyNull {
    /// Database sequences as they are, searched against the index. Everything real stays
    /// in; related pairs are cut off where the curve starts to rise relative to the same
    /// queries shuffled, which are searched alongside (Altschul's tutorial, approach i;
    /// Collins et al. 1988; Pearson 1998). Costs two calibration searches.
    Database,
    /// Each query is a database sequence with its residues shuffled: the independent-letter
    /// model BLAST's tables are fitted on (Altschul & Gish 1996).
    Shuffled,
    /// Each query is a database sequence shuffled so that every dipeptide (each pair of
    /// neighbouring residues) occurs as often as in the original (Altschul & Erickson 1985;
    /// sampled as a random Eulerian path, Kandel et al. 1996, the uShuffle k = 2 method).
    /// Keeps the rate at which a hydrophobic residue follows a hydrophobic one, and with it
    /// the lengths of hydrophobic runs, which a plain shuffle destroys.
    ShuffledDipeptide,
    /// Each query is a database sequence read back to front.
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

/// One fitted (lambda, K), stored in the index under `ka_calibration` and looked up at
/// search time by (mismatch_penalty, xdrop). The seed length and alphabet are the index's.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct KaCalibration {
    pub mismatch_penalty: f64,
    pub xdrop: f64,
    pub null: DecoyNull,
    /// Seed of the query sampler, so the fit can be reproduced.
    pub seed: u64,
    /// Calibration queries searched, and their residues added up (L in ln(K L N)).
    pub n_queries: usize,
    pub query_residues: u64,
    /// The database size the E-value uses for n: the index's k-mer count stands in for its
    /// residue count.
    pub database_kmers: u64,
    /// Regions the calibration queries produced, at any score.
    pub n_regions: usize,
    /// Chance that two positions drawn from the sampled database sequences share a class
    /// (a of the database against itself), and the closed-form lambda at that a, for
    /// reporting.
    pub match_probability: f64,
    pub lambda_analytic: f64,
    /// Minus the slope of ln(count) against x = lambda_pair S, per nat. 1 means the
    /// closed-form per-pair lambda has the right scale; a search multiplies every pair's
    /// lambda by this.
    pub slope: f64,
    pub k: f64,
    /// Width of one x bin in nats (`BIN_WIDTH`); `score_lo`, `score_hi`, `bend_score` and
    /// the survival curves are in bins, so bin b covers x in [b w, (b + 1) w).
    pub bin_width: f64,
    /// Bins the line was fitted on, inclusive.
    pub score_lo: i64,
    pub score_hi: i64,
    /// First score bin above the fit that sat above the line: the start of the homolog
    /// bend. None when the fit ran out of counts first.
    pub bend_score: Option<i64>,
    /// Root mean square of the fit residuals in ln count.
    pub rms_residual: f64,
    /// Regions with score >= s, for every s from the smallest score seen up to the largest,
    /// for plotting the curve the fit was read from.
    pub survival: Vec<(i64, u64)>,
    /// For the `Database` null: the same queries shuffled, searched the same way, and the
    /// slope of that curve over the fit window. The ratio of the two curves is what decides
    /// where the fit stops (`fit_scores_with_reference`). Empty and None for other nulls.
    pub reference_survival: Vec<(i64, u64)>,
    pub reference_lambda: Option<f64>,
    /// How the reference queries were made (`Shuffled` or `ShuffledDipeptide`).
    pub reference: Option<DecoyNull>,
}

impl KaCalibration {
    /// The factor a search multiplies every pair's closed-form lambda by: the fitted slope.
    pub fn r_database(&self) -> f64 {
        self.slope
    }

    /// The fit window in nats of x = lambda_pair S.
    pub fn x_range(&self) -> (f64, f64) {
        (self.score_lo as f64 * self.bin_width, (self.score_hi + 1) as f64 * self.bin_width)
    }

    /// Number of score bins the line was fitted on.
    pub fn n_fit_points(&self) -> i64 {
        self.score_hi - self.score_lo + 1
    }
}

/// Width of one bin of the normalised score x = lambda_pair S, in nats. Half a nat is about
/// one raw score unit at lambda 0.45.
pub const BIN_WIDTH: f64 = 0.5;

/// Fewest regions a score bin needs to enter the fit; below this the Poisson noise in
/// ln count (about 1 / sqrt(count)) is larger than the effects being fitted.
pub const MIN_BIN_COUNT: u64 = 30;

/// The fit uses at most this many score bins, the highest ones below the homolog excess,
/// so the slope is read as close to the decision tail as that excess allows.
pub const FIT_WINDOW: i64 = 8;

/// Fewest bins a fit is accepted on.
pub const MIN_FIT_POINTS: i64 = 4;

/// A bin whose ln count sits more than this many Poisson standard deviations above the
/// line fitted to the bins below it starts the homolog excess.
const BEND_SIGMAS: f64 = 2.0;

/// Slack added to the bend test, in ln count, so a bin one region above the line at large
/// counts does not end the fit.
const BEND_SLACK: f64 = 0.05;

/// Slope and intercept of ln(regions with score S) against S, and where it was read.
#[derive(Debug, Clone, PartialEq)]
pub struct ScoreFit {
    pub lambda: f64,
    /// ln(K L N (1 - e^-lambda)): the line's value at score 0.
    pub ln_intercept: f64,
    pub score_lo: i64,
    pub score_hi: i64,
    pub bend_score: Option<i64>,
    pub rms_residual: f64,
    pub survival: Vec<(i64, u64)>,
    /// Shuffled-sequence reference curve and its slope over the same window, when the fit
    /// was censored against one (`fit_scores_with_reference`); empty and None otherwise.
    pub reference_survival: Vec<(i64, u64)>,
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

/// Fit ln(regions with score S) against S, stopping below the homolog excess.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_k_theory_matches_textbook_plus_minus_one_case() {
        // Steps +1 with p = 0.3 and -1 with q = 0.7: K = (q - p)^2 / q (Ewens & Grant,
        // Statistical Methods in Bioinformatics, the BLAST random-walk chapter).
        let k = karlin_altschul_k_theory(0.3, 1.0).unwrap();
        assert!((k - 0.16 / 0.7).abs() < 1e-12, "{k}");
    }

    #[test]
    fn test_k_theory_balanced_hp_at_penalty_two() {
        // a = 0.5, C = 2: lambda = ln(golden ratio), e^lambda = phi, e^-2lambda = 1/phi^2.
        // E[X e^(lambda X)] = phi/2 - 1/phi^2 and (1 - e^-lambda) = 1 - 1/phi.
        let phi = (1.0 + 5f64.sqrt()) / 2.0;
        let expected = (phi / 2.0 - 1.0 / (phi * phi)) * (1.0 - 1.0 / phi);
        let k = karlin_altschul_k_theory(0.5, 2.0).unwrap();
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

    /// Exact exponential counts: 4000 e^(-0.4 (s - 12)) regions at score >= s for
    /// s = 12..=30, so 4000 (1 - e^-0.4) e^(-0.4 (s - 12)) at score s. Then 400 extra
    /// regions at score 27, a homolog bump. The fit reads the slope off the bins below the
    /// bump, to within rounding of the counts to whole regions, and the bump does not reach
    /// the bins below it the way it would on the survival curve.
    #[test]
    fn test_fit_scores_recovers_slope_and_stops_at_homolog_excess() {
        let at_least = |x: i64| (4000.0 * (-0.4 * (x - 12) as f64).exp()).round() as i64;
        let mut scores = Vec::new();
        for s in 12..=30 {
            let exactly = (at_least(s) - at_least(s + 1)).max(0) as usize;
            scores.extend(std::iter::repeat_n(s as f64, exactly));
        }
        let clean = fit_scores(&scores).unwrap();
        // Bins 13..=21 hold at least 30 regions; the top 8 of them are the window.
        assert_eq!((clean.score_lo, clean.score_hi, clean.bend_score), (14, 21, None));
        assert!((clean.lambda - 0.4).abs() < 0.01, "{}", clean.lambda);
        let expected_intercept = (4000.0 * (1.0 - (-0.4f64).exp())).ln() + 0.4 * 12.0;
        assert!((clean.ln_intercept - expected_intercept).abs() < 0.1, "{}", clean.ln_intercept);
        assert!(clean.rms_residual < 0.03, "{}", clean.rms_residual);
        assert_eq!(clean.survival[0], (12, scores.len() as u64));

        scores.extend(std::iter::repeat_n(27.0, 400));
        let bent = fit_scores(&scores).unwrap();
        assert_eq!((bent.score_lo, bent.score_hi, bent.bend_score), (14, 21, Some(27)));
        assert!((bent.lambda - clean.lambda).abs() < 1e-12);
        assert_eq!(bent.survival[0], (12, scores.len() as u64));
    }

    /// Reference: 400,000 e^(-0.4 (s - 12)) regions at score >= s. Real: 1.5x that (a K
    /// difference, a constant ratio) plus relatives coming in as a ramp from score 24,
    /// 150 e^(-0.05 (s - 24)) per bin. A step test would not see a ramp; the ratio test
    /// stops two bins after it starts and the slope is read below it, 4% low from the two
    /// ramp bins inside the window.
    #[test]
    fn test_fit_scores_with_reference_stops_where_the_ratio_rises() {
        let at_least = |x: i64, scale: f64| {
            (scale * 400_000.0 * (-0.4 * (x - 12) as f64).exp()).round() as i64
        };
        let (mut reference, mut real, mut plain) = (Vec::new(), Vec::new(), Vec::new());
        for s in 12..=44 {
            let per_bin =
                |scale: f64| (at_least(s, scale) - at_least(s + 1, scale)).max(0) as usize;
            reference.extend(std::iter::repeat_n(s as f64, per_bin(1.0)));
            plain.extend(std::iter::repeat_n(s as f64, per_bin(1.5)));
            let ramp = if s >= 24 {
                (150.0 * (-0.05 * (s - 24) as f64).exp()).round() as usize
            } else {
                0
            };
            real.extend(std::iter::repeat_n(s as f64, per_bin(1.5) + ramp));
        }
        let fit = fit_scores_with_reference(&real, &reference).unwrap();
        assert_eq!((fit.bend_score, fit.score_lo, fit.score_hi), (Some(26), 18, 25), "{fit:?}");
        assert!((fit.lambda - 0.385).abs() < 0.005, "{}", fit.lambda);
        assert!((fit.reference_lambda.unwrap() - 0.4).abs() < 0.005, "{:?}", fit.reference_lambda);
        assert_eq!(fit.reference_survival[0], (12, reference.len() as u64));

        // Without the ramp the fit runs to the count floor; the constant ratio is harmless.
        let fit = fit_scores_with_reference(&plain, &reference).unwrap();
        assert_eq!((fit.bend_score, fit.score_lo, fit.score_hi), (None, 26, 33), "{fit:?}");
        assert!((fit.lambda - 0.4).abs() < 0.005, "{}", fit.lambda);
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
        // Three usable bins after the skipped seed bin: too few.
        let mut scores = vec![12.0; 100];
        scores.extend(vec![13.0; 60]);
        scores.extend(vec![14.0; 40]);
        scores.extend(vec![15.0; 31]);
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

    #[test]
    fn test_calibration_reads_its_fit_window_in_x() {
        let fit = KaCalibration {
            mismatch_penalty: 2.0,
            xdrop: 8.0,
            null: DecoyNull::Database,
            seed: 1,
            n_queries: 25,
            query_residues: 9288,
            database_kmers: 8340,
            n_regions: 9838,
            match_probability: 0.5,
            lambda_analytic: 0.481,
            slope: 0.806,
            k: 0.0115,
            bin_width: BIN_WIDTH,
            score_lo: 15,
            score_hi: 22,
            bend_score: None,
            rms_residual: 0.086,
            survival: vec![(12, 9838)],
            reference_survival: Vec::new(),
            reference_lambda: None,
            reference: Some(DecoyNull::ShuffledDipeptide),
        };
        // The slope of ln(count) against lambda_pair S is the lambda scale itself; the
        // window is bins 15..=22 of width 0.5, so x from 7.5 up to but not including 11.5.
        assert_eq!(fit.lambda_scale(), 0.806);
        assert_eq!(fit.n_fit_points(), 8);
        assert_eq!(fit.x_range(), (7.5, 11.5));
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

    /// Human BCL2 (P10415, tests/testdata/fasta/bcl2.fasta). Same length, same first and
    /// last residue, every dipeptide as often as before, a different order, and the same
    /// order again from the same seed.
    #[test]
    fn test_dipeptide_shuffle_keeps_every_dipeptide_count() {
        let bcl2 = "MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK";
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
