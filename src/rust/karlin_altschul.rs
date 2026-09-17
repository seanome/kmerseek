//! Karlin-Altschul lambda and K for the E-value of an extended region, fitted on the index
//! itself: the closed forms for an unseeded ungapped search, and the survival-curve fit on
//! calibration searches that kmerseek actually uses.
//!
//! E = K m n e^(-lambda S) counts the regions with score >= S expected between an
//! unrelated query of m residues and a database of n residues. The closed forms assume
//! independent positions and count every high-scoring segment. kmerseek only finds a
//! region that contains an exact k-mer seed, extends it with an X-drop, and works on real
//! proteins whose hydrophobic runs and helix and strand periodicity are not independent
//! positions. So lambda and K are read off the search itself: a few hundred database
//! sequences are searched against the index, and ln(count of regions with score S) is a
//! straight line in S with slope -lambda and intercept ln(K L N (1 - e^-lambda)), L the
//! calibration residues and N the database residues (Altschul & Gish 1996; Pearson 1998).
//! Homologs and hydrophobic runs lift the counts at high scores; the fit stops below that.
//! The line is fitted to the count at each score, not the count at or above it, because a
//! plateau of homolog hits far up the score axis adds a constant to every survival count
//! below it and would flatten the slope there too.

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
    /// in; related pairs are cut off as the upward bend of the survival curve (Altschul's
    /// tutorial, approach i; Collins et al. 1988; Pearson 1998).
    Database,
    /// Each query is a database sequence with its residues shuffled: the independent-letter
    /// model BLAST's tables are fitted on (Altschul & Gish 1996).
    Shuffled,
    /// Each query is a database sequence read back to front.
    Reversed,
}

impl std::fmt::Display for DecoyNull {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DecoyNull::Database => write!(f, "database"),
            DecoyNull::Reversed => write!(f, "reversed"),
            DecoyNull::Shuffled => write!(f, "shuffled"),
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
    /// (a of the database against itself), and the closed-form lambda at that a.
    pub match_probability: f64,
    pub lambda_analytic: f64,
    /// Slope and intercept of the survival line: the lambda and K a search uses.
    pub lambda: f64,
    pub k: f64,
    /// Score bins the line was fitted on, inclusive.
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
}

impl KaCalibration {
    /// How far the fitted lambda sits from the independent-positions lambda at the
    /// database's own composition. A search scales every pair's closed-form lambda by
    /// this, so a pair with the database's composition gets the fitted slope and a pair of
    /// two hydrophobic sequences still gets lambda 0.
    pub fn lambda_scale(&self) -> f64 {
        if self.lambda_analytic > 0.0 {
            self.lambda / self.lambda_analytic
        } else {
            1.0
        }
    }

    /// Number of score bins the line was fitted on.
    pub fn n_fit_points(&self) -> i64 {
        self.score_hi - self.score_lo + 1
    }
}

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
/// Bins with fewer than `MIN_BIN_COUNT` regions are ignored. The lowest bin is skipped too:
/// it holds every bare seed and says nothing about extension. Starting from the next
/// `MIN_FIT_POINTS` bins, the line is extended one bin at a time upward and a bin joins
/// while its ln count is within `BEND_SIGMAS / sqrt(count) + BEND_SLACK` above the line
/// fitted so far (below it is fine: the seed requirement makes the true curve concave).
/// The first bin above that is where homologs start to show and where the fit stops. The
/// line is then read off the highest `FIT_WINDOW` bins that made it in, the ones nearest
/// the scores that decide a hit. None when fewer than `MIN_FIT_POINTS` bins qualify.
pub fn fit_scores(scores: &[f64]) -> Option<ScoreFit> {
    let usable: Vec<(i64, u64)> =
        bin_counts(scores).into_iter().skip(1).filter(|&(_, c)| c >= MIN_BIN_COUNT).collect();
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
    })
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
    }
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

    #[test]
    fn test_sampler_is_deterministic_and_distinct() {
        let a = SplitMix64::new(7).sample_indices(1000, 25);
        let b = SplitMix64::new(7).sample_indices(1000, 25);
        assert_eq!(a, b);
        assert_eq!(a.len(), 25);
        assert!(a.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(SplitMix64::new(7).sample_indices(10, 25), (0..10).collect::<Vec<_>>());
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
