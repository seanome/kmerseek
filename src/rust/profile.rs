//! Per-step timers for `kmerseek search`, switched on by setting `KMERSEEK_PROFILE`.
//!
//! Each step adds its elapsed time to one global counter. Steps that run inside rayon's
//! parallel loops are summed over all threads, so their totals are thread-seconds and can
//! exceed the run's wall time. With `KMERSEEK_PROFILE` unset, a timer costs one cached bool
//! check and reads no clock.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;
use std::time::Instant;

/// One timed part of a search, in the order a query passes through them.
#[derive(Clone, Copy, Debug)]
pub enum Step {
    /// Opening the index and reading its target list and inverted index into memory.
    LoadIndex,
    /// `prepare_query`: the query's k-mer set, TF-IDF and position prefix sums.
    PrepareQuery,
    /// Looking up each query k-mer in the inverted index and collecting candidate targets.
    CandidateLookup,
    /// Reading a candidate's signature from RocksDB and deserializing it (first use only).
    SignatureLoad,
    /// The target's k-mer set, its intersection with the query's, and the cheap filters.
    Intersection,
    /// Similarity statistics and exact matched regions (`find_matched_regions`).
    MatchedRegions,
    /// Extending regions past their exact seeds (only with extension on).
    Extension,
    /// The per-region Poisson score, enrichment and IDF sums.
    RegionScoring,
    /// Karlin-Altschul bits, E-values and region chaining (only with extension on).
    KarlinAltschul,
    /// Database-wide metrics filled into a kept result (coverage score, k-mer frequencies).
    ResultMetrics,
    /// Sorting one query's results.
    SortResults,
    /// Writing result rows to the CSV.
    WriteCsv,
}

const STEPS: [Step; 12] = [
    Step::LoadIndex,
    Step::PrepareQuery,
    Step::CandidateLookup,
    Step::SignatureLoad,
    Step::Intersection,
    Step::MatchedRegions,
    Step::Extension,
    Step::RegionScoring,
    Step::KarlinAltschul,
    Step::ResultMetrics,
    Step::SortResults,
    Step::WriteCsv,
];

static NANOS: [AtomicU64; STEPS.len()] = [const { AtomicU64::new(0) }; STEPS.len()];
static CALLS: [AtomicU64; STEPS.len()] = [const { AtomicU64::new(0) }; STEPS.len()];

/// True when `KMERSEEK_PROFILE` is set; read once per process.
pub fn enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var_os("KMERSEEK_PROFILE").is_some())
}

/// Adds the time from its creation to its drop to one step's total.
pub struct Timer {
    step: Step,
    start: Option<Instant>,
}

impl Timer {
    /// Starts timing `step`. Reads no clock when profiling is off.
    pub fn start(step: Step) -> Self {
        Self { step, start: enabled().then(Instant::now) }
    }
}

impl Drop for Timer {
    fn drop(&mut self) {
        if let Some(start) = self.start {
            let i = self.step as usize;
            NANOS[i].fetch_add(start.elapsed().as_nanos() as u64, Ordering::Relaxed);
            CALLS[i].fetch_add(1, Ordering::Relaxed);
        }
    }
}

/// Prints every step's total seconds and call count to stderr, as tab-separated lines that
/// start with `profile`, so they can be pulled out of a log with `grep ^profile`.
pub fn report() {
    if !enabled() {
        return;
    }
    eprintln!("profile\tstep\tseconds\tcalls");
    for step in STEPS {
        let i = step as usize;
        let seconds = NANOS[i].load(Ordering::Relaxed) as f64 / 1e9;
        let calls = CALLS[i].load(Ordering::Relaxed);
        eprintln!("profile\t{step:?}\t{seconds:.3}\t{calls}");
    }
}
