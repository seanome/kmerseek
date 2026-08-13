use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::path::Path;

use anyhow::Result;
use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use statrs::distribution::{DiscreteCDF, Poisson};

use crate::errors::IndexResult;
use crate::index::ProteomeIndex;
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
    /// Maximum region-scoped Poisson p-value required to keep a match, applied to the best
    /// region in the pair. A heuristic cutoff on a ranking score, not a statistically
    /// calibrated significance threshold; see `MatchedRegion::poisson_pvalue`.
    pub max_region_pvalue: f64,
}

impl Default for SearchFilters {
    /// Accepts every result `compare()` produces. Note the p-value caps must be `f64::INFINITY`,
    /// not 1.0: `query_poisson_pvalue` is 1.0 whenever there's no database frequency context
    /// (e.g. `set_query_frequencies` was never called), which is common, so a cap of 1.0
    /// combined with `compare()`'s strict-less-than keep check would wrongly reject those
    /// results.
    fn default() -> Self {
        Self {
            threshold: 0.0,
            min_shared_kmers: 0,
            max_query_pvalue: f64::INFINITY,
            max_region_pvalue: f64::INFINITY,
        }
    }
}

impl SearchFilters {
    /// A pair is kept when either scope clears its cap. Both are not required.
    ///
    /// Requiring both would bring back the problem this PR fixes: a real sub-protein domain
    /// match diluted into insignificance by the rest of the protein. BCL2/CED9 at k=15 has a
    /// whole-query p-value of 0.99 and a region p-value of 0.0007. Requiring both to pass
    /// would discard the sub-protein domain match that region scoring exists to surface.
    fn pvalues_pass(&self, query_pvalue: f64, best_region_pvalue: Option<f64>) -> bool {
        let query_passes = query_pvalue < self.max_query_pvalue;
        let region_passes =
            best_region_pvalue.is_some_and(|pvalue| pvalue < self.max_region_pvalue);
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
    pub region_poisson_pvalue: f64,
    pub region_enrichment: f64,
}

impl SearchResultCsv {
    /// Create a CSV row from a SearchResult and a MatchedRegion
    /// WHY: Every CSV row must have matched region data. This method combines the SearchResult
    /// similarity metrics with a specific matched region to create one CSV row. Each SearchResult
    /// will produce multiple CSV rows (one per matched region), with all similarity metrics
    /// repeated for each region.
    pub fn from_result_and_region(result: &SearchResult, region: &MatchedRegion) -> Self {
        // See the matching debug_assert in ProteinSearcher::compare: a region shorter than
        // ksize should never exist, and saturating_sub would otherwise hide that as a silent
        // region_n_shared_kmers = 1 instead of a loud failure.
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
            region_n_shared_kmers: region.length.saturating_sub(result.ksize) + 1,
            region_expected_shared_kmers: region.expected_shared_kmers,
            region_poisson_pvalue: region.poisson_pvalue,
            region_enrichment: region.enrichment,
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

    // One of "protein", "dayhoff", or "hp"
    pub moltype: MolType,

    /// Length of the match
    pub length: u32,

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

    /// Poisson p-value scoped to this region: the probability of seeing at least `n_shared`
    /// shared k-mers if matches happened at random, given the rate `expected_shared_kmers`.
    /// Uses the region's own k-mer count instead of the whole protein's. 1.0 without DB
    /// context.
    ///
    /// Treat this as a heuristic score for ranking candidate regions against each other, not
    /// as a calibrated significance estimate. Two separate problems keep it from being a real
    /// p-value:
    ///
    /// 1. `n_shared` is not an independent observation. It is `region_length - ksize + 1`,
    ///    arithmetic on the region's own length, and the region's length is exactly what
    ///    `find_matched_regions` chose by keeping the longest gapless run of shared k-mers. The
    ///    test is being applied to the same quantity that defined the region, which is close to
    ///    circular. Multiplying by `region_search_space` (how many positions a region could
    ///    have started at) and `db_n_targets` (how many targets were searched) into an E-value,
    ///    `poisson_pvalue * region_search_space * db_n_targets`, corrects for having picked the
    ///    best-looking window out of many candidate windows. It does not fix this problem.
    ///
    /// 2. The k-mers being counted overlap by `ksize - 1` residues, so they are not independent
    ///    trials the way the Poisson model assumes. Five overlapping k-mers spanning a single
    ///    19-residue stretch are closer to one piece of evidence, observed five times, than to
    ///    five separate pieces of evidence. Treating them as independent understates how likely
    ///    a run this long is to appear by chance.
    ///
    /// A properly calibrated version of this statistic would model the length of the longest
    /// gapless run directly (an extreme-value distribution, the same kind of model behind
    /// BLAST's E-values), account for the k-mer overlap, and be checked empirically against
    /// a decoy database. None of that is implemented here. Use this p-value to prioritize which
    /// regions to look at first, not to make a significance claim about any single region.
    pub poisson_pvalue: f64,

    /// Fold-enrichment scoped to this region: n_shared / expected_shared_kmers. 0.0 without DB
    /// context or when expected_shared_kmers is 0.
    pub enrichment: f64,
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

/// Observed over expected. 0.0 when there is no expectation to divide by, which reads as
/// "not computable" rather than the +inf the division would produce.
fn fold_enrichment(observed: u32, expected: f64) -> f64 {
    if expected > 0.0 {
        observed as f64 / expected
    } else {
        0.0
    }
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
}

impl SearchStats {
    /// Calculate search statistics from a proteome index
    pub fn from_index(index: &ProteomeIndex) -> Self {
        let signatures = index.get_signatures();
        let total_signatures = signatures.len();

        // Count k-mer frequencies across all signatures
        let mut kmer_frequencies: HashMap<u64, usize> = HashMap::new();

        for signature in signatures.iter() {
            let mins = signature.value().signature().minhash.mins();
            for min in mins {
                *kmer_frequencies.entry(min).or_insert(0) += 1;
            }
        }

        // Calculate IDF values
        let idf: HashMap<u64, f64> = kmer_frequencies
            .iter()
            .map(|(&kmer, &freq)| {
                let idf_value = (total_signatures as f64 / freq as f64).ln();
                (kmer, idf_value)
            })
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
}

impl ProteinSearcher {
    /// Create a new protein searcher from an index
    pub fn new(index: ProteomeIndex) -> Self {
        let stats = SearchStats::from_index(&index);
        let (target_list, inverted_index) = Self::build_search_structures(&index);
        let db_n_kmers = stats.kmer_frequencies.values().sum();
        Self {
            index,
            stats,
            target_list,
            inverted_index,
            sig_cache: DashMap::new(),
            query_kmer_frequencies: None,
            total_queries: 0,
            db_n_kmers,
        }
    }

    /// Load a searcher from a saved index.
    ///
    /// Fast path: if the index was built with a recent version of kmerseek (which saves a
    /// pre-built search cache), this method opens the DB without loading all signatures into
    /// memory. Signatures are then loaded on demand during search via `get_signature_by_md5()`.
    ///
    /// Slow path (backward compat): for older databases without a search cache, falls back to
    /// loading all signatures into memory and building the inverted index at startup.
    pub fn load<P: AsRef<Path>>(path: P) -> IndexResult<Self> {
        // Open DB minimally: read metadata only, leave signatures DashMap empty
        let index = ProteomeIndex::open_for_search(&path)?;

        // Fast path: pre-built search cache exists - no need to load all signatures
        if let Some(cache) = index.load_search_cache()? {
            let total_signatures = cache.target_list.len();
            let idf: HashMap<u64, f64> = cache
                .kmer_frequencies
                .iter()
                .map(|(&kmer, &freq)| {
                    let idf_value = (total_signatures as f64 / freq as f64).ln();
                    (kmer, idf_value)
                })
                .collect();
            let stats =
                SearchStats { total_signatures, idf, kmer_frequencies: cache.kmer_frequencies };
            let db_n_kmers = stats.kmer_frequencies.values().sum();
            eprintln!(
                "Loaded search cache: {} targets, {} k-mers indexed",
                total_signatures,
                cache.inverted_index.len()
            );
            return Ok(Self {
                index,
                stats,
                target_list: cache.target_list,
                inverted_index: cache.inverted_index,
                sig_cache: DashMap::new(),
                query_kmer_frequencies: None,
                total_queries: 0,
                db_n_kmers,
            });
        }

        // Slow path: old DB without search cache - load all signatures and build structures
        eprintln!(
            "No search cache found; loading all signatures (run `kmerseek index` to rebuild)"
        );
        index.load_state()?;
        let stats = SearchStats::from_index(&index);
        let (target_list, inverted_index) = Self::build_search_structures(&index);
        let db_n_kmers = stats.kmer_frequencies.values().sum();
        Ok(Self {
            index,
            stats,
            target_list,
            inverted_index,
            sig_cache: DashMap::new(),
            query_kmer_frequencies: None,
            total_queries: 0,
            db_n_kmers,
        })
    }

    /// Build an ordered target list and inverted k-mer index from the index.
    ///
    /// WHY: The inverted index maps each k-mer hash to the set of target signatures that
    /// contain it. This allows search_one to skip the vast majority of targets that share
    /// no k-mers with the query, reducing search from O(Q×T) to O(Q×candidates) where
    /// candidates << T for most real queries. Building this once at load time amortizes
    /// the cost across all subsequent searches.
    fn build_search_structures(index: &ProteomeIndex) -> (Vec<String>, HashMap<u64, Vec<u32>>) {
        let mut target_list: Vec<String> = Vec::new();
        let mut inverted_index: HashMap<u64, Vec<u32>> = HashMap::new();

        for entry in index.get_signatures().iter() {
            let idx = target_list.len() as u32;
            target_list.push(entry.key().clone());
            for min in entry.value().signature().minhash.mins() {
                inverted_index.entry(min).or_default().push(idx);
            }
        }

        (target_list, inverted_index)
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
        }
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
        candidate_set
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
            .collect()
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
        // Skip self-matches by comparing MD5 sums
        // WHY: In all-vs-all searches, we don't want to compare a signature against itself.
        // MD5 sum is a unique identifier for each signature, so comparing MD5 sums is the
        // most reliable way to detect self-matches. This is idiomatic Rust - we use early
        // returns to avoid unnecessary computation when we know the result will be invalid.
        if query.sketch.signature().md5sum == target.signature().md5sum {
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
        // The p-value checks cannot be hoisted up here: a pair is kept when either scope
        // clears, and the region scope isn't known until the regions exist. So p-value
        // filtering happens after the result is built, and pairs that fail the query scope now
        // pay for region-finding before being rejected.
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

        // Rescope the same Poisson test to each matched region individually, so a tight local
        // match doesn't get diluted by the whole protein's k-mer count.
        let ksize = query.sketch.protein_ksize() as usize;
        for region in result.matched_regions.iter_mut() {
            let lambda = self.region_expectation(query.sketch, region.start, region.end, ksize);
            // find_matched_regions never emits a region shorter than ksize (its length is
            // consecutive_count + ksize - 1, and consecutive_count >= 1), so this can't
            // actually underflow. debug_assert catches it loudly if that invariant is ever
            // broken, instead of saturating_sub silently turning a bug into n_shared = 1.
            debug_assert!(
                region.length >= ksize as u32,
                "region shorter than ksize: length={}, ksize={ksize}",
                region.length
            );
            let n_shared = region.length.saturating_sub(ksize as u32) + 1;

            region.expected_shared_kmers = lambda;
            region.poisson_pvalue = poisson_survival(n_shared, lambda);
            region.enrichment = fold_enrichment(n_shared, lambda);
        }

        // Either scope clearing its cap keeps the pair - see SearchFilters::pvalues_pass.
        let best_region_pvalue = result
            .matched_regions
            .iter()
            .map(|region| region.poisson_pvalue)
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        if !filters.pvalues_pass(query_poisson_pvalue, best_region_pvalue) {
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

    /// Expected shared k-mers by chance within one region: for every query k-mer whose start
    /// position falls inside the region, sum its database frequency (occurrences across the
    /// database divided by the number of signatures).
    ///
    /// The window is `[start, end - ksize + 1)`, not the region's full span: a k-mer belongs to
    /// the region only if it fits entirely inside. That is what makes the count equal
    /// `length - ksize + 1` at scaled=1.
    fn region_expectation(&self, query: &ProteinSketch, start: u32, end: u32, ksize: usize) -> f64 {
        let window_start = start as usize;
        // A k-mer at p covers [p, p + ksize), so it fits inside [start, end) only when
        // p <= end - ksize. A span shorter than ksize holds no whole k-mer at all - saturating
        // here would wrongly admit position 0.
        let window_end = match (end as usize).checked_sub(ksize) {
            Some(last_start) => last_start + 1,
            None => window_start,
        };
        query
            .kmer_positions()
            .iter()
            .map(|(hashval, positions)| {
                let freq = self.stats.kmer_frequencies.get(hashval).copied().unwrap_or(1) as f64
                    / self.stats.total_signatures as f64;
                let n_in_window =
                    positions.iter().filter(|&&p| p >= window_start && p < window_end).count();
                freq * n_in_window as f64
            })
            .sum()
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
    // chosen after the fact (the best-looking gapless run is kept), so a per-region p-value
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
/// let query = ProteinSketch::from_protein_sequence("query", "ATCGATCG", 10, 1, "hp").unwrap();
/// let target = ProteinSketch::from_protein_sequence("target", "ATCGATCG", 10, 1, "hp").unwrap();
/// let result = calculate_similarity(&query, &target);
/// ```
#[must_use]
pub fn calculate_similarity(query: &ProteinSketch, target: &ProteinSketch) -> Option<SearchResult> {
    let query_mins = query.mins_as_set();
    let target_mins = target.mins_as_set();
    let intersection: HashSet<u64> = query_mins.intersection(&target_mins).cloned().collect();
    calculate_similarity_from_precomputed(query, &query_mins, target, &target_mins, &intersection)
}

/// Find all consecutive matched regions of k-mer overlap between a query and target sequences
///
/// WHY: This is a standalone function because it doesn't require any state from ProteinSearcher.
/// It only operates on the sketches and intersection provided. This makes it easier to test and
/// more reusable. This is idiomatic Rust - functions that don't need state should be standalone.
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

    if query_target_pairs.is_empty() {
        return Vec::new();
    }

    // Sort pairs by query position, then by target position
    // WHY: This allows us to efficiently find consecutive regions in both query and target.
    query_target_pairs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

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
                        expected_shared_kmers: 0.0,
                        poisson_pvalue: 1.0,
                        enrichment: 0.0,
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

        // Validate that moltype sequences match (they should since they share the same k-mers).
        if query_moltype_seq != target_moltype_seq {
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
            expected_shared_kmers: 0.0,
            poisson_pvalue: 1.0,
            enrichment: 0.0,
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
            moltype: MolType::new("hp").unwrap(),
            length: 6,
            expected_shared_kmers: 2.0,
            poisson_pvalue: 0.05,
            enrichment: 1.5,
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
            moltype: "hp".to_string(),
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

        let row = SearchResultCsv::from_result_and_region(&result, &region);

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
        // Region-scoped stat columns: length 6, ksize 5 -> 6 - 5 + 1 = 2 shared k-mers.
        assert_eq!(row.region_n_shared_kmers, 2);
        assert_eq!(row.region_expected_shared_kmers, 2.0);
        assert_eq!(row.region_poisson_pvalue, 0.05);
        assert_eq!(row.region_enrichment, 1.5);
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
        ProteinSketch::from_protein_sequence(&name, &seq, 12, 1, "hp").unwrap()
    }

    #[fixture]
    fn ced9_sketch_k15() -> ProteinSketch {
        let (name, seq) = read_first_fasta_record(TEST_CED9_FASTA).unwrap();
        ProteinSketch::from_protein_sequence(&name, &seq, 15, 1, "hp").unwrap()
    }

    #[fixture]
    fn bcl2_sketch_k12() -> ProteinSketch {
        let (name, seq) = read_first_fasta_record(TEST_BLC2_FASTA).unwrap();
        ProteinSketch::from_protein_sequence(&name, &seq, 12, 1, "hp").unwrap()
    }

    #[fixture]
    fn bcl2_sketch_k15() -> ProteinSketch {
        let (name, seq) = read_first_fasta_record(TEST_BLC2_FASTA).unwrap();
        ProteinSketch::from_protein_sequence(&name, &seq, 15, 1, "hp").unwrap()
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
            15,    // ksize - k=15 is where BCL2/CED9 have good HP overlap
            1,     // scaled
            "hp",  // moltype
            false, // store_raw_sequences
        )?;

        target_index.process_fasta(target_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        // Create searcher
        let searcher = ProteinSearcher::new(target_index);

        // Create query index (BCL2)
        let query_index_path = temp_path.join("query_index");
        let query_index = ProteomeIndex::new(
            &query_index_path,
            15,    // ksize
            1,     // scaled
            "hp",  // moltype
            false, // store_raw_sequences
        )?;

        query_index.process_fasta(query_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        // Get query signatures
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
    /// either alone passing keeps the pair) the two p-value scopes. Every other field is left
    /// at its permissive default in each case, so the field under test is what causes the
    /// rejection, not some other, stricter field.
    #[test]
    fn test_search_filters_reject_candidates() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path();

        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 15, 1, "hp", false)?;
        target_index.process_fasta(
            TEST_CED9_FASTA,
            DEFAULT_PROGRESS_INTERVAL,
            DEFAULT_BATCH_SIZE,
        )?;
        let searcher = ProteinSearcher::new(target_index);

        let query_index_path = temp_path.join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, 15, 1, "hp", false)?;
        query_index.process_fasta(
            TEST_BLC2_FASTA,
            DEFAULT_PROGRESS_INTERVAL,
            DEFAULT_BATCH_SIZE,
        )?;
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

        // Both p-values are in [0.0, 1.0] and the checks are strict <, so capping both scopes at
        // 0.0 rejects everything. Both must be capped: either one passing keeps the pair.
        let pvalue_filtered = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.0,
                max_region_pvalue: 0.0,
                ..SearchFilters::default()
            },
        )?;
        assert!(pvalue_filtered.is_empty(), "capping both p-value scopes at 0.0 rejects all");

        Ok(())
    }

    /// The two p-value scopes combine with OR: a match is kept if either one passes. BCL2/CED9
    /// at k=15 is a weak whole-query match (p ~ 0.99) carrying one strong region (p ~ 0.0007).
    /// Capping only the query scope must not discard it, and capping only the region scope
    /// must keep it.
    #[test]
    fn test_pvalue_scopes_combine_with_or() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 15, 1, "hp", true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(target_index);

        let query_index_path = temp_dir.path().join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, 15, 1, "hp", true)?;
        query_index.process_fasta(TEST_CED9_FASTA, 0, DEFAULT_BATCH_SIZE)?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        let find_bcl2 =
            |results: &[SearchResult]| results.iter().any(|r| r.target_name.contains("BCL2_HUMAN"));

        // Whole-query scope alone rejects it: its query p-value is ~0.99, far above 0.05.
        let query_only = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.05,
                max_region_pvalue: 0.0,
                ..SearchFilters::default()
            },
        )?;
        assert!(!find_bcl2(&query_only), "query scope alone should reject the diluted BCL2 match");

        // Region scope alone keeps it, on the strength of the one significant region.
        let region_only = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.0,
                max_region_pvalue: 0.05,
                ..SearchFilters::default()
            },
        )?;
        assert!(find_bcl2(&region_only), "region scope should keep the sub-protein domain hit");

        // Both at the CLI default: OR means the region rescues it.
        let both = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.05,
                max_region_pvalue: 0.05,
                ..SearchFilters::default()
            },
        )?;
        assert!(find_bcl2(&both), "default OR semantics should surface the hit");

        Ok(())
    }

    /// `ProteinSearcher::new()` (used by most tests) keeps signatures in an in-memory DashMap,
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
        let target_index = ProteomeIndex::new(&target_index_path, 12, 1, "hp", false)?;
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

    #[test]
    fn test_find_matched_regions_multiple() -> Result<()> {
        // 14 is the minimum k-mersize that finds multiple match regions from Delilah's analyses
        let ksize = 12;
        let scaled = 1;
        let moltype = "hp";

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
        let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp", false)?;

        target_index.process_fasta(&target_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        // Create searcher
        let searcher = ProteinSearcher::new(target_index);

        // Create query FASTA
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">query\nATCGATCGATCGATCG")?;

        let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;

        query_index.process_fasta(&query_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

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
        let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp", false)?;
        target_index.process_fasta(&target_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        let searcher = ProteinSearcher::new(target_index);

        let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;
        query_index.process_fasta(&query_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

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
        let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp", false)?;
        target_index.process_fasta(&target_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        let searcher = ProteinSearcher::new(target_index);

        // Create query
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">query\nATCGATCGATCGATCG")?;

        let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;
        query_index.process_fasta(&query_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

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
        let index = ProteomeIndex::new(&temp_path, 10, 5, "hp", false)?;

        let query = ProteinSketch::new("test", 10, 5, "hp")?;
        let stats = SearchStats {
            total_signatures: 100,
            idf: HashMap::new(),
            kmer_frequencies: HashMap::new(),
        };

        let (target_list, inverted_index) = ProteinSearcher::build_search_structures(&index);
        let searcher = ProteinSearcher {
            index,
            stats,
            target_list,
            inverted_index,
            sig_cache: DashMap::new(),
            query_kmer_frequencies: None,
            total_queries: 0,
            db_n_kmers: 0,
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
        let moltype = "hp";

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
        let searcher = ProteinSearcher::new(target_index);

        // Create query index from CED9 in a temporary directory
        // WHY: new_with_auto_filename creates the database next to the input file, which causes
        // RocksDB lock conflicts when tests run in parallel. Using a temporary directory ensures
        // each test run has its own isolated database path.
        let query_index_path = temp_path.join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, ksize, scaled, moltype, true)?;

        query_index.process_fasta(TEST_CED9_FASTA, 0, DEFAULT_BATCH_SIZE)?;

        // Get query signatures
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
        assert_eq!(bcl2_result.moltype, moltype);

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
        assert!(
            bcl2_result.query_tfidf == 565.119680433367,
            "query_tfidf should be 565.119680433367, got {}",
            bcl2_result.query_tfidf
        );

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
        let moltype = "hp";

        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;

        let mut searcher = ProteinSearcher::new(target_index);

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

    /// Checks the CSV's `region_n_shared_kmers` shortcut (`region.length - ksize + 1`, see
    /// `SearchResultCsv::from_result_and_region`) against an independent count of real k-mer
    /// positions, for every named HP alphabet.
    ///
    /// That formula is only correct because scaled=1 means FracMinHash keeps every k-mer (no
    /// downsampling), so "how many k-mers are in this span" reduces to arithmetic on the span's
    /// length. Rather than trusting that reasoning, this test recomputes the count a different
    /// way: for each region, it walks `kmer_positions` (the sketch's own record of where each
    /// retained k-mer starts) and counts how many positions fall inside the region's span,
    /// then asserts that matches the formula's answer. Repeated for every named HP
    /// alphabet (including hp_thomas_dill_no_c) since each partitions residues into H/P
    /// differently, and the position bookkeeping has to hold for all of them, not just one.
    #[test]
    fn test_region_shared_kmer_count_exact_at_scaled_one_all_alphabets() {
        use crate::hp_alphabets::HpAlphabet;

        // Real sequence (already used elsewhere in this codebase for HP-alphabet regression
        // tests, see test_hp_encoding.rs) - self-hit so every alphabet reliably finds a region.
        let seq = "MKTAYIAKQRFLVSNSQLAGKRILVTQADTFMGPTLCEVFAEMG";
        let ksize = 8;

        for alpha in HpAlphabet::all_named() {
            let moltype = alpha.to_moltype();
            let sketch = ProteinSketch::from_protein_sequence("self", seq, ksize, 1, &moltype)
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
                    n_in_region as u32,
                    region.length - ksize + 1,
                    "{moltype}: region {:?} shared k-mer count mismatch",
                    region.subseq
                );
            }
        }
    }

    /// The same 5 shared k-mers that read as a weak whole-protein match (containment ~0.019
    /// against CED9's 266 k-mers) are the entire signal inside their own 19aa region. Runs a
    /// real database search (needed for region.poisson_pvalue's DB context) and checks the
    /// region-scoped Poisson test independently: recomputes lambda by hand from the searcher's
    /// own background frequencies, restricted to the region's span, and checks it against
    /// region.expected_shared_kmers/region.poisson_pvalue rather than trusting the same code
    /// path that produced them.
    #[test]
    fn test_region_poisson_pvalue_independently_recomputed() -> Result<()> {
        let ksize = 15;
        let scaled = 1;
        let moltype = "hp";

        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(target_index);

        let query_index_path = temp_dir.path().join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, ksize, scaled, moltype, true)?;
        query_index.process_fasta(TEST_CED9_FASTA, 0, DEFAULT_BATCH_SIZE)?;
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
        assert_relative_eq!(region.poisson_pvalue, pvalue_by_hand, epsilon = 1e-12);

        // The region-scoped null (over ~5 background-frequency k-mers) is a much smaller number
        // than the whole-protein null (over all 266 of CED9's k-mers), so the two p-values are
        // computed from different lambdas and shouldn't coincide.
        assert_ne!(region.poisson_pvalue, bcl2_result.query_poisson_pvalue);

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
        let moltype = "hp";

        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(target_index);

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
                solo_region.poisson_pvalue,
                paired_region.poisson_pvalue,
                epsilon = 1e-12
            );
        }

        Ok(())
    }

    /// The two p-value scopes are combined with OR, so the full truth table matters: a single
    /// scope clearing is enough, and only both failing rejects. Exercised directly here because
    /// through `search()` the degenerate caps are hard to reach.
    #[rstest]
    // query passes, region fails -> kept on the query scope
    #[case(0.001, Some(0.9), 0.05, 0.05, true)]
    // query fails, region passes -> kept on the region scope (the BCL2/CED9 shape)
    #[case(0.99, Some(0.0007), 0.05, 0.05, true)]
    // both pass
    #[case(0.001, Some(0.0007), 0.05, 0.05, true)]
    // both fail -> rejected
    #[case(0.99, Some(0.9), 0.05, 0.05, false)]
    // no regions at all: the region disjunct is vacuously false, query alone decides
    #[case(0.001, None, 0.05, 0.05, true)]
    #[case(0.99, None, 0.05, 0.05, false)]
    // a cap of 0.0 can never be cleared (the check is a strict <), which is how the deprecated
    // --max-pvalue alias reduces to whole-query filtering
    #[case(0.001, Some(0.0), 0.05, 0.0, true)]
    #[case(0.99, Some(0.0007), 0.05, 0.0, false)]
    // infinite caps accept anything, including the p = 1.0 of a no-DB-context result
    #[case(1.0, Some(1.0), f64::INFINITY, f64::INFINITY, true)]
    fn test_pvalues_pass_truth_table(
        #[case] query_pvalue: f64,
        #[case] best_region_pvalue: Option<f64>,
        #[case] max_query_pvalue: f64,
        #[case] max_region_pvalue: f64,
        #[case] expected: bool,
    ) {
        let filters =
            SearchFilters { max_query_pvalue, max_region_pvalue, ..SearchFilters::default() };
        assert_eq!(filters.pvalues_pass(query_pvalue, best_region_pvalue), expected);
    }

    /// The Poisson survival function is the shared engine behind both p-value scopes. The
    /// degenerate inputs return 1.0 (no evidence) instead of erroring or producing NaN, which
    /// a real search cannot reach but a caller can.
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
        let index = ProteomeIndex::new(&index_path, ksize, 1, "hp", true)?;
        index.process_fasta(TEST_CED9_FASTA, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(index);

        let (name, sequence) = read_first_fasta_record(TEST_CED9_FASTA)?;
        let sketch = ProteinSketch::from_protein_sequence(&name, &sequence, ksize, 1, "hp")?;

        // A window holding one k-mer: [0, 0 + 1) after the ksize adjustment.
        let single = searcher.region_expectation(&sketch, 0, ksize, ksize as usize);
        assert_relative_eq!(single, 1.0, epsilon = 1e-12);

        // Widening the region by one residue admits one more k-mer start.
        let double = searcher.region_expectation(&sketch, 0, ksize + 1, ksize as usize);
        assert_relative_eq!(double, 2.0, epsilon = 1e-12);

        // A span shorter than k contains no whole k-mer, so there is nothing to expect.
        let too_short = searcher.region_expectation(&sketch, 0, ksize - 1, ksize as usize);
        assert_eq!(too_short, 0.0);

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
        let index = ProteomeIndex::new(&index_path, ksize, 1, "hp", true)?;
        index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        assert_eq!(index.signature_count(), 25);

        let searcher = ProteinSearcher::new(index);
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
    /// direction. BCL2A1 vs ASPP2/TP53BP2 at k=9, in the same 25-sequence fixture database, is
    /// an overwhelming whole-protein match (115 shared k-mers scattered across 333 short
    /// regions) with no single region concentrated enough to pass on its own. Its strongest
    /// region only reaches p=0.0956, above the 0.05 default cap. This shows the OR only needs
    /// one scope to hold, in either direction.
    #[test]
    fn test_query_scope_alone_keeps_a_diffuse_match_with_no_standout_region() -> Result<()> {
        let ksize = 9;
        let scaled = 1;
        let moltype = "hp";

        let temp_dir = TempDir::new()?;
        let target_index_path = temp_dir.path().join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, ksize, scaled, moltype, true)?;
        target_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let searcher = ProteinSearcher::new(target_index);

        let query_index_path = temp_dir.path().join("query_index");
        let query_index = ProteomeIndex::new(&query_index_path, ksize, scaled, moltype, true)?;
        query_index.process_fasta(TEST_FASTA_GZ, 0, DEFAULT_BATCH_SIZE)?;
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        fn find_hit(results: &[SearchResult]) -> Option<&SearchResult> {
            results
                .iter()
                .find(|r| r.query_name.contains("B2LA1") && r.target_name.contains("ASPP2"))
        }

        // Unfiltered, to inspect the pair's raw numbers.
        let all_results = searcher.search(&query_signatures, &SearchFilters::default())?;
        let hit = find_hit(&all_results).expect("BCL2A1 vs ASPP2/TP53BP2 should be found");

        assert_eq!(hit.n_intersecting_hashes, 115);
        assert_relative_eq!(hit.query_poisson_pvalue, 2.356_930_483e-7, epsilon = 1e-15);
        assert!(hit.query_poisson_pvalue < 0.05, "whole-query scope should clearly pass");

        let best_region_pvalue = hit
            .matched_regions
            .iter()
            .map(|region| region.poisson_pvalue)
            .fold(f64::INFINITY, f64::min);
        assert_relative_eq!(best_region_pvalue, 0.095_589_196_102_397_8, epsilon = 1e-12);
        assert!(best_region_pvalue >= 0.05, "no single region should clear the default cap");

        // Region scope alone: nothing to rescue it, since no region is significant on its own.
        let region_only = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.0,
                max_region_pvalue: 0.05,
                ..SearchFilters::default()
            },
        )?;
        assert!(
            find_hit(&region_only).is_none(),
            "region scope alone should reject a match with no standout region"
        );

        // Query scope alone: the diffuse whole-protein signal is sufficient by itself.
        let query_only = searcher.search(
            &query_signatures,
            &SearchFilters {
                max_query_pvalue: 0.05,
                max_region_pvalue: 0.0,
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
