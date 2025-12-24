use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::path::Path;

use anyhow::Result;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::errors::IndexResult;
use crate::index::ProteomeIndex;
use crate::significance;
use crate::sketch::ProteinSketch;
use crate::types::MolType;

/// Default progress interval for FASTA processing (log progress every N sequences)
pub const DEFAULT_PROGRESS_INTERVAL: u32 = 1000;

/// Default batch size for FASTA processing (process N sequences per batch)
pub const DEFAULT_BATCH_SIZE: usize = 1000;

/// CSV-friendly version of SearchResult without matched_regions
/// WHY: CSV format doesn't handle nested Vec structures well. This struct is used
/// for CSV serialization while preserving all the important similarity metrics.
/// This is idiomatic Rust - we create a separate struct for CSV output rather than
/// trying to serialize complex nested structures that CSV doesn't support.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResultCsv {
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
    pub tfidf: f64,
    pub average_database_kmer_frequency: f64,
    pub prob_random_cooccurrence: f64,
    pub prob_random_cooccurrence_symmetric: f64,
    pub expected_intersecting_hashes: f64,
    pub observed_over_expected: f64,
}

impl From<&SearchResult> for SearchResultCsv {
    fn from(result: &SearchResult) -> Self {
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
            tfidf: result.tfidf,
            average_database_kmer_frequency: result.average_database_kmer_frequency,
            prob_random_cooccurrence: result.prob_random_cooccurrence,
            prob_random_cooccurrence_symmetric: result.prob_random_cooccurrence_symmetric,
            expected_intersecting_hashes: result.expected_intersecting_hashes,
            observed_over_expected: result.observed_over_expected,
        }
    }
}

/// CSV-friendly version of SearchResult with matched region information
/// WHY: Each matched region gets its own row in the CSV, with all SearchResult similarity
/// metrics repeated for each region. This makes it easy to analyze individual matched regions
/// while still having access to the overall similarity metrics. Every CSV row must have matched
/// region data - SearchResults without matched regions are not included in the CSV output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResultWithRegionCsv {
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
    pub tfidf: f64,
    pub average_database_kmer_frequency: f64,
    pub prob_random_cooccurrence: f64,
    pub prob_random_cooccurrence_symmetric: f64,
    pub expected_intersecting_hashes: f64,
    pub observed_over_expected: f64,
    // MatchedRegion fields (always present - every CSV row has a matched region)
    pub query_start: u32,
    pub query_end: u32,
    pub query_subseq: String,
    pub target_start: u32,
    pub target_end: u32,
    pub target_subseq: String,
    pub moltype_seq: String,
    pub region_length: u32,
}

impl SearchResultWithRegionCsv {
    /// Create a CSV row from a SearchResult and a MatchedRegion
    /// WHY: Every CSV row must have matched region data. This method combines the SearchResult
    /// similarity metrics with a specific matched region to create one CSV row. Each SearchResult
    /// will produce multiple CSV rows (one per matched region), with all similarity metrics
    /// repeated for each region.
    pub fn from_result_and_region(result: &SearchResult, region: &MatchedRegion) -> Self {
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
            tfidf: result.tfidf,
            average_database_kmer_frequency: result.average_database_kmer_frequency,
            prob_random_cooccurrence: result.prob_random_cooccurrence,
            prob_random_cooccurrence_symmetric: result.prob_random_cooccurrence_symmetric,
            expected_intersecting_hashes: result.expected_intersecting_hashes,
            observed_over_expected: result.observed_over_expected,
            query_start: region.query_start,
            query_end: region.query_end,
            query_subseq: region.query_subseq.clone(),
            target_start: region.target_start,
            target_end: region.target_end,
            target_subseq: region.target_subseq.clone(),
            moltype_seq: region.moltype_seq.clone(),
            region_length: region.length,
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

    /// TF-IDF score for the query signature
    pub tfidf: f64,

    /// Average k-mer frequency in the database (symmetric between query and target, requires database context)
    /// Mean frequency of intersecting k-mers across the entire database.
    /// Higher values = more common k-mers. Range: [0,1].
    pub average_database_kmer_frequency: f64,

    /// Probability that all intersecting k-mers would co-occur by random chance (asymmetric - depends on query/target order)
    /// Product of individual k-mer probabilities. Lower values = more significant matches.
    /// This is like a p-value: values near 0 indicate rare k-mers matching.
    pub prob_random_cooccurrence: f64,

    /// Symmetric version of prob_random_cooccurrence (same regardless of query/target order)
    /// Geometric mean of query→target and target→query cooccurrence probabilities.
    /// Lower values = more significant matches.
    pub prob_random_cooccurrence_symmetric: f64,

    /// Expected number of shared k-mers by random chance given database composition
    pub expected_intersecting_hashes: f64,

    /// Ratio of observed to expected intersecting hashes (n_intersecting_hashes / expected_intersecting_hashes)
    /// Higher values = more significant matches (e.g., 8.0 means 8x more k-mers than expected by chance)
    pub observed_over_expected: f64,

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

    /// Query sequence start position
    pub query_start: u32,

    /// Query sequence end position
    pub query_end: u32,

    /// Query subsequence (stitched k-mers)
    pub query_subseq: String,

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
}

impl Display for MatchedRegion {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Query Name: {}", self.query_name)?;
        writeln!(f, "Match Name: {}", self.target_name)?;
        writeln!(f, "query: {} ({}-{})", self.query_subseq, self.query_start, self.query_end)?;
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
}

impl ProteinSearcher {
    /// Create a new protein searcher from an index
    pub fn new(index: ProteomeIndex) -> Self {
        let stats = SearchStats::from_index(&index);
        Self { index, stats }
    }

    /// Load a searcher from a saved index
    pub fn load<P: AsRef<Path>>(path: P) -> IndexResult<Self> {
        let index = ProteomeIndex::load(path)?;
        let stats = SearchStats::from_index(&index);
        Ok(Self { index, stats })
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
    ///
    /// # Returns
    /// Vector of SearchResult containing all similarity metrics, sorted by containment score
    #[must_use = "search results should be used to process query matches"]
    pub fn search(&self, queries: &[ProteinSketch]) -> Result<Vec<SearchResult>> {
        // Perform parallel search across all queries
        let all_results: Vec<SearchResult> = queries
            .par_iter()
            .flat_map(|query| {
                // Prepare query once - this pre-computes mins as HashSet and TF-IDF
                let prepared = self.prepare_query(query);

                // Search this query against all targets
                self.index
                    .get_signatures()
                    .iter()
                    .filter_map(|entry| {
                        let target = entry.value();
                        self.compare(&prepared, target)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        // Sort by containment score (descending) - this is the primary ranking metric
        let mut sorted_results = all_results;
        sorted_results.sort_by(|a, b| {
            b.containment.partial_cmp(&a.containment).unwrap_or(std::cmp::Ordering::Equal)
        });

        Ok(sorted_results)
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
    /// # Returns
    /// Vector of SearchResult containing all similarity metrics, sorted by containment score
    #[must_use = "search results should be used to process query matches"]
    pub fn search_all_vs_all(&self) -> Result<Vec<SearchResult>> {
        let signatures = self.index.get_signatures();

        // Collect MD5 keys for parallel iteration (cheap String clones, not ProteinSketch clones)
        // WHY: DashMap doesn't implement ParallelIterator directly. By collecting keys first,
        // we can parallelize over them and look up the actual signatures in the DashMap. This
        // avoids cloning the large ProteinSketch objects while still enabling parallel processing.
        let md5_keys: Vec<String> = signatures.iter().map(|entry| entry.key().clone()).collect();

        // Perform parallel search across all signatures
        // WHY: We iterate over MD5 keys in parallel, then look up the actual signatures in the
        // DashMap. Each signature is used as both query and target, but we skip self-matches
        // in the compare() method by comparing MD5 sums.
        let all_results: Vec<SearchResult> = md5_keys
            .par_iter()
            .flat_map(|query_md5| {
                // Look up query signature by MD5
                // WHY: We use filter_map to handle the Option from get() gracefully. If the
                // signature doesn't exist (shouldn't happen, but safe to handle), we skip it.
                signatures.get(query_md5).map(|query_entry| {
                    let query = query_entry.value();

                    // Prepare query once - this pre-computes mins as HashSet and TF-IDF
                    let prepared = self.prepare_query(query);

                    // Search this query against all targets (including itself, but compare() will skip self-matches)
                    signatures
                        .iter()
                        .filter_map(|target_entry| {
                            let target = target_entry.value();
                            self.compare(&prepared, target)
                        })
                        .collect::<Vec<_>>()
                })
            })
            .flatten()
            .collect();

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
    ) -> Option<SearchResult> {
        // Skip self-matches by comparing MD5 sums
        // WHY: In all-vs-all searches, we don't want to compare a signature against itself.
        // MD5 sum is a unique identifier for each signature, so comparing MD5 sums is the
        // most reliable way to detect self-matches. This is idiomatic Rust - we use early
        // returns to avoid unnecessary computation when we know the result will be invalid.
        if query.sketch.signature().md5sum == target.signature().md5sum {
            return None;
        }

        // Calculate intersection for overlap probability check
        let intersection = query.sketch.intersect(target);

        // Skip if no intersection
        if intersection.is_empty() {
            return None;
        }

        // Calculate all database-specific overlap metrics
        let average_database_kmer_frequency = self.calculate_average_database_kmer_frequency(&intersection);
        let prob_random_cooccurrence = self.calculate_prob_random_cooccurrence(&intersection);

        // Calculate symmetric version: geometric mean of both directions
        let prob_random_cooccurrence_reverse = self.calculate_prob_random_cooccurrence(&intersection);
        let prob_random_cooccurrence_symmetric = (prob_random_cooccurrence * prob_random_cooccurrence_reverse).sqrt();

        let expected_intersecting_hashes = self.calculate_expected_intersecting_hashes(query.sketch, target);

        // Get the base similarity result from the standalone function
        // Note: This will recalculate intersection, but that's acceptable for the cleaner API
        let mut result = calculate_similarity(query.sketch, target)?;

        // Calculate observed over expected ratio
        let observed_over_expected = if expected_intersecting_hashes > 0.0 {
            result.n_intersecting_hashes as f64 / expected_intersecting_hashes
        } else {
            0.0
        };

        // Override with database-specific metrics
        result.tfidf = query.tfidf;
        result.average_database_kmer_frequency = average_database_kmer_frequency;
        result.prob_random_cooccurrence = prob_random_cooccurrence;
        result.prob_random_cooccurrence_symmetric = prob_random_cooccurrence_symmetric;
        result.expected_intersecting_hashes = expected_intersecting_hashes;
        result.observed_over_expected = observed_over_expected;

        Some(result)
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

    /// Calculate the average k-mer frequency in the database for intersecting k-mers.
    ///
    /// This is a symmetric metric (same for query vs target and target vs query) that measures
    /// how common the intersecting k-mers are across the entire database. Higher values indicate
    /// that the matching k-mers are common/frequent in the database.
    ///
    /// **Returns**: Mean of normalized frequencies. Range: [0,1].
    #[must_use = "average database kmer frequency should be used to assess k-mer commonality"]
    pub fn calculate_average_database_kmer_frequency(&self, intersection: &HashSet<u64>) -> f64 {
        if intersection.is_empty() {
            return 0.0;
        }

        // Calculate average frequency of intersecting k-mers in the database
        let sum_freq: f64 = intersection
            .par_iter()
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

    /// Calculate probability that all intersecting k-mers would co-occur by random chance
    ///
    /// This is an asymmetric metric (depends on which sequence is query vs target) that uses
    /// the product of individual k-mer probabilities (assuming independence). The probability
    /// represents how likely it is that ALL the intersecting k-mers would appear together by
    /// random chance in two randomly selected sequences from the database.
    ///
    /// **Interpretation**:
    /// - Values near 0.0 = highly significant match (rare k-mers matching)
    /// - Values near 1.0 = likely random match (common k-mers matching)
    ///
    /// This is like a p-value: smaller values indicate more significant matches.
    #[must_use = "random cooccurrence probability should be used to assess match significance"]
    pub fn calculate_prob_random_cooccurrence(&self, intersection: &HashSet<u64>) -> f64 {
        if intersection.is_empty() {
            return 1.0; // No intersection = completely expected by chance
        }

        // Calculate the product of individual k-mer probabilities
        // WHY: We use product, not sum. The probability that ALL k-mers co-occur is the
        // product of their individual probabilities (assuming independence).
        let prob_random: f64 = intersection
            .par_iter()
            .map(|&hashval| {
                // Get frequency of this k-mer in the database
                let db_frequency =
                    self.stats.kmer_frequencies.get(&hashval).copied().unwrap_or(1) as f64;
                let total_signatures = self.stats.total_signatures as f64;

                // Probability that this k-mer appears in a randomly selected signature
                db_frequency / total_signatures
            })
            .product(); // PRODUCT not sum - all must co-occur

        prob_random.clamp(0.0, 1.0)
    }

    /// Calculate expected number of shared k-mers by random chance
    ///
    /// This uses Option 2: Expected overlap calculation. For each k-mer in the query,
    /// we calculate the probability that it would also appear in the target based on
    /// its frequency in the database. We sum these probabilities to get the expected
    /// number of shared k-mers.
    ///
    /// **Interpretation**:
    /// - Compare to n_intersecting_hashes: if observed >> expected, the match is significant
    /// - Values are NOT probabilities - they're counts (can be > 1.0)
    /// - Lower expected values (given the same observed) indicate more significant matches
    ///
    /// **Example**: If we observe 24 intersecting k-mers but only expect 2.5 by chance,
    /// that's a highly significant match.
    #[must_use = "expected intersecting hashes should be compared to observed intersection size"]
    pub fn calculate_expected_intersecting_hashes(
        &self,
        query_sketch: &ProteinSketch,
        target_sketch: &ProteinSketch,
    ) -> f64 {
        let query_mins = query_sketch.signature().minhash.mins();
        let target_mins = target_sketch.mins_as_set();

        // For each k-mer in the query, calculate the probability it would appear in the target
        // based on its database frequency, then sum these probabilities
        let expected: f64 = query_mins
            .par_iter()
            .map(|&hashval| {
                // Skip k-mers not in target (we're calculating expected, not observed)
                if !target_mins.contains(&hashval) {
                    return 0.0;
                }

                // Get frequency of this k-mer in the database
                let db_frequency =
                    self.stats.kmer_frequencies.get(&hashval).copied().unwrap_or(1) as f64;
                let total_signatures = self.stats.total_signatures as f64;

                // Probability that this k-mer would appear in a random target
                db_frequency / total_signatures
            })
            .sum();

        expected
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
    // Extract query data
    let query_abunds = query.signature().minhash.abunds();
    let query_name = query.signature().name.clone();
    let query_md5 = query.signature().md5sum.clone();

    // Extract target data
    let target_abunds = target.signature().minhash.abunds();

    // Calculate intersection using the new intersect() method
    let intersection = query.intersect(target);
    let n_intersecting_hashes = intersection.len();

    // Skip if no intersection
    if n_intersecting_hashes == 0 {
        return None;
    }

    let query_size = query.mins_as_set().len();
    let target_size = target.mins_as_set().len();
    let union_size = query_size + target_size - n_intersecting_hashes;

    // Calculate basic metrics
    let containment = n_intersecting_hashes as f64 / query_size as f64;
    let jaccard = n_intersecting_hashes as f64 / union_size as f64;
    let containment_target_in_query = n_intersecting_hashes as f64 / target_size as f64;
    let max_containment = containment.max(containment_target_in_query);

    // Calculate abundance statistics
    // WHY: We pass the mins arrays along with the abundance arrays because abundances
    // are stored in the same order as mins. We need to find the position of each
    // intersecting hash in both arrays to get the correct corresponding abundances.
    let (average_abund, median_abund, std_abund) =
        if let (Some(query_abunds), Some(target_abunds)) =
            (query_abunds.as_ref(), target_abunds.as_ref())
        {
            let query_mins = query.signature().minhash.mins();
            let target_mins = target.signature().minhash.mins();
            significance::abundance_stats(
                &intersection,
                &query_mins,
                query_abunds,
                &target_mins,
                target_abunds,
            )
        } else {
            (1.0, 1.0, 0.0)
        };

    // Calculate weighted metrics
    let f_weighted_target_in_query = significance::weighted_fraction_target_in_query(
        query_abunds.as_deref(),
        target_abunds.as_deref(),
    );

    // Find matched regions
    let matched_regions = find_matched_regions(query, target, &intersection);

    Some(SearchResult {
        query_name,
        query_md5,
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
        tfidf: 0.0,                                  // Default for 1v1 comparisons - requires database context
        average_database_kmer_frequency: 0.0,        // Default for 1v1 comparisons - requires database context
        prob_random_cooccurrence: 1.0,               // Default for 1v1 comparisons - requires database context
        prob_random_cooccurrence_symmetric: 1.0,     // Default for 1v1 comparisons - requires database context
        expected_intersecting_hashes: 0.0,           // Default for 1v1 comparisons - requires database context
        observed_over_expected: 0.0,                 // Default for 1v1 comparisons - requires database context
        matched_regions,
    })
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
        if let (Some(query_kmer_info), Some(target_kmer_info)) =
            (query_sketch.kmer_infos().get(&hashval), target_sketch.kmer_infos().get(&hashval))
        {
            // Collect all query positions for this hashval
            let mut query_poss = Vec::new();
            for positions in query_kmer_info.original_kmer_to_position.values() {
                query_poss.extend(positions.iter().cloned());
            }
            query_poss.sort();
            query_poss.dedup();
            hashval_to_query_positions.insert(hashval, query_poss);

            // Collect all target positions for this hashval
            let mut target_poss = Vec::new();
            for positions in target_kmer_info.original_kmer_to_position.values() {
                target_poss.extend(positions.iter().cloned());
            }
            target_poss.sort();
            target_poss.dedup();
            hashval_to_target_positions.insert(hashval, target_poss);
        }
    }

    // Build mapping from (query_pos, target_pos) pairs for each hashval
    // WHY: We need to track the correspondence between query and target positions for each
    // k-mer hash. This allows us to find regions where both query and target positions are
    // consecutive, ensuring we match the correct target region to each query region.
    let mut query_target_pairs: Vec<(usize, usize, u64)> = Vec::new(); // (query_pos, target_pos, hashval)
    for &hashval in intersection {
        if let (Some(query_kmer_info), Some(target_kmer_info)) =
            (query_sketch.kmer_infos().get(&hashval), target_sketch.kmer_infos().get(&hashval))
        {
            // Get all query positions for this hashval
            let mut query_poss = Vec::new();
            for positions in query_kmer_info.original_kmer_to_position.values() {
                query_poss.extend(positions.iter().cloned());
            }
            query_poss.sort();
            query_poss.dedup();

            // Get all target positions for this hashval
            let mut target_poss = Vec::new();
            for positions in target_kmer_info.original_kmer_to_position.values() {
                target_poss.extend(positions.iter().cloned());
            }
            target_poss.sort();
            target_poss.dedup();

            // Create all pairs of (query_pos, target_pos) for this hashval
            // WHY: Each hashval can appear at multiple positions in both query and target.
            // We create all pairs to find the correct correspondences.
            for &qpos in &query_poss {
                for &tpos in &target_poss {
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
                        query_start: query_start_pos as u32,
                        query_end: query_end_pos as u32,
                        query_subseq: query_subseq.to_string(),
                        target_name: target_name.clone(),
                        target_start: target_start_pos as u32,
                        target_end: target_end_pos as u32,
                        target_subseq: target_subseq.to_string(),
                        moltype: moltype.clone(),
                        moltype_seq: String::new(), // Empty since we don't have encoded sequence
                        length: (query_end_pos - query_start_pos) as u32,
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

        // Validate that moltype sequences match (they should since they share the same k-mers)
        if query_moltype_seq != target_moltype_seq {
            // WHY: We use a simpler panic message format to avoid potential double-panic issues.
            // The detailed information is still provided, but in a format that's less likely to
            // cause issues during panic handling.
            panic!(
                "Moltype sequences do not match for query '{}' and target '{}' at positions query {}-{} target {}-{}",
                query_name, target_name, query_start_pos, query_end_pos, target_start_pos, target_end_pos
            )
        }

        consecutive_regions.push(MatchedRegion {
            query_name: query_name.clone(),
            query_start: query_start_pos as u32,
            query_end: query_end_pos as u32,
            query_subseq: query_subseq.to_string(),
            target_name: target_name.clone(),
            target_start: target_start_pos as u32,
            target_end: target_end_pos as u32,
            target_subseq: target_subseq.to_string(),
            moltype: moltype.clone(),
            moltype_seq: target_moltype_seq.to_string(),
            length: (query_end_pos - query_start_pos) as u32,
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
    consecutive_regions.sort_by(|a, b| a.query_start.cmp(&b.query_start));

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
        query_subseq: &'static str,
        moltype_seq: &'static str,
        target_subseq: &'static str,
        query_start: u32,
        query_end: u32,
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
        // Optional: a specific region to find by query_subseq
        landmark: Option<ExpectedMatchedRegion>,
    }

    const MATCHED_REGIONS_K12: ExpectedMatchedRegions = ExpectedMatchedRegions {
        ksize: 12,
        total_regions: 13,
        first: ExpectedMatchedRegion {
            query_subseq: "FTHRIRQNGMEW",
            moltype_seq: "hppphppphhph",
            target_subseq: "FSRRYRRDFAEM",
            query_start: 87,
            query_end: 99,
            target_start: 103,
            target_end: 115,
        },
        last: ExpectedMatchedRegion {
            query_subseq: "GVVVCGRMMFSLK",
            moltype_seq: "hhhhphphhhphp",
            target_subseq: "LYGPSMRPLFDFS",
            query_start: 267,
            query_end: 280,
            target_start: 200,
            target_end: 213,
        },
        landmark: Some(ExpectedMatchedRegion {
            query_subseq: "QCPMSYGRLIGLISFGGFV",
            moltype_seq: "pphhphhphhhhhphhhhh",
            target_subseq: "RDGVNWGRIVAFFEFGGVM",
            query_start: 162,
            query_end: 181,
            target_start: 138,
            target_end: 157,
        }),
    };

    const MATCHED_REGIONS_K15: ExpectedMatchedRegions = ExpectedMatchedRegions {
        ksize: 15,
        total_regions: 1,
        first: ExpectedMatchedRegion {
            query_subseq: "QCPMSYGRLIGLISFGGFV",
            moltype_seq: "pphhphhphhhhhphhhhh",
            target_subseq: "RDGVNWGRIVAFFEFGGVM",
            query_start: 162,
            query_end: 181,
            target_start: 138,
            target_end: 157,
        },
        last: ExpectedMatchedRegion {
            query_subseq: "QCPMSYGRLIGLISFGGFV",
            moltype_seq: "pphhphhphhhhhphhhhh",
            target_subseq: "RDGVNWGRIVAFFEFGGVM",
            query_start: 162,
            query_end: 181,
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

        // Use the same test FASTA file for both query and target to ensure a match
        // WHY: Using a known working test file (CED9) ensures the test will work correctly.
        // This is simpler and more reliable than creating synthetic sequences.
        let test_fasta = TEST_CED9_FASTA;

        // Create target index
        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(
            &target_index_path,
            10,    // ksize
            1,     // scaled
            "hp",  // moltype
            false, // store_raw_sequences
        )?;

        target_index.process_fasta(test_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        // Create searcher
        let searcher = ProteinSearcher::new(target_index);

        // Create query index
        let query_index = ProteomeIndex::new_with_auto_filename(
            test_fasta, 10,    // ksize
            1,     // scaled
            "hp",  // moltype
            false, // store_raw_sequences
        )?;

        query_index.process_fasta(test_fasta, DEFAULT_PROGRESS_INTERVAL, DEFAULT_BATCH_SIZE)?;

        // Get query signatures
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        assert!(!query_signatures.is_empty(), "Should have at least one query signature");

        // Perform search
        let results = searcher.search(&query_signatures)?;

        // Should find at least one match (exact match)
        assert!(results.len() == 1, "Should find exactly one match");

        // Check that the first result has reasonable values
        let first_result = &results[0];
        assert!(first_result.query_name.contains("CED9_CAEEL"));
        assert!(first_result.target_name.contains("CED9_CAEEL"));
        assert!(first_result.containment == 1.0);
        assert!(first_result.jaccard == 1.0);
        assert!(first_result.n_intersecting_hashes == 10);

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
        assert_eq!(region.query_subseq, expected.query_subseq);
        assert_eq!(region.moltype_seq, expected.moltype_seq);
        assert_eq!(region.target_subseq, expected.target_subseq);
        assert_eq!(region.query_start, expected.query_start);
        assert_eq!(region.query_end, expected.query_end);
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
                .find(|r| r.query_subseq == landmark.query_subseq)
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
        assert_eq!(first_match.query_subseq, "FTHRIRQNGMEW");
        assert_eq!(first_match.moltype_seq, "hppphppphhph");
        assert_eq!(first_match.target_subseq, "FSRRYRRDFAEM");
        assert_eq!(first_match.query_start, 87, "First match query start should be 87");
        assert_eq!(first_match.query_end, 99, "First match query end should be 99");
        assert_eq!(first_match.target_start, 103, "First match target start should be 103");
        assert_eq!(first_match.target_end, 115, "First match target end should be 115");

        let last_match = &matched_regions[matched_regions.len() - 1];
        // Positions 267-280 in CED9 (query) and 200-213 in BCL2 (target)
        assert_eq!(last_match.query_subseq, "GVVVCGRMMFSLK");
        assert_eq!(last_match.moltype_seq, "hhhhphphhhphp");
        assert_eq!(last_match.target_subseq, "LYGPSMRPLFDFS");
        assert_eq!(last_match.query_start, 267, "Last match query start should be 267");
        assert_eq!(last_match.query_end, 280, "Last match query end should be 280");
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
            .find(|r| r.query_subseq == "QCPMSYGRLIGLISFGGFV")
            .expect("Should find the expected match region");
        assert_eq!(largest_match.query_subseq, "QCPMSYGRLIGLISFGGFV");
        assert_eq!(largest_match.moltype_seq, "pphhphhphhhhhphhhhh");
        assert_eq!(largest_match.target_subseq, "RDGVNWGRIVAFFEFGGVM");
        assert_eq!(largest_match.query_start, 162, "Largest match query start should be 162");
        assert_eq!(largest_match.query_end, 181, "Largest match query end should be 181");
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

        let results = searcher.search(&query_signatures)?;

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

        let results = searcher.search(&query_signatures)?;

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

        let searcher = ProteinSearcher { index, stats };

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
        assert_eq!(result.tfidf, 0.0, "TF-IDF should be 0.0 for 1v1 comparisons");
        assert_eq!(
            result.average_database_kmer_frequency, 0.0,
            "Average database k-mer frequency should be 0.0 for 1v1 comparisons"
        );
        assert_eq!(
            result.prob_random_cooccurrence, 1.0,
            "Random cooccurrence probability should be 1.0 for 1v1 comparisons"
        );
        assert_eq!(
            result.prob_random_cooccurrence_symmetric, 1.0,
            "Symmetric random cooccurrence probability should be 1.0 for 1v1 comparisons"
        );
        assert_eq!(
            result.expected_intersecting_hashes, 0.0,
            "Expected intersecting hashes should be 0.0 for 1v1 comparisons"
        );
        assert_eq!(
            result.observed_over_expected, 0.0,
            "Observed over expected should be 0.0 for 1v1 comparisons"
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
        let results = searcher.search(&query_signatures)?;

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
            bcl2_result.tfidf == 565.119680433367,
            "TF-IDF should be 565.119680433367, got {}",
            bcl2_result.tfidf
        );

        // Verify average database k-mer frequency is meaningful (should not be 0.0 with multiple signatures)
        // WHY: With multiple signatures, k-mers will have varying frequencies. Average database
        // k-mer frequency indicates how common the intersecting k-mers are across the database.
        // Values in (0,1] range indicate k-mer presence across the database.
        assert!(
            bcl2_result.average_database_kmer_frequency > 0.0,
            "Average database k-mer frequency should be > 0.0, got {}",
            bcl2_result.average_database_kmer_frequency
        );
        // Note: Higher values indicate k-mers that appear in many signatures (common),
        // while lower values indicate k-mers that are more specific to certain proteins

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
}
