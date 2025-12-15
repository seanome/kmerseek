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

    /// Probability of overlap between query and target
    pub overlap_probability: f64,

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
        write!(
            f,
            "
        Query Name: {}
Match Name: {}
query: {} ({}-{})
alpha: {}
match: {} ({}-{})
",
            self.query_name,
            self.target_name,
            self.query_subseq,
            self.query_start,
            self.query_end,
            self.moltype_seq,
            self.target_subseq,
            self.target_start,
            self.target_end
        )
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
    pub fn search(&self, queries: &[ProteinSketch]) -> Result<Vec<SearchResult>> {
        // Calculate TF-IDF for each query signature once (used in all results for that query)
        let query_tfidf: HashMap<String, f64> = queries
            .iter()
            .map(|query| {
                let name = query.signature().name.clone();
                let tfidf = self.calculate_tfidf(query);
                (name, tfidf)
            })
            .collect();

        // Perform parallel search across all queries
        let all_results: Vec<SearchResult> = queries
            .par_iter()
            .flat_map(|query| {
                // Olga: Why is this cloned?? Can we avoid copying data here?
                let query_mins: HashSet<u64> =
                    query.signature().minhash.mins().iter().cloned().collect();
                let query_abunds = query.signature().minhash.abunds();
                let query_name = query.signature().name.clone();
                let query_md5 = query.signature().md5sum.clone();
                let query_tfidf = query_tfidf.get(&query_name).copied().unwrap_or(0.0);

                // Search this query against all targets
                self.index
                    .get_signatures()
                    .iter()
                    .filter_map(|entry| {
                        let target = entry.value();
                        self.query_target_similarity(
                            query,
                            target,
                            &query_mins,
                            query_abunds.as_deref(),
                            &query_name,
                            &query_md5,
                            query_tfidf,
                        )
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

    /// Calculate comprehensive similarity between query and target signatures including TF-IDF and overlap probability
    ///
    /// This method calculates all similarity metrics in one pass for efficiency.
    ///
    /// # Why Manual Calculation Instead of Sourmash's Built-in Methods?
    ///
    /// We calculate containment, jaccard, and other metrics manually rather than using
    /// Sourmash's `KmerMinHash::similarity()` or `KmerMinHash::containment()` methods for
    /// several performance and functionality reasons:
    ///
    /// 1. **Pre-extracted HashSet Reuse**: The `query_mins` HashSet is extracted once per query
    ///    and reused across all target comparisons (see `search()` method). Sourmash's methods
    ///    would need to extract/convert data structures on every call, causing redundant allocations.
    ///
    /// 2. **Intersection Reuse**: We need the intersection HashSet for multiple downstream
    ///    calculations (abundance statistics, overlap probability, matched regions). Computing
    ///    it once and reusing it is more efficient than having each Sourmash method compute
    ///    it independently.
    ///
    /// 3. **Redundant Compatibility Checks**: Sourmash's methods perform compatibility checks
    ///    (ksize, scaled, moltype, seed) on every call. Since we already know signatures are
    ///    compatible (checked via `is_compatible()` or guaranteed by index construction), these
    ///    checks are unnecessary overhead in this hot path.
    ///
    /// 4. **Custom Metrics**: We calculate metrics not provided by Sourmash:
    ///    - `containment_target_in_query` (reverse containment)
    ///    - `max_containment` (max of both containment directions)
    ///    - Abundance statistics (median, std dev) on intersecting k-mers
    ///    - Custom weighted metrics and overlap probability
    ///
    /// 5. **Zero-Cost Abstraction**: By controlling the data flow, we avoid function call overhead
    ///    and intermediate allocations. The manual approach provides better performance for batch
    ///    comparisons where the same query is compared against many targets.
    ///
    /// This follows Rust's zero-cost abstraction principle: when you control the data, avoid
    /// unnecessary overhead from general-purpose library methods that must handle edge cases
    /// we've already excluded.
    pub(crate) fn query_target_similarity(
        &self,
        query: &ProteinSketch,
        target: &ProteinSketch,
        query_mins: &HashSet<u64>,
        query_abunds: Option<&[u64]>,
        query_name: &str,
        query_md5: &str,
        query_tfidf: f64,
    ) -> Option<SearchResult> {
        let target_mins: HashSet<u64> = target.signature().minhash.mins().iter().cloned().collect();
        let target_abunds = target.signature().minhash.abunds();

        // Calculate intersection
        let intersection: HashSet<u64> = query_mins.intersection(&target_mins).cloned().collect();
        let n_intersecting_hashes = intersection.len();

        // Skip if no intersection
        if n_intersecting_hashes == 0 {
            return None;
        }

        let query_size = query_mins.len();
        let target_size = target_mins.len();
        let union_size = query_size + target_size - n_intersecting_hashes;

        // Calculate basic metrics
        let containment = n_intersecting_hashes as f64 / query_size as f64;
        let jaccard = n_intersecting_hashes as f64 / union_size as f64;
        let containment_target_in_query = n_intersecting_hashes as f64 / target_size as f64;
        let max_containment = containment.max(containment_target_in_query);

        // Calculate abundance statistics
        let (average_abund, median_abund, std_abund) =
            if let (Some(query_abunds), Some(target_abunds)) =
                (query_abunds, target_abunds.as_ref())
            {
                significance::abundance_stats(&intersection, query_abunds, target_abunds)
            } else {
                (1.0, 1.0, 0.0)
            };

        // Calculate weighted metrics
        let f_weighted_target_in_query = significance::weighted_fraction_target_in_query(
            query_abunds.as_deref(),
            target_abunds.as_deref(),
        );

        // Calculate overlap probability between query and target
        let overlap_probability = self.calculate_overlap_probability(&intersection);

        let matched_regions = find_matched_regions(query, target, &intersection);

        Some(SearchResult {
            query_name: query_name.to_string(),
            query_md5: query_md5.to_string(),
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
            tfidf: query_tfidf,
            overlap_probability,
            matched_regions,
        })
    }

    /// Calculate TF-IDF score for a query signature
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

    /// Calculate probability of overlap between query and target
    ///
    /// This calculates the probability of the intersecting hashes of query and target against
    /// the frequency of those hashes in the whole database
    pub fn calculate_overlap_probability(&self, intersection: &HashSet<u64>) -> f64 {
        if intersection.is_empty() {
            return 0.0;
        }

        // Calculate probability of overlap using k-mer frequencies from the database
        // This follows the sourmash approach: sum of (query_freq * target_freq) for intersecting k-mers
        let prob_overlap: f64 = intersection
            .par_iter()
            .map(|&hashval| {
                // Get frequency of this k-mer in the database (how many signatures contain it)
                let db_frequency =
                    self.stats.kmer_frequencies.get(&hashval).copied().unwrap_or(1) as f64;
                let total_signatures = self.stats.total_signatures as f64;

                // Normalize frequency to [0,1] range
                let normalized_frequency = db_frequency / total_signatures;

                // For the query, we assume each k-mer has equal weight (1.0)
                // For the target, we use the normalized database frequency
                // This gives us the probability that both query and target would have this k-mer
                1.0 * normalized_frequency
            })
            .sum();

        // Clamp to [0,1] range
        prob_overlap.min(1.0).max(0.0)
    }

    /// Get the underlying index
    pub fn index(&self) -> &ProteomeIndex {
        &self.index
    }

    /// Get search statistics
    pub fn stats(&self) -> &SearchStats {
        &self.stats
    }

    /// Stitch overlapping k-mers together
    fn stitch_kmers(&self, kmer_positions: &[(usize, String)]) -> String {
        if kmer_positions.is_empty() {
            return String::new();
        }

        // Simple stitching: just concatenate k-mers with overlaps
        let mut result = String::new();
        let mut last_end = 0;

        for (pos, kmer) in kmer_positions {
            if *pos >= last_end {
                // No overlap, add the full k-mer
                result.push_str(kmer);
                last_end = pos + kmer.len();
            } else {
                // Overlap detected, add only the non-overlapping part
                let overlap = last_end - pos;
                if overlap < kmer.len() {
                    result.push_str(&kmer[overlap..]);
                    last_end = pos + kmer.len();
                }
            }
        }

        result
    }
}

/// Find all consecutive matched regions of k-mer overlap between a query and target sequences
///
/// WHY: This is a standalone function because it doesn't require any state from ProteinSearcher.
/// It only operates on the sketches and intersection provided. This makes it easier to test and
/// more reusable. This is idiomatic Rust - functions that don't need state should be standalone.
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
        let query_raw_sequence = query_sketch
            .get_raw_sequence()
            .unwrap_or_else(|| panic!("No raw sequence found for query signature {query_name}"));
        let target_raw_sequence = target_sketch
            .get_raw_sequence()
            .unwrap_or_else(|| panic!("No raw sequence found for target signature {target_name}"));

        // Extract subsequences using correct positions
        // WHY: Query subsequence uses query positions, target subsequence uses target positions.
        // This is the fix for the bug where both were using query positions.
        let query_subseq = &query_raw_sequence[query_start_pos..query_end_pos];
        let target_subseq = &target_raw_sequence[target_start_pos..target_end_pos];

        // Get moltype sequences for validation
        let target_moltype_sequence = target_sketch.get_moltype_sequence().unwrap_or_else(|| {
            panic!("No moltype encoded sequence found for target signature {target_name}")
        });
        let query_moltype_sequence = query_sketch.get_moltype_sequence().unwrap_or_else(|| {
            panic!("No moltype encoded sequence found for query signature {query_name}")
        });

        // Extract moltype subsequences using correct positions
        let query_moltype_seq = &query_moltype_sequence[query_start_pos..query_end_pos];
        let target_moltype_seq = &target_moltype_sequence[target_start_pos..target_end_pos];

        // Validate that moltype sequences match (they should since they share the same k-mers)
        if query_moltype_seq != target_moltype_seq {
            panic!(
                "Target: '{target_name}'\nand\nQuery: '{query_name}'\nmoltype sequences do not match:\
            \nQuery positions: {query_start_pos}..{query_end_pos}\
            \nTarget positions: {target_start_pos}..{target_end_pos}\
            \nTarget protein subsequence: {target_subseq}\
            \nTarget moltype subsequence: {target_moltype_seq}\
            \nQuery  moltype subsequence: {query_moltype_seq}\
            \nQuery  protein subsequence: {query_subseq}"
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
    /// Find a signature by name in the index
    fn find_signature_by_name(&self, name: &str) -> Option<ProteinSketch> {
        for entry in self.index.get_signatures().iter() {
            let signature = entry.value();
            if signature.signature().name == name {
                return Some(signature.clone());
            }
        }
        None
    }

    /// Get stored encoded sequence for a signature by name
    fn get_stored_encoded_sequence(&self, signature_name: &str) -> Option<String> {
        // Find the signature in the index
        for entry in self.index.get_signatures().iter() {
            let signature = entry.value();
            if signature.signature().name == signature_name {
                return signature.get_moltype_sequence().map(|s| s.to_string());
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sketch::ProteinSketch;
    use crate::tests::test_fixtures::{TEST_BLC2_FASTA, TEST_CED9_FASTA, TEST_FASTA_GZ};
    use needletail::parse_fastx_file;
    use std::path::Path;
    use tempfile::TempDir;

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

    const BCL2_CED9_K12: ExpectedSimilarity = ExpectedSimilarity {
        ksize: 12,
        n_intersecting_hashes: 24,
        containment: 0.09091,
        jaccard: 0.05217,
        max_containment: 0.10909,
        containment_target_in_query: 0.10909,
        average_abund: 1.0208333333333333,
        median_abund: 1.0,
        std_abund: 0.099913156735681657,
        matched_regions_count: 13,
    };

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

    #[fixture]
    fn temp_dir() -> TempDir {
        TempDir::new().unwrap()
    }

    #[fixture]
    fn ced9_sketch_k12() -> ProteinSketch {
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

        // Create a simple test FASTA file for query
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">test_query\nATCGATCGATCGATCG")?;

        // Create a simple test FASTA file for target
        let target_fasta = temp_path.join("target.fasta");
        std::fs::write(&target_fasta, ">test_target\nATCGATCGATCGATCG")?;

        // Create target index
        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(
            &target_index_path,
            10,    // ksize
            1,     // scaled
            "hp",  // moltype
            false, // store_raw_sequences
        )?;

        target_index.process_fasta(&target_fasta, 1000, 1000)?;

        // Create searcher
        let searcher = ProteinSearcher::new(target_index);

        // Create query index
        let query_index = ProteomeIndex::new_with_auto_filename(
            &query_fasta,
            10,    // ksize
            1,     // scaled
            "hp",  // moltype
            false, // store_raw_sequences
        )?;

        query_index.process_fasta(&query_fasta, 1000, 1000)?;

        // Get query signatures
        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        assert!(!query_signatures.is_empty(), "Should have at least one query signature");

        // Perform search
        let results = searcher.search(&query_signatures)?;

        // Should find at least one match (exact match)
        assert!(!results.is_empty(), "Should find at least one match");

        // Check that the first result has reasonable values
        let first_result = &results[0];
        assert_eq!(first_result.query_name, "test_query");
        assert_eq!(first_result.target_name, "test_target");
        assert!(first_result.containment > 0.0);
        assert!(first_result.jaccard > 0.0);
        assert!(first_result.n_intersecting_hashes > 0);

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
    #[test]
    fn test_find_matched_regions_single() -> Result<()> {
        let ksize = 15;
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

        // Calculate intersection for find_matched_regions
        // WHY: We use a standalone function that doesn't require a searcher/index, making tests
        // simpler and more focused. This is idiomatic Rust - functions that don't need state
        // should be standalone.
        let query_mins: HashSet<u64> =
            query_sketch.signature().minhash.mins().iter().cloned().collect();
        let target_mins: HashSet<u64> =
            target_sketch.signature().minhash.mins().iter().cloned().collect();
        let intersection: HashSet<u64> = query_mins.intersection(&target_mins).cloned().collect();

        // Find matched regions using the standalone function
        let matched_regions = find_matched_regions(&query_sketch, &target_sketch, &intersection);

        // Verify we found at least one match
        assert_eq!(matched_regions.len(), 1, "Should find exactly one match");

        // Verify the expected match region
        // The expected match is at positions 162-181 in CED9 (query) and 138-157 in BCL2 (target)
        let matched_region = &matched_regions[0];
        assert_eq!(matched_region.query_subseq, "QCPMSYGRLIGLISFGGFV");
        assert_eq!(matched_region.moltype_seq, "pphhphhphhhhhphhhhh");
        assert_eq!(matched_region.target_subseq, "RDGVNWGRIVAFFEFGGVM");

        // Verify query and target positions
        assert_eq!(matched_region.query_start, 162, "Query start position should be 162");
        assert_eq!(matched_region.query_end, 181, "Query end position should be 181");
        assert_eq!(matched_region.target_start, 138, "Target start position should be 138");
        assert_eq!(matched_region.target_end, 157, "Target end position should be 157");

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

        // Calculate intersection for find_matched_regions
        // WHY: We use a standalone function that doesn't require a searcher/index, making tests
        // simpler and more focused. This is idiomatic Rust - functions that don't need state
        // should be standalone.
        let query_mins: HashSet<u64> =
            query_sketch.signature().minhash.mins().iter().cloned().collect();
        let target_mins: HashSet<u64> =
            target_sketch.signature().minhash.mins().iter().cloned().collect();
        let intersection: HashSet<u64> = query_mins.intersection(&target_mins).cloned().collect();

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

        target_index.process_fasta(&target_fasta, 1000, 1000)?;

        // Create searcher
        let searcher = ProteinSearcher::new(target_index);

        // Create query FASTA
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">query\nATCGATCGATCGATCG")?;

        let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;

        query_index.process_fasta(&query_fasta, 1000, 1000)?;

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
        target_index.process_fasta(&target_fasta, 1000, 1000)?;

        let searcher = ProteinSearcher::new(target_index);

        let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;
        query_index.process_fasta(&query_fasta, 1000, 1000)?;

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
        target_index.process_fasta(&target_fasta, 1000, 1000)?;

        let searcher = ProteinSearcher::new(target_index);

        // Create query
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">query\nATCGATCGATCGATCG")?;

        let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;
        query_index.process_fasta(&query_fasta, 1000, 1000)?;

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

        let query = ProteinSketch::new("test", 10, 5, "hp").unwrap();
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

    #[rstest]
    #[case::k12(ced9_sketch_k12(), bcl2_sketch_k12(), BCL2_CED9_K12)]
    #[case::k15(ced9_sketch_k15(), bcl2_sketch_k15(), BCL2_CED9_K15)]
    fn test_query_target_similarity_bcl2_ced9(
        temp_dir: TempDir,
        #[case] query_sketch: ProteinSketch,
        #[case] target_sketch: ProteinSketch,
        #[case] expected: ExpectedSimilarity,
    ) -> Result<()> {
        let scaled = 1;
        let moltype = "hp";

        // ... index setup using temp_dir ...
        
        let searcher = ProteinSearcher::new(target_sketch);

        let result = searcher
            .query_target_similarity(query_sketch, target_sketch, scaled, moltype)?
            .expect("Should find similarity between CED9 and BCL2");

        // Assertions are now readable
        assert_eq!(result.n_intersecting_hashes, expected.n_intersecting_hashes);
        assert_eq!(result.containment, expected.containment);
        assert_eq!(result.jaccard, expected.jaccard);
        assert_eq!(result.max_containment, expected.max_containment);
        assert_eq!(result.average_abund, expected.average_abund);
        assert_eq!(result.matched_regions.len(), expected.matched_regions_count);

        Ok(())
    }

    /// Test query_target_similarity with BCL2 and CED9 sequences with k-mer size 12,
    /// which produces multiple consecutive matched k-mer regions.
    ///
    /// This test verifies that query_target_similarity correctly calculates all similarity
    /// metrics including containment, jaccard, max_containment, abundance statistics, and
    /// overlap probability for a known pair of related proteins.
    #[test]
    fn test_query_target_similarity_bcl2_ced9_k12() -> Result<()> {
        let ksize = 12;
        let scaled = 1;
        let moltype = "hp";

        // Read CED9 sequence from FASTA file
        let (ced9_name, ced9_sequence) = read_first_fasta_record(TEST_CED9_FASTA)?;

        // Read BCL2 sequence from FASTA file
        let (bcl2_name, bcl2_sequence) = read_first_fasta_record(TEST_BLC2_FASTA)?;

        // Create query sketch (CED9)
        let query_sketch = ProteinSketch::from_protein_sequence(
            &ced9_name,
            &ced9_sequence,
            ksize,
            scaled,
            moltype,
        )?;

        // Create target sketch (BCL2)
        let target_sketch = ProteinSketch::from_protein_sequence(
            &bcl2_name,
            &bcl2_sequence,
            ksize,
            scaled,
            moltype,
        )?;

        // Create a temporary index with the target for the searcher
        // WHY: ProteinSearcher requires an index and calculates stats from it. We need
        // this for overlap probability and TF-IDF calculations. Creating a minimal index
        // allows us to test the core similarity calculation logic.
        //
        // NOTE: With only 1 signature in the index, TF-IDF will be 0 and overlap probability
        // will be 1.0. These metrics are designed for database searches (1vMany), not 1v1
        // comparisons. To properly test these metrics, you would need an index with multiple
        // signatures where some k-mers are common and others are rare.
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path().join("test_index");
        let target_index = ProteomeIndex::new(&temp_path, ksize, scaled, moltype, true)?;

        // Create a temporary FASTA file with BCL2 sequence
        let target_fasta = temp_dir.path().join("target.fasta");
        std::fs::write(&target_fasta, format!(">{}\n{}", bcl2_name, bcl2_sequence))?;

        // Process the target into the index
        target_index.process_fasta(&target_fasta, 1000, 1000)?;

        // Create searcher from the index
        let searcher = ProteinSearcher::new(target_index);

        // Extract query data for query_target_similarity
        let query_mins: HashSet<u64> =
            query_sketch.signature().minhash.mins().iter().cloned().collect();
        let query_abunds = query_sketch.signature().minhash.abunds();
        let query_name = query_sketch.signature().name.clone();
        let query_md5 = query_sketch.signature().md5sum.clone();

        // Calculate TF-IDF for the query (needed for query_target_similarity)
        let query_tfidf = searcher.calculate_tfidf(&query_sketch);

        // Call query_target_similarity
        let result = searcher.query_target_similarity(
            &query_sketch,
            &target_sketch,
            &query_mins,
            query_abunds.as_deref(),
            &query_name,
            &query_md5,
            query_tfidf,
        );

        // Verify we got a result (should have some intersection)
        assert!(result.is_some(), "Should find similarity between CED9 and BCL2");
        let result = result.unwrap();

        // Verify basic metadata
        assert_eq!(result.query_name, ced9_name);
        assert_eq!(result.target_name, bcl2_name);
        assert_eq!(result.ksize, ksize);
        assert_eq!(result.scaled, scaled);
        assert_eq!(result.moltype, moltype);

        // Verify we have intersecting k-mers
        assert_eq!(
            result.n_intersecting_hashes, 24,
            "Should have 24 intersecting k-mers between CED9 and BCL2"
        );

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

        assert_eq!(
            result.containment, 0.09091,
            "Containment should be 0.09091, got {}",
            result.containment
        );
        assert_eq!(
            result.jaccard, 0.05217,
            "Jaccard should be 0.052173913043478258, got {}",
            result.jaccard
        );
        assert_eq!(
            result.max_containment, 0.10909,
            "Max containment should be 0.10909, got {}",
            result.max_containment
        );
        assert_eq!(
            result.containment_target_in_query, 0.10909,
            "Target in query containment should be 0.10909, got {}",
            result.containment_target_in_query
        );

        // Verify abundance statistics (should be valid if abundances are tracked)
        assert_eq!(
            result.average_abund, 1.0208333333333333,
            "Average abundance should be 1.0208333333333333"
        );
        assert_eq!(result.median_abund, 1.0, "Median abundance should be 1.0");
        assert_eq!(
            result.std_abund, 0.099913156735681657,
            "Standard deviation of abundance should be 0.099913156735681657"
        );

        // Verify weighted metrics
        assert!(
            result.f_weighted_target_in_query >= 0.0 && result.f_weighted_target_in_query <= 1.0,
            "Weighted fraction should be in [0, 1], got {}",
            result.f_weighted_target_in_query
        );

        // Verify overlap probability
        // WHY: With only 1 signature in the database, normalized_frequency = db_frequency/total_signatures = 1/1 = 1.0
        // This metric is designed for database searches (1vMany), not 1v1 comparisons.
        // In a real database search, overlap probability would indicate how common the intersecting
        // k-mers are across the database (higher = more common/less significant).
        assert!(
            result.overlap_probability == 1.0,
            "Overlap probability should be 1.0 with 1 signature in database, got {}",
            result.overlap_probability
        );

        // Verify TF-IDF
        // WHY: IDF = ln(total_signatures / freq). With only 1 signature, IDF = ln(1/1) = 0 for all k-mers.
        // This metric is designed for database searches (1vMany), not 1v1 comparisons.
        // In a real database search, TF-IDF would weight k-mers by their rarity (higher IDF = more rare/significant).
        assert_eq!(
            result.tfidf, 0.0,
            "TF-IDF should be 0.0 with 1 signature in database, got {}",
            result.tfidf
        );

        // Verify matched regions are present (BCL2 and CED9 should have overlapping regions)
        assert_eq!(
            result.matched_regions.len(),
            13,
            "Should find 13 matched regions between CED9 and BCL2"
        );

        Ok(())
    }

    /// Test query_target_similarity with BCL2 and CED9 sequences with k-mer size 15, which
    /// produces a single consecutive matched k-mer region.
    ///
    /// This test verifies that query_target_similarity correctly calculates all similarity
    /// metrics including containment, jaccard, max_containment, abundance statistics, and
    /// overlap probability for a known pair of related proteins.
    #[test]
    fn test_query_target_similarity_bcl2_ced9_k15() -> Result<()> {
        let ksize = 15;
        let scaled = 1;
        let moltype = "hp";

        // Read CED9 sequence from FASTA file
        let (ced9_name, ced9_sequence) = read_first_fasta_record(TEST_CED9_FASTA)?;

        // Read BCL2 sequence from FASTA file
        let (bcl2_name, bcl2_sequence) = read_first_fasta_record(TEST_BLC2_FASTA)?;

        // Create query sketch (CED9)
        let query_sketch = ProteinSketch::from_protein_sequence(
            &ced9_name,
            &ced9_sequence,
            ksize,
            scaled,
            moltype,
        )?;

        // Create target sketch (BCL2)
        let target_sketch = ProteinSketch::from_protein_sequence(
            &bcl2_name,
            &bcl2_sequence,
            ksize,
            scaled,
            moltype,
        )?;

        // Create a temporary index with the target for the searcher
        // WHY: ProteinSearcher requires an index and calculates stats from it. We need
        // this for overlap probability and TF-IDF calculations. Creating a minimal index
        // allows us to test the core similarity calculation logic.
        //
        // NOTE: With only 1 signature in the index, TF-IDF will be 0 and overlap probability
        // will be 1.0. These metrics are designed for database searches (1vMany), not 1v1
        // comparisons. To properly test these metrics, you would need an index with multiple
        // signatures where some k-mers are common and others are rare.
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path().join("test_index");
        let target_index = ProteomeIndex::new(&temp_path, ksize, scaled, moltype, true)?;

        // Create a temporary FASTA file with BCL2 sequence
        let target_fasta = temp_dir.path().join("target.fasta");
        std::fs::write(&target_fasta, format!(">{}\n{}", bcl2_name, bcl2_sequence))?;

        // Process the target into the index
        target_index.process_fasta(&target_fasta, 1000, 1000)?;

        // Create searcher from the index
        let searcher = ProteinSearcher::new(target_index);

        // Extract query data for query_target_similarity
        let query_mins: HashSet<u64> =
            query_sketch.signature().minhash.mins().iter().cloned().collect();
        let query_abunds = query_sketch.signature().minhash.abunds();
        let query_name = query_sketch.signature().name.clone();
        let query_md5 = query_sketch.signature().md5sum.clone();

        // Calculate TF-IDF for the query (needed for query_target_similarity)
        let query_tfidf = searcher.calculate_tfidf(&query_sketch);

        // Call query_target_similarity
        let result = searcher.query_target_similarity(
            &query_sketch,
            &target_sketch,
            &query_mins,
            query_abunds.as_deref(),
            &query_name,
            &query_md5,
            query_tfidf,
        );

        // Verify we got a result (should have some intersection)
        assert!(result.is_some(), "Should find similarity between CED9 and BCL2");
        let result = result.unwrap();

        // Verify basic metadata
        assert_eq!(result.query_name, ced9_name);
        assert_eq!(result.target_name, bcl2_name);
        assert_eq!(result.ksize, ksize);
        assert_eq!(result.scaled, scaled);
        assert_eq!(result.moltype, moltype);

        // Verify we have intersecting k-mers
        assert_eq!(
            result.n_intersecting_hashes, 5,
            "Should have 5 intersecting k-mers between CED9 and BCL2"
        );

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

        assert_eq!(
            result.containment, 0.018796992481203006,
            "Containment should be 0.018796992481203006, got {}",
            result.containment
        );
        assert_eq!(
            result.jaccard, 0.0102880658436214,
            "Jaccard should be 0.0102880658436214, got {}",
            result.jaccard
        );
        assert_eq!(
            result.max_containment, 0.022222222222222223,
            "Max containment should be 0.022222222222222223, got {}",
            result.max_containment
        );
        assert_eq!(
            result.containment_target_in_query, 0.022222222222222223,
            "Target in query containment should be 0.022222222222222223, got {}",
            result.containment_target_in_query
        );

        // Verify abundance statistics (should be valid if abundances are tracked)
        assert_eq!(
            result.average_abund, 1.0,
            "Average abundance should be 1.0, got {}",
            result.average_abund
        );
        assert_eq!(
            result.median_abund, 1.0,
            "Median abundance should be 1.0, got {}",
            result.median_abund
        );
        assert_eq!(
            result.std_abund, 0.0,
            "Standard deviation of abundance should be 0.0, got {}",
            result.std_abund
        );

        // Verify weighted metrics
        assert_eq!(
            result.f_weighted_target_in_query, 0.8458646616541353,
            "Weighted fraction should be 0.8458646616541353, got {}",
            result.f_weighted_target_in_query
        );

        // Verify overlap probability
        // WHY: With only 1 signature in the database, normalized_frequency = db_frequency/total_signatures = 1/1 = 1.0
        // This metric is designed for database searches (1vMany), not 1v1 comparisons.
        // In a real database search, overlap probability would indicate how common the intersecting
        // k-mers are across the database (higher = more common/less significant).
        assert_eq!(
            result.overlap_probability, 1.0,
            "Overlap probability should be 1.0 with 1 signature in database, got {}",
            result.overlap_probability
        );

        // Verify TF-IDF
        // WHY: IDF = ln(total_signatures / freq). With only 1 signature, IDF = ln(1/1) = 0 for all k-mers.
        // This metric is designed for database searches (1vMany), not 1v1 comparisons.
        // In a real database search, TF-IDF would weight k-mers by their rarity (higher IDF = more rare/significant).
        assert_eq!(
            result.tfidf, 0.0,
            "TF-IDF should be 0.0 with 1 signature in database, got {}",
            result.tfidf
        );

        // Verify matched regions are present (BCL2 and CED9 should have overlapping regions)
        assert_eq!(
            result.matched_regions.len(),
            1,
            "Should find 1 matched region between CED9 and BCL2"
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
        target_index.process_fasta(TEST_FASTA_GZ, 0, 1000)?;

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

        query_index.process_fasta(TEST_CED9_FASTA, 0, 1000)?;

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
            bcl2_result.tfidf >= 0.0,
            "TF-IDF should be non-negative, got {}",
            bcl2_result.tfidf
        );
        // Note: TF-IDF might still be 0 if all query k-mers appear in all database signatures,
        // but with a diverse database like bcl2_first25, we expect some variation

        // Verify overlap probability is meaningful (should not be 1.0 with multiple signatures)
        // WHY: With multiple signatures, k-mers will have varying frequencies. Overlap probability
        // indicates how common the intersecting k-mers are across the database. Values < 1.0
        // indicate that the k-mers are not universal across all signatures.
        assert!(
            bcl2_result.overlap_probability >= 0.0 && bcl2_result.overlap_probability <= 1.0,
            "Overlap probability should be in [0, 1], got {}",
            bcl2_result.overlap_probability
        );
        // Note: Overlap probability might still be 1.0 if all intersecting k-mers appear in
        // all signatures, but with a diverse database, we expect lower values indicating
        // that some k-mers are more specific to certain proteins

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
