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
    fn query_target_similarity(
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
    use crate::tests::test_fixtures::{TEST_BLC2_FASTA, TEST_CED9_FASTA};
    use needletail::parse_fastx_file;
    use std::path::Path;
    use tempfile::TempDir;

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
}
