use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::Result;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::errors::IndexResult;
use crate::index::ProteomeIndex;
use crate::signature::ProteinSignature;

/// Search result for a single query-target pair
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SearchResult {
    /// Query sequence name
    pub query_name: String,
    /// Query sequence MD5 hash
    pub query_md5: String,
    /// Target sequence name
    pub match_name: String,
    /// Containment score (intersection / query_size)
    pub containment: f64,
    /// Number of intersecting k-mers
    pub intersect_hashes: usize,
    /// K-mer size used
    pub ksize: u32,
    /// Scaled factor used
    pub scaled: u32,
    /// Molecular type (hp, dayhoff, protein)
    pub moltype: String,
    /// Target sequence MD5 hash
    pub match_md5: String,
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
    /// Query containment ANI (Average Nucleotide Identity)
    pub query_containment_ani: f64,
    /// Match containment ANI
    pub match_containment_ani: f64,
    /// Average containment ANI
    pub average_containment_ani: f64,
    /// Maximum containment ANI
    pub max_containment_ani: f64,
    /// Number of weighted found k-mers
    pub n_weighted_found: usize,
    /// Total weighted hashes
    pub total_weighted_hashes: usize,
    /// Containment of target in query
    pub containment_target_in_query: f64,
    /// Weighted fraction of target in query
    pub f_weighted_target_in_query: f64,
}

/// Detailed search result with k-mer information for stitched output
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetailedSearchResult {
    /// Target sequence name
    pub match_name: String,
    /// Query sequence name
    pub query_name: String,
    /// Query sequence start position
    pub query_start: u32,
    /// Query sequence end position
    pub query_end: u32,
    /// Query sequence (stitched k-mers)
    pub query: String,
    /// Target sequence start position
    pub match_start: u32,
    /// Target sequence end position
    pub match_end: u32,
    /// Target sequence (stitched k-mers)
    pub r#match: String,
    /// Encoded sequence (hp/dayhoff/protein encoding)
    pub encoded: String,
    /// Length of the match
    pub length: u32,
    /// Formatted output for stderr display
    pub to_print: String,
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

/// Represents a matching region between query and target sequences
#[derive(Debug, Clone)]
struct MatchingRegion {
    query_start: usize,
    query_end: usize,
    match_start: usize,
    match_end: usize,
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

    /// Search a single query signature against the database
    pub fn search_single(&self, query: &ProteinSignature) -> Result<Vec<SearchResult>> {
        let query_mins: HashSet<u64> = query.signature().minhash.mins().iter().cloned().collect();
        let query_abunds = query.signature().minhash.abunds();
        let query_name = query.signature().name.clone();
        let query_md5 = query.signature().md5sum.clone();

        let results: Vec<SearchResult> = self
            .index
            .get_signatures()
            .iter()
            .filter_map(|entry| {
                let target = entry.value();
                self.calculate_similarity(
                    query,
                    target,
                    &query_mins,
                    query_abunds.as_deref(),
                    &query_name,
                    &query_md5,
                )
            })
            .collect();

        // Sort by containment score (descending)
        let mut sorted_results = results;
        sorted_results.sort_by(|a, b| {
            b.containment.partial_cmp(&a.containment).unwrap_or(std::cmp::Ordering::Equal)
        });

        Ok(sorted_results)
    }

    /// Search multiple query signatures against the database
    pub fn search_multiple(&self, queries: &[ProteinSignature]) -> Result<Vec<SearchResult>> {
        let all_results: Vec<SearchResult> = queries
            .par_iter()
            .flat_map(|query| self.search_single(query).unwrap_or_default())
            .collect();

        // Sort by containment score (descending)
        let mut sorted_results = all_results;
        sorted_results.sort_by(|a, b| {
            b.containment.partial_cmp(&a.containment).unwrap_or(std::cmp::Ordering::Equal)
        });

        Ok(sorted_results)
    }

    /// Calculate similarity between query and target signatures
    fn calculate_similarity(
        &self,
        query: &ProteinSignature,
        target: &ProteinSignature,
        query_mins: &HashSet<u64>,
        query_abunds: Option<&[u64]>,
        query_name: &str,
        query_md5: &str,
    ) -> Option<SearchResult> {
        let target_mins: HashSet<u64> = target.signature().minhash.mins().iter().cloned().collect();
        let target_abunds = target.signature().minhash.abunds();

        // Calculate intersection
        let intersection: HashSet<u64> = query_mins.intersection(&target_mins).cloned().collect();
        let intersect_hashes = intersection.len();

        // Skip if no intersection
        if intersect_hashes == 0 {
            return None;
        }

        let query_size = query_mins.len();
        let target_size = target_mins.len();
        let union_size = query_size + target_size - intersect_hashes;

        // Calculate basic metrics
        let containment = intersect_hashes as f64 / query_size as f64;
        let jaccard = intersect_hashes as f64 / union_size as f64;
        let containment_target_in_query = intersect_hashes as f64 / target_size as f64;
        let max_containment = containment.max(containment_target_in_query);

        // Calculate abundance statistics
        let (average_abund, median_abund, std_abund) =
            if let (Some(query_abunds), Some(target_abunds)) =
                (query_abunds, target_abunds.as_ref())
            {
                self.calculate_abundance_stats(&intersection, query_abunds, target_abunds)
            } else {
                (1.0, 1.0, 0.0)
            };

        // Calculate ANI (Average Nucleotide Identity) - simplified version
        let query_containment_ani = self.calculate_ani(containment, query_size);
        let match_containment_ani = self.calculate_ani(containment_target_in_query, target_size);
        let average_containment_ani = (query_containment_ani + match_containment_ani) / 2.0;
        let max_containment_ani = query_containment_ani.max(match_containment_ani);

        // Calculate weighted metrics
        let (n_weighted_found, total_weighted_hashes, f_weighted_target_in_query) =
            self.calculate_weighted_metrics(&intersection, query_abunds, target_abunds.as_deref());

        Some(SearchResult {
            query_name: query_name.to_string(),
            query_md5: query_md5.to_string(),
            match_name: target.signature().name.clone(),
            containment,
            intersect_hashes,
            ksize: query.protein_ksize(),
            scaled: query.signature().minhash.scaled(),
            moltype: query.moltype().to_string(),
            match_md5: target.signature().md5sum.clone(),
            jaccard,
            max_containment,
            average_abund,
            median_abund,
            std_abund,
            query_containment_ani,
            match_containment_ani,
            average_containment_ani,
            max_containment_ani,
            n_weighted_found,
            total_weighted_hashes,
            containment_target_in_query,
            f_weighted_target_in_query,
        })
    }

    /// Calculate abundance statistics for intersecting k-mers
    fn calculate_abundance_stats(
        &self,
        intersection: &HashSet<u64>,
        query_abunds: &[u64],
        target_abunds: &[u64],
    ) -> (f64, f64, f64) {
        let mut abunds = Vec::new();

        // Get abundances for intersecting k-mers
        for (i, &_min) in intersection.iter().enumerate() {
            if let (Some(&query_abund), Some(&target_abund)) =
                (query_abunds.get(i), target_abunds.get(i))
            {
                abunds.push((query_abund + target_abund) as f64 / 2.0);
            }
        }

        if abunds.is_empty() {
            return (1.0, 1.0, 0.0);
        }

        // Calculate statistics
        let sum: f64 = abunds.iter().sum();
        let average = sum / abunds.len() as f64;

        abunds.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median = if abunds.len() % 2 == 0 {
            (abunds[abunds.len() / 2 - 1] + abunds[abunds.len() / 2]) / 2.0
        } else {
            abunds[abunds.len() / 2]
        };

        let variance: f64 =
            abunds.iter().map(|&x| (x - average).powi(2)).sum::<f64>() / abunds.len() as f64;
        let std_dev = variance.sqrt();

        (average, median, std_dev)
    }

    /// Calculate ANI (Average Nucleotide Identity) from containment
    fn calculate_ani(&self, containment: f64, _size: usize) -> f64 {
        // Simplified ANI calculation based on containment
        // This is a rough approximation - in practice, ANI calculation is more complex
        if containment <= 0.0 {
            0.0
        } else {
            // Use a logarithmic relationship for ANI
            let ani = 1.0 - (-containment.ln()).exp();
            ani.clamp(0.0, 1.0)
        }
    }

    /// Calculate weighted metrics
    fn calculate_weighted_metrics(
        &self,
        intersection: &HashSet<u64>,
        query_abunds: Option<&[u64]>,
        target_abunds: Option<&[u64]>,
    ) -> (usize, usize, f64) {
        let n_weighted_found = intersection.len();
        let total_weighted_hashes = intersection.len(); // Simplified

        let f_weighted_target_in_query =
            if let (Some(query_abunds), Some(target_abunds)) = (query_abunds, target_abunds) {
                let query_weight: f64 = query_abunds.iter().sum::<u64>() as f64;
                let target_weight: f64 = target_abunds.iter().sum::<u64>() as f64;
                if query_weight > 0.0 {
                    target_weight / query_weight
                } else {
                    0.0
                }
            } else {
                n_weighted_found as f64 / total_weighted_hashes as f64
            };

        (n_weighted_found, total_weighted_hashes, f_weighted_target_in_query)
    }

    /// Calculate TF-IDF score for a query signature
    pub fn calculate_tfidf(&self, query: &ProteinSignature) -> f64 {
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
    pub fn calculate_overlap_probability(
        &self,
        query: &ProteinSignature,
        target: &ProteinSignature,
    ) -> f64 {
        let query_size = query.signature().minhash.mins().len();
        let target_size = target.signature().minhash.mins().len();
        let total_kmers = self.stats.total_signatures;

        // Simplified probability calculation
        // This is a rough approximation - actual probability calculation would be more complex
        if query_size == 0 || target_size == 0 {
            return 0.0;
        }

        let expected_overlap = (query_size * target_size) as f64 / total_kmers as f64;
        let actual_overlap = self.calculate_intersection_size(query, target) as f64;

        if expected_overlap <= 0.0 {
            0.0
        } else {
            // Use a Poisson-like probability calculation
            let ratio = actual_overlap / expected_overlap;
            if ratio > 1.0 {
                1.0 - (-ratio).exp()
            } else {
                ratio
            }
        }
    }

    /// Calculate the size of intersection between two signatures
    fn calculate_intersection_size(
        &self,
        query: &ProteinSignature,
        target: &ProteinSignature,
    ) -> usize {
        let query_mins: HashSet<u64> = query.signature().minhash.mins().iter().cloned().collect();
        let target_mins: HashSet<u64> = target.signature().minhash.mins().iter().cloned().collect();

        query_mins.intersection(&target_mins).count()
    }

    /// Get the underlying index
    pub fn index(&self) -> &ProteomeIndex {
        &self.index
    }

    /// Get search statistics
    pub fn stats(&self) -> &SearchStats {
        &self.stats
    }

    /// Search with k-mer extraction and stitching (always enabled)
    pub fn search_with_kmer_extraction(
        &self,
        queries: &[ProteinSignature],
    ) -> Result<Vec<DetailedSearchResult>> {
        let mut detailed_results = Vec::new();

        for query in queries {
            let query_mins: HashSet<u64> =
                query.signature().minhash.mins().iter().cloned().collect();
            let _query_name = query.signature().name.clone();

            // Find matches for this query
            let matches: Vec<_> = self
                .index
                .get_signatures()
                .iter()
                .filter_map(|entry| {
                    let target = entry.value().clone();
                    let target_mins: HashSet<u64> =
                        target.signature().minhash.mins().iter().cloned().collect();
                    let intersection: HashSet<u64> =
                        query_mins.intersection(&target_mins).cloned().collect();

                    if intersection.is_empty() {
                        None
                    } else {
                        Some((target, intersection))
                    }
                })
                .collect();

            // For each match, create detailed results
            for (target, intersection) in matches {
                if let Some(detailed_result) =
                    self.create_detailed_result(query, &target, &intersection)
                {
                    detailed_results.push(detailed_result);
                }
            }
        }

        // Sort by query_start, query_end
        detailed_results.sort_by(|a, b| {
            a.query_start.cmp(&b.query_start).then_with(|| a.query_end.cmp(&b.query_end))
        });

        Ok(detailed_results)
    }

    /// Create a detailed search result with stitched k-mers
    fn create_detailed_result(
        &self,
        query: &ProteinSignature,
        target: &ProteinSignature,
        intersection: &HashSet<u64>,
    ) -> Option<DetailedSearchResult> {
        let query_name = query.signature().name.clone();
        let match_name = target.signature().name.clone();

        // Try to get raw sequences first (most accurate)
        if let (Some(query_seq), Some(target_seq)) =
            (query.get_raw_sequence(), target.get_raw_sequence())
        {
            return self.create_detailed_result_from_sequences(
                &query_name,
                &match_name,
                query_seq,
                target_seq,
                intersection,
                query,
            );
        }

        // Fall back to k-mer stitching
        self.create_detailed_result_from_kmers(
            &query_name,
            &match_name,
            query,
            target,
            intersection,
        )
    }

    /// Create detailed result from raw sequences (most accurate)
    fn create_detailed_result_from_sequences(
        &self,
        query_name: &str,
        match_name: &str,
        query_seq: &str,
        target_seq: &str,
        intersection: &HashSet<u64>,
        query_signature: &ProteinSignature,
    ) -> Option<DetailedSearchResult> {
        // Find the best matching region based on k-mer positions
        let matching_regions =
            self.find_matching_regions_with_signatures(query_signature, match_name, intersection)?;

        // Extract the matching regions from the sequences
        let query_region = &query_seq[matching_regions.query_start..matching_regions.query_end];
        let target_region = &target_seq[matching_regions.match_start..matching_regions.match_end];

        // Get the encoded sequence for the query region
        // Always generate the encoded sequence from the extracted query_region to ensure correct length
        let encoded_seq = self.encode_sequence_hp(query_region);

        let to_print = format!(
            "---\nQuery Name: {}\nMatch Name: {}\nquery: {} ({}-{})\nalpha: {}\nmatch: {} ({}-{})\n",
            query_name, match_name, query_region, matching_regions.query_start, matching_regions.query_end,
            encoded_seq, target_region, matching_regions.match_start, matching_regions.match_end
        );

        Some(DetailedSearchResult {
            match_name: match_name.to_string(),
            query_name: query_name.to_string(),
            query_start: matching_regions.query_start as u32,
            query_end: matching_regions.query_end as u32,
            query: query_region.to_string(),
            match_start: matching_regions.match_start as u32,
            match_end: matching_regions.match_end as u32,
            r#match: target_region.to_string(),
            encoded: encoded_seq,
            length: (matching_regions.query_end - matching_regions.query_start) as u32,
            to_print,
        })
    }

    /// Create detailed result by stitching k-mers together
    fn create_detailed_result_from_kmers(
        &self,
        query_name: &str,
        match_name: &str,
        query: &ProteinSignature,
        target: &ProteinSignature,
        intersection: &HashSet<u64>,
    ) -> Option<DetailedSearchResult> {
        // Get k-mer information for intersecting k-mers
        let mut query_kmers = Vec::new();
        let mut target_kmers = Vec::new();

        for &hashval in intersection {
            if let (Some(query_kmer_info), Some(target_kmer_info)) =
                (query.kmer_infos().get(&hashval), target.kmer_infos().get(&hashval))
            {
                // Collect k-mer positions - we'll extract the actual k-mers from raw sequences
                for &pos in &query_kmer_info.positions {
                    query_kmers.push((pos, String::new())); // Will be filled from raw sequence
                }
                for &pos in &target_kmer_info.positions {
                    target_kmers.push((pos, String::new())); // Will be filled from raw sequence
                }
            }
        }

        if query_kmers.is_empty() || target_kmers.is_empty() {
            return None;
        }

        // Sort by position
        query_kmers.sort_by_key(|(pos, _)| *pos);
        target_kmers.sort_by_key(|(pos, _)| *pos);

        // Stitch k-mers together (simplified approach)
        let query_stitched = self.stitch_kmers(&query_kmers);
        let target_stitched = self.stitch_kmers(&target_kmers);

        let query_start = query_kmers.first().map(|(pos, _)| *pos as u32).unwrap_or(0);
        let query_end =
            query_kmers.last().map(|(pos, kmer)| (*pos + kmer.len()) as u32).unwrap_or(0);
        let match_start = target_kmers.first().map(|(pos, _)| *pos as u32).unwrap_or(0);
        let match_end =
            target_kmers.last().map(|(pos, kmer)| (*pos + kmer.len()) as u32).unwrap_or(0);

        let length = query_end - query_start;

        // Try to get the stored encoded sequence first, otherwise generate it
        let encoded_seq = if let Some(stored_encoded) = self.get_stored_encoded_sequence(query_name)
        {
            // Use the stored encoded sequence for the matching region
            if query_start < stored_encoded.len() as u32 && query_end <= stored_encoded.len() as u32
            {
                stored_encoded[query_start as usize..query_end as usize].to_string()
            } else {
                self.encode_sequence_hp(&query_stitched)
            }
        } else {
            // Fallback: generate encoded sequence
            self.encode_sequence_hp(&query_stitched)
        };

        let to_print = format!(
            "---\nQuery Name: {}\nMatch Name: {}\nquery: {} ({}-{})\nalpha: {}\nmatch: {} ({}-{})\n",
            query_name, match_name, query_stitched, query_start, query_end, encoded_seq, target_stitched, match_start, match_end
        );

        Some(DetailedSearchResult {
            match_name: match_name.to_string(),
            query_name: query_name.to_string(),
            query_start,
            query_end,
            query: query_stitched,
            match_start,
            match_end,
            r#match: target_stitched,
            encoded: encoded_seq,
            length,
            to_print,
        })
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

    /// Find the best matching region based on k-mer positions with both signatures
    fn find_matching_regions_with_signatures(
        &self,
        query_signature: &ProteinSignature,
        match_name: &str,
        intersection: &HashSet<u64>,
    ) -> Option<MatchingRegion> {
        // Find the target signature in the index
        let target_sig = self.find_signature_by_name(match_name)?;

        // Get the k-mer size
        let ksize = query_signature.protein_ksize() as usize;

        // Collect all positions for intersecting k-mers from both signatures
        let mut query_positions = Vec::new();
        let mut target_positions = Vec::new();

        for &hashval in intersection {
            if let (Some(query_kmer_info), Some(target_kmer_info)) =
                (query_signature.kmer_infos().get(&hashval), target_sig.kmer_infos().get(&hashval))
            {
                query_positions.extend(&query_kmer_info.positions);
                target_positions.extend(&target_kmer_info.positions);
            }
        }

        if query_positions.is_empty() || target_positions.is_empty() {
            return None;
        }

        // Sort positions
        query_positions.sort();
        target_positions.sort();

        // Find the best matching region by looking for the longest contiguous sequence
        // of overlapping k-mers. For now, use a simpler approach: find the region
        // with the most k-mers in a reasonable window size.

        // Use a sliding window approach to find the best matching region
        let mut best_query_start = 0;
        let mut best_query_end = 0;
        let mut best_target_start = 0;
        let mut best_target_end = 0;
        let mut max_kmer_count = 0;

        // Try different window sizes to find the best match
        for window_size in [16, 32, 48, 64, 96, 128] {
            for &query_start in &query_positions {
                let query_end = query_start + window_size;

                // Count how many k-mers fall within this window
                let kmer_count = query_positions
                    .iter()
                    .filter(|&&pos| pos >= query_start && pos + ksize <= query_end)
                    .count();

                if kmer_count > max_kmer_count && kmer_count >= 2 {
                    max_kmer_count = kmer_count;
                    best_query_start = query_start;
                    best_query_end = query_end;

                    // Find corresponding target region
                    // For simplicity, use the first target position that overlaps
                    if let Some(&target_start) = target_positions.first() {
                        best_target_start = target_start;
                        best_target_end = target_start + window_size;
                    }
                }
            }
        }

        // If we found a good match, use it; otherwise fall back to the simple approach
        if max_kmer_count >= 2 {
            Some(MatchingRegion {
                query_start: best_query_start,
                query_end: best_query_end,
                match_start: best_target_start,
                match_end: best_target_end,
            })
        } else {
            // Fallback: use the first and last positions
            let query_start = *query_positions.first()?;
            let query_end = *query_positions.last()? + ksize;
            let target_start = *target_positions.first()?;
            let target_end = *target_positions.last()? + ksize;

            Some(MatchingRegion {
                query_start,
                query_end,
                match_start: target_start,
                match_end: target_end,
            })
        }
    }

    /// Find a signature by name in the index
    fn find_signature_by_name(&self, name: &str) -> Option<ProteinSignature> {
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
                return signature.get_encoded_sequence().map(|s| s.to_string());
            }
        }
        None
    }

    /// Encode a protein sequence using HP encoding (hydrophobic/polar)
    fn encode_sequence_hp(&self, sequence: &str) -> String {
        use sourmash::encodings::aa_to_hp;

        let mut encoded = String::with_capacity(sequence.len());

        for byte in sequence.bytes() {
            let hp_char = aa_to_hp(byte);
            encoded.push(match hp_char {
                b'h' => 'h',
                b'p' => 'p',
                _ => 'h', // Default to hydrophobic for unknown characters
            });
        }

        encoded
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::signature::ProteinSignature;
    use tempfile::TempDir;

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
        let results = searcher.search_multiple(&query_signatures)?;

        // Should find at least one match (exact match)
        assert!(!results.is_empty(), "Should find at least one match");

        // Check that the first result has reasonable values
        let first_result = &results[0];
        assert_eq!(first_result.query_name, "test_query");
        assert_eq!(first_result.match_name, "test_target");
        assert!(first_result.containment > 0.0);
        assert!(first_result.jaccard > 0.0);
        assert!(first_result.intersect_hashes > 0);

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

    /// Test overlap probability calculation
    #[test]
    fn test_overlap_probability() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let temp_path = temp_dir.path();

        // Create target FASTA
        let target_fasta = temp_path.join("target.fasta");
        std::fs::write(&target_fasta, ">target\nATCGATCGATCGATCG")?;

        let target_index_path = temp_path.join("target_index");
        let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp", false)?;

        target_index.process_fasta(&target_fasta, 1000, 1000)?;

        let searcher = ProteinSearcher::new(target_index);

        // Create query FASTA
        let query_fasta = temp_path.join("query.fasta");
        std::fs::write(&query_fasta, ">query\nATCGATCGATCGATCG")?;

        let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;

        query_index.process_fasta(&query_fasta, 1000, 1000)?;

        let query_signatures: Vec<_> =
            query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

        let target_signatures: Vec<_> =
            searcher.index().get_signatures().iter().map(|entry| entry.value().clone()).collect();

        // Calculate overlap probability
        let prob =
            searcher.calculate_overlap_probability(&query_signatures[0], &target_signatures[0]);
        assert!((0.0..=1.0).contains(&prob), "Probability should be between 0 and 1");

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

        let results = searcher.search_multiple(&query_signatures)?;

        if !results.is_empty() {
            let result = &results[0];

            // Test that all required fields are present and have reasonable values
            assert!(!result.query_name.is_empty());
            assert!(!result.query_md5.is_empty());
            assert!(!result.match_name.is_empty());
            assert!(!result.match_md5.is_empty());
            assert!(!result.moltype.is_empty());

            assert!(result.containment >= 0.0 && result.containment <= 1.0);
            assert!(result.jaccard >= 0.0 && result.jaccard <= 1.0);
            assert!(result.max_containment >= 0.0 && result.max_containment <= 1.0);
            assert!(result.intersect_hashes > 0);
            assert!(result.ksize > 0);
            assert!(result.scaled > 0);

            // Test abundance statistics
            assert!(result.average_abund >= 0.0);
            assert!(result.median_abund >= 0.0);
            assert!(result.std_abund >= 0.0);

            // Test ANI values
            assert!(result.query_containment_ani >= 0.0 && result.query_containment_ani <= 1.0);
            assert!(result.match_containment_ani >= 0.0 && result.match_containment_ani <= 1.0);
            assert!(result.average_containment_ani >= 0.0 && result.average_containment_ani <= 1.0);
            assert!(result.max_containment_ani >= 0.0 && result.max_containment_ani <= 1.0);

            // Test weighted metrics
            assert!(result.n_weighted_found > 0);
            assert!(result.total_weighted_hashes > 0);
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

        let results = searcher.search_multiple(&query_signatures)?;

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

        let query = ProteinSignature::new("test", 10, 5, "hp").unwrap();
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
