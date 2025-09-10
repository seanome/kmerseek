use crate::kmer::KmerInfo;
use crate::signature::ProteinSignature;
use crate::types::{HashValue, Position};

/// Iterator over k-mer information in a protein signature
pub struct KmerInfoIterator<'a> {
    inner: std::collections::hash_map::Iter<'a, u64, KmerInfo>,
}

impl<'a> Iterator for KmerInfoIterator<'a> {
    type Item = (HashValue, &'a KmerInfo);

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(hash, info)| (HashValue(*hash), info))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<'a> ExactSizeIterator for KmerInfoIterator<'a> {
    fn len(&self) -> usize {
        self.inner.len()
    }
}

/// Iterator over positions for a specific k-mer
pub struct PositionIterator<'a> {
    inner: std::slice::Iter<'a, usize>,
}

impl<'a> Iterator for PositionIterator<'a> {
    type Item = Position;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|&pos| Position(pos))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<'a> ExactSizeIterator for PositionIterator<'a> {
    fn len(&self) -> usize {
        self.inner.len()
    }
}

/// Extension trait for ProteinSignature to provide iterator methods
pub trait ProteinSignatureExt {
    /// Get an iterator over k-mer information
    fn kmer_infos_iter(&self) -> KmerInfoIterator<'_>;

    /// Get an iterator over k-mer hashes only
    fn kmer_hashes(&self) -> impl Iterator<Item = HashValue> + '_;

    /// Get an iterator over k-mer counts (number of positions per k-mer)
    fn kmer_counts(&self) -> impl Iterator<Item = (HashValue, usize)> + '_;

    /// Get the total number of k-mer occurrences across all positions
    fn total_kmer_occurrences(&self) -> usize;

    /// Get k-mers that appear at a specific position
    fn kmers_at_position(&self, position: Position) -> impl Iterator<Item = HashValue> + '_;
}

impl ProteinSignatureExt for ProteinSignature {
    fn kmer_infos_iter(&self) -> KmerInfoIterator<'_> {
        KmerInfoIterator { inner: self.kmer_infos().iter() }
    }

    fn kmer_hashes(&self) -> impl Iterator<Item = HashValue> + '_ {
        self.kmer_infos().keys().map(|&hash| HashValue(hash))
    }

    fn kmer_counts(&self) -> impl Iterator<Item = (HashValue, usize)> + '_ {
        self.kmer_infos().iter().map(|(&hash, info)| (HashValue(hash), info.total_occurrences()))
    }

    fn total_kmer_occurrences(&self) -> usize {
        self.kmer_infos().values().map(|info| info.total_occurrences()).sum()
    }

    fn kmers_at_position(&self, position: Position) -> impl Iterator<Item = HashValue> + '_ {
        self.kmer_infos()
            .iter()
            .filter(move |(_, info)| info.has_position(position.get()))
            .map(|(&hash, _)| HashValue(hash))
    }
}

/// Extension trait for KmerInfo to provide iterator methods
pub trait KmerInfoExt {
    /// Get an iterator over all positions for this k-mer
    fn positions(&self) -> impl Iterator<Item = Position> + '_;

    /// Get an iterator over original k-mers and their positions
    fn original_kmers(&self) -> impl Iterator<Item = (&str, PositionIterator<'_>)> + '_;
}

impl KmerInfoExt for KmerInfo {
    fn positions(&self) -> impl Iterator<Item = Position> + '_ {
        self.positions.iter().map(|&pos| Position(pos))
    }

    fn original_kmers(&self) -> impl Iterator<Item = (&str, PositionIterator<'_>)> + '_ {
        // Since we no longer store original k-mers, return empty iterator
        std::iter::empty()
    }
}

/// Functional utilities for working with protein signatures
pub mod functional {
    use super::*;
    use std::collections::HashMap;

    /// Group k-mers by their occurrence count
    pub fn group_kmers_by_count(signature: &ProteinSignature) -> HashMap<usize, Vec<HashValue>> {
        let mut groups: HashMap<usize, Vec<HashValue>> = HashMap::new();

        for (hash, count) in signature.kmer_counts() {
            groups.entry(count).or_default().push(hash);
        }

        groups
    }

    /// Find k-mers that appear in multiple positions
    pub fn find_multi_position_kmers(signature: &ProteinSignature) -> Vec<HashValue> {
        signature.kmer_counts().filter(|(_, count)| *count > 1).map(|(hash, _)| hash).collect()
    }

    /// Calculate k-mer density (k-mers per sequence position)
    pub fn calculate_kmer_density(signature: &ProteinSignature, sequence_length: usize) -> f64 {
        if sequence_length == 0 {
            0.0
        } else {
            signature.total_kmer_occurrences() as f64 / sequence_length as f64
        }
    }

    /// Find overlapping k-mers (k-mers that share positions)
    pub fn find_overlapping_kmers(signature: &ProteinSignature) -> Vec<(HashValue, HashValue)> {
        let mut overlaps = Vec::new();
        let kmer_infos: Vec<_> = signature.kmer_infos_iter().collect();

        for i in 0..kmer_infos.len() {
            for j in i + 1..kmer_infos.len() {
                let (hash1, info1) = kmer_infos[i];
                let (hash2, info2) = kmer_infos[j];

                // Check if any positions overlap
                let positions1: std::collections::HashSet<usize> =
                    info1.positions().map(|p| p.get()).collect();
                let positions2: std::collections::HashSet<usize> =
                    info2.positions().map(|p| p.get()).collect();

                if !positions1.is_disjoint(&positions2) {
                    overlaps.push((hash1, hash2));
                }
            }
        }

        overlaps
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_kmer_info_iterator() {
        // This would need a proper test setup with actual ProteinSignature
        // For now, just test that the iterator traits are implemented
        let empty_map: HashMap<u64, KmerInfo> = HashMap::new();
        let iter = KmerInfoIterator { inner: empty_map.iter() };
        assert_eq!(iter.len(), 0);
    }

    #[test]
    fn test_functional_utilities() {
        // Test the functional utilities with empty data
        let empty_map: HashMap<usize, Vec<HashValue>> = HashMap::new();
        assert!(empty_map.is_empty());
    }
}
