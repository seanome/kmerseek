use crate::sketch::ProteinSketch;
use crate::types::{HashValue, Position};

/// Iterator over (hashval, positions) entries in a protein sketch
pub struct KmerPositionsIterator<'a> {
    inner: std::collections::hash_map::Iter<'a, u64, Vec<usize>>,
}

impl<'a> Iterator for KmerPositionsIterator<'a> {
    type Item = (HashValue, &'a Vec<usize>);

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(hash, positions)| (HashValue(*hash), positions))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<'a> ExactSizeIterator for KmerPositionsIterator<'a> {
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

/// Extension trait for ProteinSketch to provide iterator methods
pub trait ProteinSketchExt {
    /// Get an iterator over (hashval, positions) entries
    fn kmer_positions_iter(&self) -> KmerPositionsIterator<'_>;

    /// Get an iterator over k-mer hashes only
    fn kmer_hashes(&self) -> impl Iterator<Item = HashValue> + '_;

    /// Get an iterator over k-mer counts (number of positions per k-mer)
    fn kmer_counts(&self) -> impl Iterator<Item = (HashValue, usize)> + '_;

    /// Get the total number of k-mer position entries across all hashes
    fn total_kmer_occurrences(&self) -> usize;

    /// Get hashes that have a position entry at the given position
    fn kmers_at_position(&self, position: Position) -> impl Iterator<Item = HashValue> + '_;
}

impl ProteinSketchExt for ProteinSketch {
    fn kmer_positions_iter(&self) -> KmerPositionsIterator<'_> {
        KmerPositionsIterator { inner: self.kmer_positions().iter() }
    }

    fn kmer_hashes(&self) -> impl Iterator<Item = HashValue> + '_ {
        self.kmer_positions().keys().map(|&hash| HashValue(hash))
    }

    fn kmer_counts(&self) -> impl Iterator<Item = (HashValue, usize)> + '_ {
        self.kmer_positions().iter().map(|(&hash, positions)| (HashValue(hash), positions.len()))
    }

    fn total_kmer_occurrences(&self) -> usize {
        self.kmer_positions().values().map(|positions| positions.len()).sum()
    }

    fn kmers_at_position(&self, position: Position) -> impl Iterator<Item = HashValue> + '_ {
        self.kmer_positions()
            .iter()
            .filter(move |(_, positions)| positions.contains(&position.get()))
            .map(|(&hash, _)| HashValue(hash))
    }
}

/// Functional utilities for working with protein signatures
pub mod functional {
    use super::*;
    use std::collections::HashMap;

    /// Group k-mers by their occurrence count
    pub fn group_kmers_by_count(signature: &ProteinSketch) -> HashMap<usize, Vec<HashValue>> {
        let mut groups: HashMap<usize, Vec<HashValue>> = HashMap::new();
        for (hash, count) in signature.kmer_counts() {
            groups.entry(count).or_default().push(hash);
        }
        groups
    }

    /// Find k-mers that appear at more than one position
    pub fn find_multi_position_kmers(signature: &ProteinSketch) -> Vec<HashValue> {
        signature.kmer_counts().filter(|(_, count)| *count > 1).map(|(hash, _)| hash).collect()
    }

    /// Calculate k-mer density (positions per sequence length)
    pub fn calculate_kmer_density(signature: &ProteinSketch, sequence_length: usize) -> f64 {
        if sequence_length == 0 {
            0.0
        } else {
            signature.total_kmer_occurrences() as f64 / sequence_length as f64
        }
    }

    /// Find hashes whose position sets overlap (share at least one position)
    pub fn find_overlapping_kmers(signature: &ProteinSketch) -> Vec<(HashValue, HashValue)> {
        let mut overlaps = Vec::new();
        let entries: Vec<_> = signature.kmer_positions_iter().collect();

        for i in 0..entries.len() {
            for j in i + 1..entries.len() {
                let (hash1, positions1) = entries[i];
                let (hash2, positions2) = entries[j];

                let set1: std::collections::HashSet<usize> = positions1.iter().copied().collect();
                let set2: std::collections::HashSet<usize> = positions2.iter().copied().collect();

                if !set1.is_disjoint(&set2) {
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
    use crate::sketch::ProteinSketch;
    use crate::types::Position;
    use std::collections::HashMap;

    // 40-residue sequence: the 10-residue motif "MKTAYIAKQR" repeated 4 times.
    // At protein k=5, scaled=1 it yields 36 k-mer occurrences over 10 distinct hashes,
    // each hash recurring across the 4 copies.
    const REPEAT_SEQ: &str = "MKTAYIAKQRMKTAYIAKQRMKTAYIAKQRMKTAYIAKQR";
    const REP_DISTINCT: usize = 10;
    const REP_TOTAL_OCC: usize = 36;
    // Hash of the first k-mer "MKTAY" (murmur64, seed 42, protein encoding).
    const REP_HASH_AT_0: u64 = 8299134706416968033;

    fn sketch() -> ProteinSketch {
        ProteinSketch::from_protein_sequence("p", REPEAT_SEQ, 5, 1, "protein").unwrap()
    }

    #[test]
    fn test_kmer_positions_iterator() {
        let empty_map: HashMap<u64, Vec<usize>> = HashMap::new();
        let iter = KmerPositionsIterator { inner: empty_map.iter() };
        assert_eq!(iter.len(), 0);
    }

    #[test]
    fn test_functional_utilities() {
        let empty_map: HashMap<usize, Vec<HashValue>> = HashMap::new();
        assert!(empty_map.is_empty());
    }

    #[test]
    fn test_protein_sketch_ext_iterators() {
        let s = sketch();
        assert_eq!(s.kmer_positions().len(), REP_DISTINCT);

        assert_eq!(s.kmer_positions_iter().len(), REP_DISTINCT);
        assert_eq!(s.kmer_hashes().count(), REP_DISTINCT);
        assert_eq!(s.kmer_counts().count(), REP_DISTINCT);
        assert_eq!(s.total_kmer_occurrences(), REP_TOTAL_OCC);

        // Exactly the "MKTAY" hash starts at position 0.
        let at_zero: Vec<HashValue> = s.kmers_at_position(Position(0)).collect();
        assert_eq!(at_zero, vec![HashValue(REP_HASH_AT_0)]);
    }

    #[test]
    fn test_functional_grouping_and_density() {
        use functional::*;
        let s = sketch();

        let groups = group_kmers_by_count(&s);
        let grouped_hashes: usize = groups.values().map(|v| v.len()).sum();
        assert_eq!(grouped_hashes, REP_DISTINCT);

        // Every distinct k-mer recurs across the 4 motif copies, so all 10 are multi-position.
        assert_eq!(find_multi_position_kmers(&s).len(), REP_DISTINCT);

        // Density = 36 occurrences / 40 residues = 0.9.
        assert_eq!(calculate_kmer_density(&s, REPEAT_SEQ.len()), 0.9);
        assert_eq!(calculate_kmer_density(&s, 0), 0.0);

        // Distinct k-mers never share a start position, so no overlaps are found,
        // but the call still exercises the overlap-detection loop.
        assert_eq!(find_overlapping_kmers(&s).len(), 0);
    }
}
