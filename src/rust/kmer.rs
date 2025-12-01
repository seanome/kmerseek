use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::ops::Deref;

/// An immutable k-mer string with pre-allocated capacity
///
/// WHY: K-mers have fixed length and are never modified after creation.
/// This struct ensures capacity is always set correctly and provides
/// type safety to distinguish k-mers from regular strings.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Kmer {
    // WHY: Using a private String ensures immutability - the struct can't be
    // modified after construction. The capacity is set to the exact k-mer
    // length to avoid unnecessary allocations.
    inner: String,
}

impl Kmer {
    /// Create a new Kmer from an already-extracted substring
    ///
    /// WHY: This constructor is for when you already have the k-mer substring.
    /// It pre-allocates the String with the exact capacity needed.
    pub fn new(subseq: &str) -> Self {
        let mut s = String::with_capacity(subseq.len());
        s.push_str(subseq);
        Self { inner: s }
    }

    /// Create a new Kmer by extracting from a sequence
    ///
    /// WHY: This constructor extracts the k-mer substring from a larger sequence
    /// and pre-allocates with the exact capacity. This is the common case when
    /// iterating through a sequence to extract k-mers.
    pub fn from_sequence(sequence: &str, start: usize, ksize: usize) -> Self {
        let end = start + ksize;
        let subseq = &sequence[start..end];
        Self::new(subseq)
    }
}

impl Deref for Kmer {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl AsRef<str> for Kmer {
    fn as_ref(&self) -> &str {
        &self.inner
    }
}

impl Hash for Kmer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.inner.hash(state);
    }
}

impl PartialEq for Kmer {
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl Eq for Kmer {}

impl PartialEq<str> for Kmer {
    fn eq(&self, other: &str) -> bool {
        self.inner == other
    }
}

impl PartialEq<&str> for Kmer {
    fn eq(&self, other: &&str) -> bool {
        self.inner == *other
    }
}

/// Represents information about a k-mer occurrence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerInfo {
    // Sourmash branchwater uses u32 for ksize so we will, too
    pub ksize: u32,
    pub hashval: u64,
    pub encoded_kmer: String,
    pub original_kmer_to_position: HashMap<String, Vec<usize>>,
}

impl KmerInfo {
    /// Creates a new KmerInfo with pre-allocated capacity
    pub fn new(hashval: u64, ksize: u32) -> Self {
        Self {
            ksize,
            hashval,
            encoded_kmer: String::with_capacity(ksize.try_into().unwrap()),
            original_kmer_to_position: HashMap::new(),
        }
    }

    /// Adds a k-mer position
    ///
    /// WHY: We use String::with_capacity(ksize) because k-mer lengths are fixed
    /// and will never change. This avoids unnecessary reallocations when creating
    /// the String key for the HashMap.
    pub fn add_position(&mut self, kmer: &str, position: usize) {
        // WHY: Pre-allocate String with exact k-mer size since k-mer length is fixed.
        // This is more efficient than to_string() which may allocate more than needed.
        let kmer_string = {
            let mut s = String::with_capacity(self.ksize as usize);
            s.push_str(kmer);
            s
        };
        self.original_kmer_to_position
            .entry(kmer_string)
            .or_insert_with(|| Vec::with_capacity(1))
            .push(position);
    }

    /// Get the number of unique original k-mers
    pub fn unique_kmer_count(&self) -> usize {
        self.original_kmer_to_position.len()
    }

    /// Get the number of k-mer occurrences
    pub fn total_occurrences(&self) -> usize {
        self.original_kmer_to_position.values().map(|positions| positions.len()).sum()
    }

    /// Check if this k-mer appears at a specific position
    pub fn has_position(&self, position: usize) -> bool {
        self.original_kmer_to_position.values().any(|positions| positions.contains(&position))
    }

    /// Get all positions where this k-mer occurs
    pub fn get_positions(&self) -> Vec<usize> {
        self.original_kmer_to_position.values().flatten().cloned().collect()
    }
}
