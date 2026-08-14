use std::ops::Deref;

/// An immutable k-mer string with pre-allocated capacity
///
/// WHY: K-mers have fixed length and are never modified after creation.
/// This struct ensures capacity is always set correctly and provides
/// type safety to distinguish k-mers from regular strings.
#[derive(Debug, Clone)]
pub struct Kmer {
    inner: String,
}

impl Kmer {
    pub fn new(subseq: &str) -> Self {
        let mut s = String::with_capacity(subseq.len());
        s.push_str(subseq);
        Self { inner: s }
    }

    pub fn from_sequence(sequence: &str, start: usize, ksize: usize) -> Self {
        Self::new(&sequence[start..start + ksize])
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

/// True if every byte in a k-mer is the same.
///
/// This is the lowest-complexity case for any alphabet: a raw amino-acid
/// homopolymer such as `"AAAAA"`, or a homopolymer run under a reduced
/// alphabet like HP, such as `"hhhhh"`. Empty input is not a homopolymer.
///
/// Case-insensitive, so it accepts both lowercase and uppercase encodings.
pub fn is_homopolymer_kmer(kmer: &[u8]) -> bool {
    match kmer.split_first() {
        None => false,
        Some((&first, rest)) => rest.iter().all(|b| b.eq_ignore_ascii_case(&first)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_homopolymer_kmer_detects_runs() {
        assert!(is_homopolymer_kmer(b"AAAAA"));
        assert!(is_homopolymer_kmer(b"hhhhh"));
        assert!(is_homopolymer_kmer(b"ppppp"));
        // Case-insensitive: matches the uppercase encoding used for custom HP alphabets.
        assert!(is_homopolymer_kmer(b"HHHHH"));
        // Mixed case of the same residue/class is still a homopolymer run.
        assert!(is_homopolymer_kmer(b"aAaA"));
        // A single-residue k-mer is trivially a homopolymer run.
        assert!(is_homopolymer_kmer(b"A"));
    }

    #[test]
    fn is_homopolymer_kmer_rejects_mixed_kmers() {
        assert!(!is_homopolymer_kmer(b"MKTAY"));
        assert!(!is_homopolymer_kmer(b"hhphh"));
        assert!(!is_homopolymer_kmer(b"ph"));
        assert!(!is_homopolymer_kmer(b""));
    }
}
