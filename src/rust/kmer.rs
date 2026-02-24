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
