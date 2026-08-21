pub mod aminoacid;
pub mod config;
pub mod encoding;
pub mod errors;
pub mod hp_alphabets;
pub mod index;
pub mod io;
pub mod iterators;
pub mod kmer;
pub mod metrics;
pub mod reduced_alphabets;
pub mod search;
pub mod signature;
pub mod significance;
pub mod sketch;
pub mod types;

#[cfg(test)]
mod tests;

// Re-export main types for easier access
pub use aminoacid::AminoAcidAmbiguity;
pub use index::ProteomeIndex;
pub use signature::SEED;
