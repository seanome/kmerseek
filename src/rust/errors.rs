use thiserror::Error;

#[derive(Debug, Error)]
pub enum IndexError {
    #[error("Database error: {0}")]
    Database(#[from] rocksdb::Error),

    #[error("Invalid molecular type: {0}")]
    InvalidMoltype(String),

    #[error("Serialization error: {0}")]
    Serialization(#[from] bincode::Error),

    #[error("Invalid amino acid '{0}' found at position {1}")]
    InvalidAminoAcid(char, usize),

    #[error("Lock acquisition failed: {0}")]
    LockError(String),

    #[error("Invalid k-mer size: {0}")]
    InvalidKsize(u32),

    #[error("No saved state found in database")]
    NoSavedState,

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("UTF-8 error: {0}")]
    Utf8(#[from] std::str::Utf8Error),

    #[error("FASTA parsing error: {0}")]
    FastaParsing(String),
}

pub type IndexResult<T> = Result<T, IndexError>;
