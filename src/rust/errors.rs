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

    #[error("Builder error: {0}")]
    BuilderError(String),

    #[error("Sourmash error: {0}")]
    SourmashError(String),

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Anyhow error: {0}")]
    AnyhowError(#[from] anyhow::Error),

    #[error("Validation error: {message}")]
    ValidationError { message: String },

    #[error("Configuration error: {field} - {message}")]
    ConfigurationError { field: String, message: String },

    #[error("Performance warning: {message}")]
    PerformanceWarning { message: String },
}

pub type IndexResult<T> = Result<T, IndexError>;

/// Extension trait for Result to add context
pub trait IndexResultExt<T> {
    /// Add context to an error
    fn with_context<F>(self, f: F) -> IndexResult<T>
    where
        F: FnOnce() -> String;

    /// Add context with a field name
    fn with_field_context(self, field: &str, message: &str) -> IndexResult<T>;
}

impl<T, E> IndexResultExt<T> for Result<T, E>
where
    E: Into<IndexError>,
{
    fn with_context<F>(self, f: F) -> IndexResult<T>
    where
        F: FnOnce() -> String,
    {
        self.map_err(|e| {
            let context = f();
            match e.into() {
                IndexError::ValidationError { message } => {
                    IndexError::ValidationError { message: format!("{}: {}", context, message) }
                }
                other => IndexError::ValidationError { message: format!("{}: {}", context, other) },
            }
        })
    }

    fn with_field_context(self, field: &str, message: &str) -> IndexResult<T> {
        self.map_err(|e| IndexError::ConfigurationError {
            field: field.to_string(),
            message: format!("{}: {}", message, e.into()),
        })
    }
}
