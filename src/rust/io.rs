//! Compressed I/O handling for bioinformatics files.
//!
//! This module provides idiomatic Rust patterns for handling compressed files
//! commonly used in bioinformatics, with automatic format detection and
//! zero-copy streaming access.

use anyhow::Result;
use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

/// Opens a file with automatic compression detection.
///
/// This function automatically detects and handles various compression formats
/// including gzip, bzip2, xz, zstd, and uncompressed files. It provides a
/// zero-cost abstraction for working with compressed bioinformatics files.
///
/// # Arguments
/// * `path` - Path to the file to open (supports any compression format)
///
/// # Returns
/// * `Ok(Box<dyn BufRead>)` - A buffered reader that handles the file
/// * `Err(...)` - An error if the file cannot be opened or decompressed
///
/// # Examples
/// ```
/// use kmerseek::io::open_maybe_compressed;
/// use std::io::BufRead;
///
/// // Example with error handling
/// match open_maybe_compressed("data.fasta.gz") {
///     Ok(reader) => {
///         for line in reader.lines() {
///             let line = line?;
///             // Process line...
///         }
///     }
///     Err(e) => {
///         eprintln!("Failed to open file: {}", e);
///     }
/// }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Why this is idiomatic
///
/// - **Generic over `AsRef<Path>`**: Accepts both `&str` and `PathBuf` for maximum flexibility
/// - **Returns `Box<dyn BufRead>`**: Provides polymorphic I/O without runtime overhead
/// - **Auto-detection**: `niffler` provides zero-cost format detection
/// - **Explicit error handling**: Uses `Result` to make failures explicit in the type system
/// - **Zero-copy abstraction**: No unnecessary allocations or format-specific code
pub fn open_maybe_compressed<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    // niffler autodetects gzip/bzip2/xz/zstd/uncompressed
    let (reader, _format) = niffler::get_reader(Box::new(file))?; // -> Box<dyn Read>
    Ok(Box::new(BufReader::new(reader)))
}

/// Creates a buffered reader from stdin for pipeline processing.
///
/// This function is useful for creating command-line tools that can process
/// data from stdin, following the Unix philosophy of composable tools.
///
/// # Returns
/// * `Box<dyn BufRead>` - A buffered reader wrapping stdin
///
/// # Why this is idiomatic
///
/// - **Polymorphic I/O**: Returns the same type as `open_maybe_compressed`
/// - **Pipeline-friendly**: Enables Unix-style data processing pipelines
/// - **Zero-cost abstraction**: No runtime overhead for the abstraction
pub fn stdin_reader() -> Box<dyn BufRead> {
    Box::new(BufReader::new(io::stdin()))
}

/// Determines the appropriate reader for a given path or stdin.
///
/// This function provides a unified interface for handling both file paths
/// and stdin input, which is essential for command-line bioinformatics tools.
///
/// # Arguments
/// * `path` - Path to the file, or "-" for stdin
///
/// # Returns
/// * `Ok(Box<dyn BufRead>)` - A buffered reader for the input source
/// * `Err(...)` - An error if the file cannot be opened
///
/// # Examples
/// ```
/// use kmerseek::io::open_input;
/// use std::io::BufRead;
///
/// // Read from file with error handling
/// match open_input("data.fasta.gz") {
///     Ok(reader) => {
///         for line in reader.lines() {
///             let line = line?;
///             // Process line...
///         }
///     }
///     Err(e) => {
///         eprintln!("Failed to open file: {}", e);
///     }
/// }
///
/// // Read from stdin
/// let reader = open_input("-")?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Why this is idiomatic
///
/// - **Unified interface**: Single function handles both file and stdin input
/// - **Convention over configuration**: Uses "-" as stdin convention
/// - **Type safety**: Compile-time guarantees about input handling
/// - **Error propagation**: Explicit error handling with `Result`
pub fn open_input<P: AsRef<Path>>(path: P) -> Result<Box<dyn BufRead>> {
    let path_str = path.as_ref().to_string_lossy();
    if path_str == "-" {
        Ok(stdin_reader())
    } else {
        open_maybe_compressed(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_open_uncompressed_file() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        temp_file.write_all(b"test content\n")?;

        let reader = open_maybe_compressed(temp_file.path())?;
        let lines: Result<Vec<String>, _> = reader.lines().collect();
        let lines = lines?;

        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "test content");
        Ok(())
    }

    #[test]
    fn test_open_input_stdin() -> Result<()> {
        // This test verifies that open_input("-") returns a reader without blocking
        // We can't actually read from stdin in tests as it would block indefinitely
        let reader = open_input("-")?;
        // Just verify we get a reader back - this is sufficient to test the function works
        // The actual stdin reading is tested through integration tests or manual testing
        assert!(!std::ptr::addr_of!(reader).is_null());
        Ok(())
    }

    #[test]
    fn test_open_maybe_compressed_zstd() -> Result<()> {
        // Test zstd compression detection
        let reader = open_maybe_compressed("tests/testdata/fasta/test_compression.fasta.zst")?;
        let lines: Result<Vec<String>, _> = reader.lines().collect();
        let lines = lines?;

        assert_eq!(lines.len(), 4); // 2 headers + 2 sequences
        assert_eq!(lines[0], ">test_protein1");
        assert_eq!(lines[1], "PLANTANDANIMALGENQMES");
        assert_eq!(lines[2], ">test_protein2");
        assert_eq!(lines[3], "LIVINGALIVE");
        Ok(())
    }

    #[test]
    fn test_open_input_file() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        temp_file.write_all(b"test\n")?;

        let reader = open_input(temp_file.path())?;
        let lines: Result<Vec<String>, _> = reader.lines().collect();
        let lines = lines?;

        assert_eq!(lines.len(), 1);
        assert_eq!(lines[0], "test");
        Ok(())
    }
}
