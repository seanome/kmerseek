use crate::errors::{IndexError, IndexResult};
use crate::types::{KmerSize, MolType, Scaled};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Configuration for the ProteomeIndex
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexConfig {
    /// K-mer size
    pub ksize: KmerSize,
    /// Scaled value for MinHash
    pub scaled: Scaled,
    /// Molecular type
    pub moltype: MolType,
    /// Database path
    pub db_path: PathBuf,
    /// Whether to store raw sequences
    pub store_raw_sequences: bool,
    /// Performance settings
    pub performance: PerformanceConfig,
    /// Memory settings
    pub memory: MemoryConfig,
}

/// Performance-related configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceConfig {
    /// Number of threads for parallel processing
    pub num_threads: Option<usize>,
    /// Batch size for processing
    pub batch_size: usize,
    /// Enable performance metrics
    pub enable_metrics: bool,
    /// Progress reporting interval
    pub progress_interval: u32,
}

/// Memory-related configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MemoryConfig {
    /// Maximum memory usage in MB
    pub max_memory_mb: Option<usize>,
    /// Enable memory monitoring
    pub enable_monitoring: bool,
    /// Cache size for k-mers
    pub kmer_cache_size: usize,
}

impl Default for PerformanceConfig {
    fn default() -> Self {
        Self {
            num_threads: None, // Use system default
            batch_size: 1000,
            enable_metrics: true,
            progress_interval: 1000,
        }
    }
}

impl Default for MemoryConfig {
    fn default() -> Self {
        Self {
            max_memory_mb: None, // No limit
            enable_monitoring: false,
            kmer_cache_size: 10000,
        }
    }
}

impl IndexConfig {
    /// Create a new configuration with validation
    pub fn new(
        ksize: u32,
        scaled: u32,
        moltype: &str,
        db_path: PathBuf,
        store_raw_sequences: bool,
    ) -> IndexResult<Self> {
        let ksize = KmerSize::new(ksize).map_err(|e| IndexError::ConfigurationError {
            field: "ksize".to_string(),
            message: e,
        })?;

        let scaled = Scaled::new(scaled).map_err(|e| IndexError::ConfigurationError {
            field: "scaled".to_string(),
            message: e,
        })?;

        let moltype = MolType::new(moltype).map_err(|e| IndexError::ConfigurationError {
            field: "moltype".to_string(),
            message: e,
        })?;

        Ok(Self {
            ksize,
            scaled,
            moltype,
            db_path,
            store_raw_sequences,
            performance: PerformanceConfig::default(),
            memory: MemoryConfig::default(),
        })
    }

    /// Validate the configuration
    pub fn validate(&self) -> IndexResult<()> {
        // Check if database path is writable
        if let Some(parent) = self.db_path.parent() {
            if !parent.exists() {
                return Err(IndexError::ConfigurationError {
                    field: "db_path".to_string(),
                    message: format!("Database directory does not exist: {:?}", parent),
                });
            }
        }

        // Validate performance settings
        if self.performance.batch_size == 0 {
            return Err(IndexError::ConfigurationError {
                field: "batch_size".to_string(),
                message: "Batch size must be greater than 0".to_string(),
            });
        }

        // Validate memory settings
        if let Some(max_memory) = self.memory.max_memory_mb {
            if max_memory == 0 {
                return Err(IndexError::ConfigurationError {
                    field: "max_memory_mb".to_string(),
                    message: "Maximum memory must be greater than 0".to_string(),
                });
            }
        }

        Ok(())
    }

    /// Get the number of threads to use
    pub fn effective_thread_count(&self) -> usize {
        self.performance
            .num_threads
            .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1))
    }

    /// Check if memory monitoring is enabled
    pub fn should_monitor_memory(&self) -> bool {
        self.memory.enable_monitoring || self.memory.max_memory_mb.is_some()
    }
}

/// Configuration builder for more complex setups
pub struct IndexConfigBuilder {
    config: IndexConfig,
}

impl IndexConfigBuilder {
    /// Create a new builder with basic configuration
    pub fn new(ksize: u32, scaled: u32, moltype: &str, db_path: PathBuf) -> IndexResult<Self> {
        let config = IndexConfig::new(ksize, scaled, moltype, db_path, false)?;
        Ok(Self { config })
    }

    /// Set whether to store raw sequences
    pub fn store_raw_sequences(mut self, store: bool) -> Self {
        self.config.store_raw_sequences = store;
        self
    }

    /// Set the number of threads
    pub fn num_threads(mut self, threads: usize) -> Self {
        self.config.performance.num_threads = Some(threads);
        self
    }

    /// Set the batch size
    pub fn batch_size(mut self, size: usize) -> Self {
        self.config.performance.batch_size = size;
        self
    }

    /// Enable or disable metrics
    pub fn enable_metrics(mut self, enable: bool) -> Self {
        self.config.performance.enable_metrics = enable;
        self
    }

    /// Set the progress interval
    pub fn progress_interval(mut self, interval: u32) -> Self {
        self.config.performance.progress_interval = interval;
        self
    }

    /// Set maximum memory usage
    pub fn max_memory_mb(mut self, memory_mb: usize) -> Self {
        self.config.memory.max_memory_mb = Some(memory_mb);
        self
    }

    /// Enable memory monitoring
    pub fn enable_memory_monitoring(mut self, enable: bool) -> Self {
        self.config.memory.enable_monitoring = enable;
        self
    }

    /// Set k-mer cache size
    pub fn kmer_cache_size(mut self, size: usize) -> Self {
        self.config.memory.kmer_cache_size = size;
        self
    }

    /// Build the final configuration
    pub fn build(self) -> IndexResult<IndexConfig> {
        self.config.validate()?;
        Ok(self.config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_config_creation() {
        let config =
            IndexConfig::new(5, 5, "protein", PathBuf::from("/tmp/test.db"), false).unwrap();

        assert_eq!(config.ksize.get(), 5);
        assert_eq!(config.scaled.get(), 5);
        assert_eq!(config.moltype.get(), "protein");
        assert_eq!(config.store_raw_sequences, false);
    }

    #[test]
    fn test_config_validation() {
        // Test invalid ksize
        let config = IndexConfig::new(
            0, // Invalid ksize
            5,
            "protein",
            PathBuf::from("/tmp/test.db"),
            false,
        );
        assert!(config.is_err());

        // Test invalid scaled value
        let config = IndexConfig::new(
            5,
            100, // Too large scaled value
            "protein",
            PathBuf::from("/tmp/test.db"),
            false,
        );
        assert!(config.is_err());

        // Test invalid moltype
        let config = IndexConfig::new(
            5,
            5,
            "invalid", // Invalid moltype
            PathBuf::from("/tmp/test.db"),
            false,
        );
        assert!(config.is_err());
    }

    #[test]
    fn test_config_builder() {
        let config = IndexConfigBuilder::new(5, 5, "protein", PathBuf::from("/tmp/test.db"))
            .unwrap()
            .store_raw_sequences(true)
            .num_threads(4)
            .batch_size(500)
            .build()
            .unwrap();

        assert_eq!(config.store_raw_sequences, true);
        assert_eq!(config.performance.num_threads, Some(4));
        assert_eq!(config.performance.batch_size, 500);
    }
}
