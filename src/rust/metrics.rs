use serde::{Deserialize, Serialize};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{Duration, Instant};

/// Performance metrics for the index
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexMetrics {
    /// Number of signatures processed
    pub signatures_processed: u64,
    /// Number of k-mers processed
    pub kmers_processed: u64,
    /// Total processing time
    pub total_processing_time: Duration,
    /// Average time per signature
    pub avg_time_per_signature: Duration,
    /// Memory usage (approximate)
    pub memory_usage_bytes: u64,
    /// Number of database operations
    pub database_operations: u64,
    /// Number of cache hits
    pub cache_hits: u64,
    /// Number of cache misses
    pub cache_misses: u64,
}

impl Default for IndexMetrics {
    fn default() -> Self {
        Self {
            signatures_processed: 0,
            kmers_processed: 0,
            total_processing_time: Duration::ZERO,
            avg_time_per_signature: Duration::ZERO,
            memory_usage_bytes: 0,
            database_operations: 0,
            cache_hits: 0,
            cache_misses: 0,
        }
    }
}

impl IndexMetrics {
    /// Calculate cache hit ratio
    pub fn cache_hit_ratio(&self) -> f64 {
        let total = self.cache_hits + self.cache_misses;
        if total == 0 {
            0.0
        } else {
            self.cache_hits as f64 / total as f64
        }
    }

    /// Calculate k-mers per second
    pub fn kmers_per_second(&self) -> f64 {
        if self.total_processing_time.as_secs_f64() == 0.0 {
            0.0
        } else {
            self.kmers_processed as f64 / self.total_processing_time.as_secs_f64()
        }
    }

    /// Calculate signatures per second
    pub fn signatures_per_second(&self) -> f64 {
        if self.total_processing_time.as_secs_f64() == 0.0 {
            0.0
        } else {
            self.signatures_processed as f64 / self.total_processing_time.as_secs_f64()
        }
    }
}

/// Thread-safe metrics collector
#[derive(Debug)]
pub struct MetricsCollector {
    signatures_processed: AtomicU64,
    kmers_processed: AtomicU64,
    database_operations: AtomicU64,
    cache_hits: AtomicU64,
    cache_misses: AtomicU64,
    start_time: Instant,
}

impl MetricsCollector {
    /// Create a new metrics collector
    pub fn new() -> Self {
        Self {
            signatures_processed: AtomicU64::new(0),
            kmers_processed: AtomicU64::new(0),
            database_operations: AtomicU64::new(0),
            cache_hits: AtomicU64::new(0),
            cache_misses: AtomicU64::new(0),
            start_time: Instant::now(),
        }
    }

    /// Record a signature being processed
    pub fn record_signature(&self) {
        self.signatures_processed.fetch_add(1, Ordering::Relaxed);
    }

    /// Record k-mers being processed
    pub fn record_kmers(&self, count: u64) {
        self.kmers_processed.fetch_add(count, Ordering::Relaxed);
    }

    /// Record a database operation
    pub fn record_database_operation(&self) {
        self.database_operations.fetch_add(1, Ordering::Relaxed);
    }

    /// Record a cache hit
    pub fn record_cache_hit(&self) {
        self.cache_hits.fetch_add(1, Ordering::Relaxed);
    }

    /// Record a cache miss
    pub fn record_cache_miss(&self) {
        self.cache_misses.fetch_add(1, Ordering::Relaxed);
    }

    /// Get current metrics
    pub fn get_metrics(&self) -> IndexMetrics {
        let total_time = self.start_time.elapsed();
        let signatures_processed = self.signatures_processed.load(Ordering::Relaxed);

        IndexMetrics {
            signatures_processed,
            kmers_processed: self.kmers_processed.load(Ordering::Relaxed),
            total_processing_time: total_time,
            avg_time_per_signature: if signatures_processed > 0 {
                Duration::from_nanos(total_time.as_nanos() as u64 / signatures_processed)
            } else {
                Duration::ZERO
            },
            memory_usage_bytes: 0, // Would need to implement memory tracking
            database_operations: self.database_operations.load(Ordering::Relaxed),
            cache_hits: self.cache_hits.load(Ordering::Relaxed),
            cache_misses: self.cache_misses.load(Ordering::Relaxed),
        }
    }
}

impl Default for MetricsCollector {
    fn default() -> Self {
        Self::new()
    }
}

/// Performance timer for measuring operation durations
pub struct PerformanceTimer {
    start: Instant,
    operation: String,
}

impl PerformanceTimer {
    /// Start timing an operation
    pub fn start(operation: &str) -> Self {
        Self { start: Instant::now(), operation: operation.to_string() }
    }

    /// Finish timing and return the duration
    pub fn finish(self) -> Duration {
        self.start.elapsed()
    }

    /// Finish timing and log the result
    pub fn finish_and_log(self) -> Duration {
        let operation = self.operation.clone();
        let duration = self.finish();
        println!("Operation '{}' took {:?}", operation, duration);
        duration
    }
}

/// Macro for easy performance timing
#[macro_export]
macro_rules! time_operation {
    ($operation:expr, $code:block) => {{
        let timer = $crate::metrics::PerformanceTimer::start($operation);
        let result = $code;
        let _duration = timer.finish();
        result
    }};
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;
    use std::time::Duration;

    #[test]
    fn test_metrics_collector() {
        let collector = MetricsCollector::new();

        collector.record_signature();
        collector.record_kmers(100);
        collector.record_cache_hit();
        collector.record_cache_miss();

        let metrics = collector.get_metrics();
        assert_eq!(metrics.signatures_processed, 1);
        assert_eq!(metrics.kmers_processed, 100);
        assert_eq!(metrics.cache_hits, 1);
        assert_eq!(metrics.cache_misses, 1);
        assert_eq!(metrics.cache_hit_ratio(), 0.5);
    }

    #[test]
    fn test_performance_timer() {
        let timer = PerformanceTimer::start("test_operation");
        thread::sleep(Duration::from_millis(10));
        let duration = timer.finish();

        assert!(duration >= Duration::from_millis(10));
    }
}
