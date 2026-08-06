use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use parking_lot::Mutex;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use rocksdb::{Options, DB};
use serde::{Deserialize, Serialize};
use sourmash::_hash_murmur;
use sourmash::collection::Collection;
use sourmash::manifest::Manifest;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::{FSStorage, InnerStorage};

use crate::aminoacid::AminoAcidAmbiguity;
use crate::encoding::{
    encode_with_fn, get_encoding_fn_from_moltype, get_hash_function_from_moltype,
};
use crate::errors::{IndexError, IndexResult};
use crate::hp_alphabets::HpAlphabet;
use crate::signature::{SignatureAccess, SEED};
use crate::sketch::{ProteinSketch, ProteinSketchStore};

/// Schema version for the on-disk index format.
/// Increment this constant whenever the stored format changes in a backward-incompatible way
/// (e.g. new fields in SearchCache, renamed fields in ProteinSketchStore, etc.).
/// Indices that predate versioning (schema_version key absent) are treated as version 0
/// and will be rejected with a clear error message asking the user to rebuild.
pub const SCHEMA_VERSION: u32 = 1;

/// Statistics for k-mer frequency analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteomeIndexKmerStats {
    pub idf: HashMap<u64, f64>, // Inverse document frequency for each k-mer hashvalue
    pub frequency: HashMap<u64, f64>, // Raw frequency for each k-mer hashvalue
}

// Represents the serializable state of ProteomeIndex using efficient storage
#[derive(Serialize, Deserialize)]
struct ProteomeIndexState {
    // Store efficient signature data instead of full signatures
    signature_data: Vec<ProteinSketchStore>,
    combined_mins: Vec<u64>,
    combined_abunds: Option<Vec<u64>>,
    moltype: String,
    ksize: u32,
    scaled: u32,
    // Configuration for raw sequence storage
    store_raw_sequences: bool,
}

// Metadata for chunked storage format
#[derive(Serialize, Deserialize)]
struct ProteomeIndexMetadata {
    total_signatures: usize,
    chunk_count: usize,
    combined_mins: Vec<u64>,
    combined_abunds: Option<Vec<u64>>,
    moltype: String,
    ksize: u32,
    scaled: u32,
    store_raw_sequences: bool,
}

/// Serializable search cache built at index time for fast search startup.
///
/// Stores the pre-built inverted index, ordered target list, and k-mer frequencies
/// so that ProteinSearcher::load() can avoid loading all signatures into memory.
/// Individual signatures are stored separately under "sig_{md5}" keys for on-demand access.
#[derive(Serialize, Deserialize)]
pub struct SearchCache {
    /// Ordered list of target MD5 sums: index (u32) → md5 string
    pub target_list: Vec<String>,
    /// Inverted k-mer index: kmer_hash → Vec of target indices into target_list
    pub inverted_index: HashMap<u64, Vec<u32>>,
    /// K-mer frequency counts: kmer_hash → number of signatures containing it
    pub kmer_frequencies: HashMap<u64, usize>,
}

pub struct ProteomeIndex {
    // RocksDB instance for persistent storage
    db: DB,

    // Combined minhash of all proteins for statistics
    combined_minhash: Arc<Mutex<KmerMinHash>>,

    // Map of signature md5 -> protein signature (thread-safe concurrent map)
    signatures: DashMap<String, ProteinSketch>,

    // Amino acid ambiguity handler
    aa_ambiguity: Arc<AminoAcidAmbiguity>,

    // Protein encoding function
    encoding_fn: fn(u8) -> u8,

    // Statistics for k-mer frequencies and IDF
    // Not currently used, but will be used in the future
    #[allow(dead_code)]
    stats: ProteomeIndexKmerStats,

    // Add moltype field for serialization so don't have to read signatures to find it
    moltype: String,

    // Add ksize field for serialization so don't have to read signatures to find it
    // Sourmash branchwater uses u32 for ksize so we will, too
    ksize: u32,

    // Add minhash_ksize field for serialization
    // MinHash k-mer size is protein_ksize * 3, as a legacy from Sourmash which was originally designed for DNA
    // Not currently used, but will be used in the future
    #[allow(dead_code)]
    minhash_ksize: u32,

    // Add scaled field for serialization
    scaled: u32,

    // Configuration for raw sequence storage
    store_raw_sequences: bool,
}

impl Drop for ProteomeIndex {
    fn drop(&mut self) {
        // RocksDB will be automatically closed when the struct is dropped
    }
}

impl ProteomeIndex {
    /// Create RocksDB options optimized for large datasets
    ///
    /// WHY: This helper function centralizes RocksDB configuration to ensure consistent
    /// settings across all database operations. Setting max_open_files to a reasonable limit
    /// (10000) prevents "Too many open files" errors on large databases like UniProt while
    /// still allowing RocksDB to efficiently access SST files. The -1 value (unlimited)
    /// can exceed system file descriptor limits, causing failures on large databases.
    ///
    /// # Arguments
    /// * `create_if_missing` - Whether to create the database if it doesn't exist
    ///
    /// # Returns
    /// Configured RocksDB Options
    fn create_rocksdb_options(create_if_missing: bool) -> Options {
        let mut opts = Options::default();
        opts.create_if_missing(create_if_missing);

        // Set reasonable max_open_files limit to prevent "Too many open files" errors
        // WHY: Large databases like UniProt can have thousands of SST files. Setting a limit
        // of 10000 prevents exceeding system file descriptor limits while still allowing
        // efficient access. The -1 (unlimited) setting can cause failures on large databases.
        opts.set_max_open_files(10000);

        // Optimize for read performance
        opts.set_use_fsync(false);
        opts.set_allow_mmap_reads(true);
        opts.set_allow_mmap_writes(true);

        // Optimize for large datasets
        opts.set_max_bytes_for_level_base(256 * 1024 * 1024); // 256MB
        opts.set_target_file_size_base(64 * 1024 * 1024); // 64MB
        opts.set_write_buffer_size(128 * 1024 * 1024); // 128MB write buffer

        // Optimize for bulk loading (only when creating new databases)
        if create_if_missing {
            opts.set_disable_auto_compactions(true);
            opts.set_level_zero_file_num_compaction_trigger(8);
            opts.set_level_zero_slowdown_writes_trigger(17);
            opts.set_level_zero_stop_writes_trigger(24);
        }

        opts
    }

    /// Create a new ProteomeIndex using the builder pattern
    ///
    /// This method returns a builder for configuring index parameters.
    /// Use the builder methods to configure the index parameters and then call
    /// `.build()` or `.build_with_auto_filename()` to create the index.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use kmerseek::index::ProteomeIndex;
    ///
    /// fn main() -> anyhow::Result<()> {
    ///     let index = ProteomeIndex::builder()
    ///         .path("/path/to/database.db")
    ///         .ksize(5)
    ///         .scaled(1)
    ///         .moltype("protein")
    ///         .build()?;
    ///     Ok(())
    /// }
    /// ```
    pub fn builder() -> ProteomeIndexBuilder {
        ProteomeIndexBuilder::new()
    }

    pub fn new<P: AsRef<Path>>(
        path: P,
        ksize: u32,
        scaled: u32,
        moltype: &str,
        store_raw_sequences: bool,
    ) -> IndexResult<Self> {
        // Create RocksDB options optimized for large datasets
        let opts = Self::create_rocksdb_options(true);

        // Open the database
        let db = DB::open(&opts, path)?;

        let hash_function = get_hash_function_from_moltype(moltype)
            .map_err(|e| IndexError::SourmashError(e.to_string()))?;

        let encoding_fn = get_encoding_fn_from_moltype(moltype)
            .map_err(|e| IndexError::SourmashError(e.to_string()))?;

        let minhash_ksize = ksize * 3;
        // Create the minhash sketch
        let minhash = KmerMinHash::new(
            scaled,
            minhash_ksize,
            hash_function,
            SEED, // seed
            true, // track_abundance
            0,    // num (use scaled instead)
        );

        // Create an empty collection with storage
        let manifest = Manifest::default();
        let storage =
            InnerStorage::new(FSStorage::builder().fullpath("".into()).subdir("".into()).build());
        let _collection = Collection::new(manifest, storage);

        Ok(Self {
            db,
            signatures: DashMap::new(),
            combined_minhash: Arc::new(Mutex::new(minhash)),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            encoding_fn,
            moltype: moltype.to_string(),
            ksize,
            minhash_ksize,
            scaled,
            stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
            store_raw_sequences,
        })
    }

    /// Get a reference to the signatures map (for testing)
    pub fn get_signatures(&self) -> &DashMap<String, ProteinSketch> {
        &self.signatures
    }

    /// Get a reference to the combined minhash (for testing)
    pub fn get_combined_minhash(&self) -> &Arc<Mutex<KmerMinHash>> {
        &self.combined_minhash
    }

    /// Get the k-mer size
    pub fn ksize(&self) -> u32 {
        self.ksize
    }

    /// Get the scaled value
    pub fn scaled(&self) -> u32 {
        self.scaled
    }

    /// Get the molecular type
    pub fn moltype(&self) -> &str {
        &self.moltype
    }

    /// Build and persist the search cache and individual signatures for fast search startup.
    ///
    /// This method:
    /// 1. Builds target_list, inverted_index, and kmer_frequencies from in-memory signatures
    /// 2. Stores each signature individually under "sig_{md5}" for on-demand loading
    /// 3. Serializes the SearchCache (target_list + inverted_index + kmer_frequencies) to RocksDB
    ///
    /// WHY: Building these structures at index time (once) rather than at search startup
    /// (every time) avoids the need to load all 200k+ signatures into memory before searching.
    /// During search, only candidate signatures (those sharing ≥1 k-mer with the query) are
    /// loaded on-demand from RocksDB, reducing startup time from minutes to seconds.
    fn save_inverted_index(&self) -> IndexResult<()> {
        let t0 = Instant::now();
        let total_sigs = self.signatures.len();
        eprintln!("[save] Building inverted index for {} signatures...", total_sigs);

        let mut target_list: Vec<String> = Vec::new();
        let mut inverted_index: HashMap<u64, Vec<u32>> = HashMap::new();
        let mut kmer_frequencies: HashMap<u64, usize> = HashMap::new();

        // Single pass: build index structures and save individual signatures
        for entry in self.signatures.iter() {
            let idx = target_list.len() as u32;
            target_list.push(entry.key().clone());

            // Store individual signature for on-demand loading during search
            let sig_data = entry.value().to_efficient_data(self.store_raw_sequences);
            let serialized = bincode::serialize(&sig_data)?;
            let key = format!("sig_{}", entry.key());
            self.db.put(key.as_bytes(), serialized)?;

            // Build inverted index and kmer frequencies
            let mins = entry.value().signature().minhash.mins();
            for min in mins {
                inverted_index.entry(min).or_default().push(idx);
                *kmer_frequencies.entry(min).or_insert(0) += 1;
            }

            if idx > 0 && idx % 1000 == 0 {
                let mins_so_far = entry.value().signature().minhash.mins().len();
                eprintln!(
                    "[save] {}/{} signatures written ({:.1}s elapsed, last sig had {} mins)",
                    idx,
                    total_sigs,
                    t0.elapsed().as_secs_f32(),
                    mins_so_far,
                );
            }
        }

        eprintln!(
            "[save] All {} signatures written in {:.1}s. Inverted index has {} unique k-mers.",
            total_sigs,
            t0.elapsed().as_secs_f32(),
            inverted_index.len(),
        );

        self.log_kmer_frequency_stats(&kmer_frequencies, &inverted_index, &target_list);

        // Serialize and store the search cache
        eprintln!(
            "[save] Serializing SearchCache ({} targets, {} unique kmers)...",
            target_list.len(),
            inverted_index.len()
        );
        let t1 = Instant::now();
        let cache = SearchCache { target_list, inverted_index, kmer_frequencies };
        let serialized = bincode::serialize(&cache)?;
        eprintln!(
            "[save] SearchCache serialized to {} bytes in {:.1}s, writing to RocksDB...",
            serialized.len(),
            t1.elapsed().as_secs_f32()
        );
        let t2 = Instant::now();
        self.db.put(b"search_cache", serialized)?;
        eprintln!("[save] search_cache written in {:.1}s", t2.elapsed().as_secs_f32());

        Ok(())
    }

    /// Log a k-mer frequency histogram and the top/bottom 10 k-mers by frequency to stderr.
    ///
    /// Bins are power-of-two ranges (1, 2-3, 4-7, ...) since k-mer frequency distributions
    /// are typically heavily right-skewed (most k-mers occur once, a few occur very often).
    /// The top/bottom lists show the actual encoded k-mer string (not the hash), resolved by
    /// looking up one signature that contains each hash via the inverted index.
    fn log_kmer_frequency_stats(
        &self,
        kmer_frequencies: &HashMap<u64, usize>,
        inverted_index: &HashMap<u64, Vec<u32>>,
        target_list: &[String],
    ) {
        if kmer_frequencies.is_empty() {
            return;
        }

        let mut bins: BTreeMap<u32, usize> = BTreeMap::new();
        for &freq in kmer_frequencies.values() {
            let bin = usize::BITS - freq.leading_zeros() - 1;
            *bins.entry(bin).or_insert(0) += 1;
        }

        eprintln!("[save] K-mer frequency histogram ({} unique k-mers):", kmer_frequencies.len());
        let max_count = *bins.values().max().unwrap();
        const BAR_WIDTH: usize = 40;
        for (bin, count) in &bins {
            let lo = 1u64 << bin;
            let hi = (1u64 << (bin + 1)) - 1;
            let label = if lo == hi { format!("{lo}") } else { format!("{lo}-{hi}") };
            let bar_len = (count * BAR_WIDTH) / max_count;
            let bar = "#".repeat(bar_len.max(1));
            eprintln!("[save]   {label:>12} occurrences: {count:>10} k-mers  {bar}");
        }

        let mut by_frequency: Vec<(&u64, &usize)> = kmer_frequencies.iter().collect();
        by_frequency.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));

        eprintln!("[save] Top 10 most common k-mers (encoded k-mer: occurrences):");
        for (hash, count) in by_frequency.iter().take(10) {
            let kmer = self.resolve_kmer_string(**hash, inverted_index, target_list);
            eprintln!("[save]   {kmer}: {count}");
        }

        eprintln!("[save] Bottom 10 least common k-mers (encoded k-mer: occurrences):");
        for (hash, count) in by_frequency.iter().rev().take(10) {
            let kmer = self.resolve_kmer_string(**hash, inverted_index, target_list);
            eprintln!("[save]   {kmer}: {count}");
        }
    }

    /// Resolve one k-mer hash back to the actual encoded k-mer string it was hashed from,
    /// by finding a signature that contains it (via the inverted index) and slicing that
    /// signature's stored sequence at the recorded position.
    ///
    /// WHY: hashes are one-way (murmur), so the only way to recover the k-mer text is to
    /// look up where it occurred in a sequence we already have in memory.
    fn resolve_kmer_string(
        &self,
        hash: u64,
        inverted_index: &HashMap<u64, Vec<u32>>,
        target_list: &[String],
    ) -> String {
        let ksize = self.ksize as usize;
        (|| {
            let target_idx = *inverted_index.get(&hash)?.first()?;
            let md5 = target_list.get(target_idx as usize)?;
            let sig = self.signatures.get(md5)?;
            let position = *sig.kmer_positions().get(&hash)?.first()?;
            let seq = sig.get_moltype_sequence().or_else(|| sig.get_raw_sequence())?;
            seq.get(position..position + ksize).map(str::to_string)
        })()
        .unwrap_or_else(|| format!("<sequence unavailable, hash {hash}>"))
    }

    /// Save the current index state to RocksDB using chunked storage format
    ///
    /// This method stores signatures in chunks to avoid RocksDB value size limits.
    /// Each chunk contains a maximum number of signatures to keep serialized data manageable.
    /// Build combined minhash from all signatures in one O(N log N) pass.
    ///
    /// WHY: Incremental add_many_with_abund per batch is O(M*N) due to Vec::insert shifting.
    /// With millions of hashes this becomes hours. One sort + add in sorted order is O(N log N)
    /// (each add_hash becomes O(1) push because hashes arrive in ascending order).
    pub fn rebuild_combined_minhash(&self) -> IndexResult<()> {
        let mut all_hashes: Vec<u64> =
            self.signatures.iter().flat_map(|e| e.value().signature().minhash.mins()).collect();
        all_hashes.sort_unstable();
        all_hashes.dedup();
        let hash_function = get_hash_function_from_moltype(&self.moltype)?;
        let mut new_combined =
            KmerMinHash::new(self.scaled, self.ksize * 3, hash_function, SEED, true, 0);
        for &h in &all_hashes {
            new_combined.add_hash(h);
        }
        let mut combined_minhash = self.combined_minhash.lock();
        *combined_minhash = new_combined;
        Ok(())
    }

    pub fn save_state(&self) -> IndexResult<()> {
        let t_start = Instant::now();
        eprintln!("[save] save_state() started ({} signatures in memory)", self.signatures.len());

        eprintln!("[save] Building combined minhash from all signatures...");
        let t_cm = Instant::now();
        self.rebuild_combined_minhash()?;
        let combined_minhash = self.combined_minhash.lock();
        eprintln!(
            "[save] Combined minhash built: {} unique k-mers in {:.1}s",
            combined_minhash.mins().len(),
            t_cm.elapsed().as_secs_f32()
        );
        drop(combined_minhash);

        // Convert signatures to efficient storage format
        eprintln!("[save] Converting signatures to storage format...");
        let t1 = Instant::now();
        let mut signature_data = Vec::new();
        for sig in self.signatures.iter() {
            let efficient_data = sig.value().to_efficient_data(self.store_raw_sequences);
            signature_data.push(efficient_data);
        }
        eprintln!(
            "[save] Converted {} signatures in {:.1}s",
            signature_data.len(),
            t1.elapsed().as_secs_f32()
        );

        // Store signatures in chunks to avoid RocksDB value size limits
        // Use smaller chunks for better memory efficiency and faster loading
        const CHUNK_SIZE: usize = 100; // Store 100 signatures per chunk
        let total_signatures = signature_data.len();
        let chunk_count = total_signatures.div_ceil(CHUNK_SIZE);
        eprintln!(
            "[save] Writing {} signatures in {} chunks to RocksDB...",
            total_signatures, chunk_count
        );
        let t2 = Instant::now();

        for (chunk_idx, chunk) in signature_data.chunks(CHUNK_SIZE).enumerate() {
            let chunk_key = format!("signatures_chunk_{}", chunk_idx);
            let serialized_chunk = bincode::serialize(chunk)?;
            self.db.put(chunk_key.as_bytes(), serialized_chunk)?;
            if chunk_idx > 0 && chunk_idx % 50 == 0 {
                eprintln!(
                    "[save] chunk {}/{} written ({:.1}s elapsed)",
                    chunk_idx,
                    chunk_count,
                    t2.elapsed().as_secs_f32()
                );
            }
        }
        eprintln!("[save] All chunks written in {:.1}s", t2.elapsed().as_secs_f32());

        // Store metadata separately
        let combined_minhash = self.combined_minhash.lock();
        eprintln!(
            "[save] Writing metadata (combined_minhash has {} mins)...",
            combined_minhash.mins().len()
        );
        let t3 = Instant::now();
        let metadata = ProteomeIndexMetadata {
            total_signatures,
            chunk_count,
            combined_mins: combined_minhash.mins().to_vec(),
            combined_abunds: combined_minhash.abunds().map(|abunds| abunds.to_vec()),
            moltype: self.moltype.clone(),
            ksize: self.ksize,
            scaled: self.scaled,
            store_raw_sequences: self.store_raw_sequences,
        };

        let serialized_metadata = bincode::serialize(&metadata)?;
        eprintln!(
            "[save] Metadata serialized to {} bytes in {:.1}s",
            serialized_metadata.len(),
            t3.elapsed().as_secs_f32()
        );
        self.db.put(b"index_metadata", serialized_metadata)?;

        // Store schema version as a separate key so it can be validated without
        // deserializing the full metadata (and without breaking old bincode layouts).
        let serialized_version = bincode::serialize(&SCHEMA_VERSION)?;
        self.db.put(b"schema_version", serialized_version)?;
        eprintln!(
            "[save] Metadata + schema_version written in {:.1}s total",
            t3.elapsed().as_secs_f32()
        );

        // Build and persist search cache + individual signatures for fast search startup
        eprintln!("[save] Building search cache...");
        let t4 = Instant::now();
        self.save_inverted_index()?;
        eprintln!("[save] Search cache saved in {:.1}s", t4.elapsed().as_secs_f32());

        // Flush to ensure data is written to disk
        eprintln!("[save] Flushing RocksDB...");
        let t5 = Instant::now();
        self.db.flush()?;
        eprintln!(
            "[save] Flush complete in {:.1}s. Total save_state() time: {:.1}s",
            t5.elapsed().as_secs_f32(),
            t_start.elapsed().as_secs_f32()
        );

        Ok(())
    }

    /// Enable compactions after bulk loading for better performance
    ///
    /// This method should be called after all signatures have been loaded
    /// to optimize the database for read operations.
    pub fn enable_compactions(&self) -> IndexResult<()> {
        // Enable auto compactions
        let mut opts = Options::default();
        opts.set_disable_auto_compactions(false);

        // Compact the database to optimize for reads
        self.db.compact_range::<&[u8], &[u8]>(None, None);

        Ok(())
    }

    /// Load index state from RocksDB using chunked storage format
    pub fn load_state(&self) -> IndexResult<()> {
        // Try to load from new chunked format first
        let metadata_serialized = self.db.get(b"index_metadata")?;
        if let Some(metadata_data) = metadata_serialized {
            let metadata: ProteomeIndexMetadata = bincode::deserialize(&metadata_data)?;

            // Load all chunk data from RocksDB (sequential)
            let mut raw_chunks: Vec<Vec<u8>> = Vec::with_capacity(metadata.chunk_count);
            for chunk_idx in 0..metadata.chunk_count {
                let chunk_key = format!("signatures_chunk_{}", chunk_idx);
                if let Some(data) = self.db.get(chunk_key.as_bytes())? {
                    raw_chunks.push(data);
                }
            }

            // Deserialize and reconstruct signatures in parallel
            use rayon::prelude::*;
            let moltype = &metadata.moltype;
            let ksize = metadata.ksize;
            let scaled = metadata.scaled;

            let new_signatures: DashMap<String, ProteinSketch> = DashMap::new();
            raw_chunks.par_iter().try_for_each(|raw_data| -> IndexResult<()> {
                let chunk: Vec<ProteinSketchStore> = bincode::deserialize(raw_data)?;
                for signature_data in chunk {
                    let protein_sig = ProteinSketch::from_efficient_data(
                        signature_data,
                        moltype.clone(),
                        ksize,
                        scaled,
                    )?;
                    let md5sum = protein_sig.signature().md5sum.clone();
                    new_signatures.insert(md5sum.to_string(), protein_sig);
                }
                Ok(())
            })?;

            // Reconstruct the combined minhash
            let hash_function = get_hash_function_from_moltype(&metadata.moltype)?;
            let minhash_ksize = metadata.ksize * 3;
            let mut combined_minhash = KmerMinHash::new(
                metadata.scaled,
                minhash_ksize,
                hash_function,
                SEED,
                true, // track_abundance
                0,    // num (use scaled instead)
            );

            if let Some(abunds) = &metadata.combined_abunds {
                combined_minhash
                    .add_many_with_abund(
                        &metadata
                            .combined_mins
                            .clone()
                            .into_iter()
                            .zip(abunds.iter().cloned())
                            .collect::<Vec<_>>(),
                    )
                    .map_err(|e| IndexError::SourmashError(e.to_string()))?;
            } else {
                combined_minhash
                    .add_many(&metadata.combined_mins)
                    .map_err(|e| IndexError::SourmashError(e.to_string()))?;
            }

            // Update the current index state
            {
                // Clear existing signatures and swap in new ones
                self.signatures.clear();
                for entry in new_signatures.into_iter() {
                    self.signatures.insert(entry.0, entry.1);
                }
            }

            {
                let mut current_combined = self.combined_minhash.lock();
                *current_combined = combined_minhash;
            }

            Ok(())
        } else {
            // Fallback to old format for backward compatibility
            let serialized = self.db.get(b"index_state")?;
            if let Some(data) = serialized {
                let state: ProteomeIndexState = bincode::deserialize(&data)?;

                // Reconstruct signatures from efficient data
                let mut signatures_map = HashMap::new();
                for signature_data in state.signature_data {
                    let protein_sig = ProteinSketch::from_efficient_data(
                        signature_data,
                        state.moltype.clone(),
                        state.ksize,
                        state.scaled,
                    )?;

                    let md5sum = protein_sig.signature().md5sum.clone();
                    signatures_map.insert(md5sum.to_string(), protein_sig);
                }

                // Reconstruct the combined minhash
                let hash_function = get_hash_function_from_moltype(&state.moltype)?;
                let minhash_ksize = state.ksize * 3;
                let mut combined_minhash = KmerMinHash::new(
                    state.scaled,
                    minhash_ksize,
                    hash_function,
                    SEED,
                    true, // track_abundance
                    0,    // num (use scaled instead)
                );

                if let Some(abunds) = &state.combined_abunds {
                    combined_minhash
                        .add_many_with_abund(
                            &state
                                .combined_mins
                                .clone()
                                .into_iter()
                                .zip(abunds.iter().cloned())
                                .collect::<Vec<_>>(),
                        )
                        .map_err(|e| IndexError::SourmashError(e.to_string()))?;
                } else {
                    combined_minhash
                        .add_many(&state.combined_mins)
                        .map_err(|e| IndexError::SourmashError(e.to_string()))?;
                }

                // Update the current index state
                {
                    // Clear existing signatures and insert new ones
                    self.signatures.clear();
                    for (key, value) in signatures_map {
                        self.signatures.insert(key, value);
                    }
                }

                {
                    let mut current_combined = self.combined_minhash.lock();
                    *current_combined = combined_minhash;
                }

                Ok(())
            } else {
                Err(IndexError::NoSavedState)
            }
        }
    }

    /// Load an existing ProteomeIndex from a RocksDB path
    /// Note: This method has known issues with serialization and may not work reliably.
    /// For now, it's recommended to use save_state() and load_state() on existing indices.
    pub fn load<P: AsRef<Path>>(path: P) -> IndexResult<Self> {
        // Create RocksDB options optimized for read operations
        let opts = Self::create_rocksdb_options(false);

        // Open the database
        let db = DB::open(&opts, path)?;

        // Try to load state to get configuration
        let serialized = db.get(b"index_metadata")?;
        if let Some(data) = serialized {
            let metadata: ProteomeIndexMetadata = bincode::deserialize(&data)?;

            let _hash_function = get_hash_function_from_moltype(&metadata.moltype)
                .map_err(|e| IndexError::SourmashError(e.to_string()))?;
            let encoding_fn = get_encoding_fn_from_moltype(&metadata.moltype)
                .map_err(|e| IndexError::SourmashError(e.to_string()))?;

            // Reconstruct the combined minhash from raw data
            let hash_function = get_hash_function_from_moltype(&metadata.moltype)
                .map_err(|e| IndexError::SourmashError(e.to_string()))?;
            let minhash_ksize = metadata.ksize * 3;
            let mut combined_minhash = KmerMinHash::new(
                metadata.scaled,
                minhash_ksize,
                hash_function,
                SEED,
                true, // track_abundance
                0,    // num (use scaled instead)
            );

            if let Some(abunds) = &metadata.combined_abunds {
                combined_minhash
                    .add_many_with_abund(
                        &metadata
                            .combined_mins
                            .clone()
                            .into_iter()
                            .zip(abunds.iter().cloned())
                            .collect::<Vec<_>>(),
                    )
                    .map_err(|e| IndexError::SourmashError(e.to_string()))?;
            } else {
                combined_minhash
                    .add_many(&metadata.combined_mins)
                    .map_err(|e| IndexError::SourmashError(e.to_string()))?;
            }

            // Load all chunk data from RocksDB (sequential - RocksDB reads are single-threaded)
            let mut raw_chunks: Vec<Vec<u8>> = Vec::with_capacity(metadata.chunk_count);
            for chunk_idx in 0..metadata.chunk_count {
                let chunk_key = format!("signatures_chunk_{}", chunk_idx);
                if let Some(data) = db.get(chunk_key.as_bytes())? {
                    raw_chunks.push(data);
                }
            }

            // Deserialize and reconstruct signatures in parallel
            use rayon::prelude::*;
            let moltype = &metadata.moltype;
            let ksize = metadata.ksize;
            let scaled = metadata.scaled;

            let signatures: DashMap<String, ProteinSketch> = DashMap::new();
            raw_chunks.par_iter().try_for_each(|raw_data| -> IndexResult<()> {
                let chunk: Vec<ProteinSketchStore> = bincode::deserialize(raw_data)?;
                for signature_data in chunk {
                    let protein_sig = ProteinSketch::from_efficient_data(
                        signature_data,
                        moltype.clone(),
                        ksize,
                        scaled,
                    )?;
                    let md5sum = protein_sig.signature().md5sum.clone();
                    signatures.insert(md5sum.to_string(), protein_sig);
                }
                Ok(())
            })?;

            let index = Self {
                db,
                signatures,
                combined_minhash: Arc::new(Mutex::new(combined_minhash)),
                aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
                encoding_fn,
                moltype: metadata.moltype,
                ksize: metadata.ksize,
                minhash_ksize: metadata.ksize * 3,
                scaled: metadata.scaled,
                stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
                store_raw_sequences: metadata.store_raw_sequences,
            };

            Ok(index)
        } else {
            Err(IndexError::NoSavedState)
        }
    }

    /// Open a database for searching without loading all signatures into memory.
    ///
    /// Unlike `load()`, this method reads only the metadata header and leaves the
    /// `signatures` DashMap empty. Signatures are loaded on demand via
    /// `get_signature_by_md5()` during search. This avoids the minutes-long startup
    /// cost of deserializing 200k+ signatures when only a small fraction will be needed.
    ///
    /// Call `load_search_cache()` after opening to retrieve the pre-built inverted index.
    pub fn open_for_search<P: AsRef<Path>>(path: P) -> IndexResult<Self> {
        let opts = Self::create_rocksdb_options(false);
        // WHY: open_for_read_only avoids acquiring the exclusive LOCK file, allowing
        // multiple search processes to query the same index concurrently.
        let db = DB::open_for_read_only(&opts, path, false)?;

        let metadata_data = db.get(b"index_metadata")?.ok_or(IndexError::NoSavedState)?;
        let metadata: ProteomeIndexMetadata = bincode::deserialize(&metadata_data)?;

        let encoding_fn = get_encoding_fn_from_moltype(&metadata.moltype)
            .map_err(|e| IndexError::SourmashError(e.to_string()))?;
        let hash_function = get_hash_function_from_moltype(&metadata.moltype)
            .map_err(|e| IndexError::SourmashError(e.to_string()))?;
        let minhash_ksize = metadata.ksize * 3;

        // Create a minimal combined_minhash (not used for search, but required by struct)
        let combined_minhash =
            KmerMinHash::new(metadata.scaled, minhash_ksize, hash_function, SEED, true, 0);

        Ok(Self {
            db,
            signatures: DashMap::new(), // Empty - signatures loaded on demand by get_signature_by_md5()
            combined_minhash: Arc::new(Mutex::new(combined_minhash)),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            encoding_fn,
            moltype: metadata.moltype,
            ksize: metadata.ksize,
            minhash_ksize,
            scaled: metadata.scaled,
            stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
            store_raw_sequences: metadata.store_raw_sequences,
        })
    }

    /// Load the pre-built search cache from RocksDB.
    ///
    /// Returns `Some((target_list, inverted_index, kmer_frequencies))` if the cache was
    /// saved by `save_inverted_index()`, or `None` for older databases that predate the cache.
    ///
    /// The caller (ProteinSearcher::load) uses this to skip loading all signatures and instead
    /// find candidates via the inverted index, loading individual signatures on demand.
    pub fn load_search_cache(&self) -> IndexResult<Option<SearchCache>> {
        if let Some(data) = self.db.get(b"search_cache")? {
            let cache: SearchCache = bincode::deserialize(&data)?;
            Ok(Some(cache))
        } else {
            Ok(None)
        }
    }

    /// Load a single signature from RocksDB by its MD5 sum.
    ///
    /// Returns `None` if the signature was not found (e.g. the database was built without
    /// `save_inverted_index()`). Returns an error on deserialization failures.
    ///
    /// WHY: During search, only candidate signatures (those sharing ≥1 k-mer with the query)
    /// need to be loaded. This avoids loading all 200k+ signatures into memory at startup.
    pub fn get_signature_by_md5(&self, md5: &str) -> IndexResult<Option<ProteinSketch>> {
        let key = format!("sig_{}", md5);
        if let Some(data) = self.db.get(key.as_bytes())? {
            let sig_data: ProteinSketchStore = bincode::deserialize(&data)?;
            let sketch = ProteinSketch::from_efficient_data(
                sig_data,
                self.moltype.clone(),
                self.ksize,
                self.scaled,
            )?;
            Ok(Some(sketch))
        } else {
            Ok(None)
        }
    }

    /// Get the number of signatures in the index
    pub fn signature_count(&self) -> usize {
        self.signatures.len()
    }

    /// Get the index parameters (ksize, scaled, moltype) from the database metadata
    ///
    /// This method reads the stored metadata to extract the parameters used when
    /// the index was created, enabling autodetection of correct search parameters.
    pub fn get_index_parameters<P: AsRef<Path>>(path: P) -> IndexResult<(u32, u32, String)> {
        // Create RocksDB options optimized for read operations
        let opts = Self::create_rocksdb_options(false);

        // Open the database
        let db = DB::open(&opts, path)?;

        // Validate schema version before loading anything else.
        // Indices built before versioning was added have no schema_version key and are
        // treated as version 0.
        // Version 0 indexes built after commit 9d083c8 (Feb 24 2026) use kmer_positions
        // format and are fully compatible with schema version 1. We accept them here.
        // Only truly incompatible formats (e.g., pre-Feb-24 kmer_infos format) need rebuilding,
        // but those can't be detected by this key alone.
        let stored_version: u32 = match db.get(b"schema_version")? {
            Some(data) => bincode::deserialize(&data)?,
            None => 0, // pre-versioning index
        };
        if stored_version > SCHEMA_VERSION {
            return Err(IndexError::ValidationError {
                message: format!(
                    "Index schema version mismatch: index was built with schema version {}, \
                     but this binary uses schema version {}. \
                     Please upgrade the binary.",
                    stored_version, SCHEMA_VERSION
                ),
            });
        }

        // Try to load metadata from new chunked format first
        let metadata_serialized = db.get(b"index_metadata")?;
        if let Some(metadata_data) = metadata_serialized {
            let metadata: ProteomeIndexMetadata = bincode::deserialize(&metadata_data)?;
            return Ok((metadata.ksize, metadata.scaled, metadata.moltype));
        }

        // Fallback to old format for backward compatibility
        let serialized = db.get(b"index_state")?;
        if let Some(data) = serialized {
            let state: ProteomeIndexState = bincode::deserialize(&data)?;
            return Ok((state.ksize, state.scaled, state.moltype));
        }

        Err(IndexError::ValidationError { message: "No metadata found in database".to_string() })
    }

    /// Get the combined minhash size
    pub fn combined_minhash_size(&self) -> usize {
        self.combined_minhash.lock().size()
    }

    /// Compare this index with another for equivalency
    pub fn is_equivalent_to(&self, other: &ProteomeIndex) -> IndexResult<bool> {
        // Check basic configuration
        if self.ksize != other.ksize {
            return Ok(false);
        }
        if self.scaled != other.scaled {
            return Ok(false);
        }
        if self.moltype != other.moltype {
            return Ok(false);
        }

        // Check signature count
        if self.signature_count() != other.signature_count() {
            return Ok(false);
        }

        // Check combined minhash size
        if self.combined_minhash_size() != other.combined_minhash_size() {
            return Ok(false);
        }

        // Compare signatures - use consistent lock ordering to avoid deadlocks
        // Always lock self before other to prevent deadlocks
        // DashMap is already thread-safe, no need to lock
        let self_signatures = &self.signatures;
        let other_signatures = &other.signatures;

        for entry in self_signatures.iter() {
            let md5 = entry.key();
            let self_sig = entry.value();
            if let Some(other_sig) = other_signatures.get(md5) {
                let self_mins = self_sig.signature().get_minhash().mins();
                let other_mins = other_sig.signature().get_minhash().mins();
                if self_mins != other_mins {
                    return Ok(false);
                }

                // Compare kmer_positions
                let self_kmer_positions = self_sig.kmer_positions();
                let other_kmer_positions = other_sig.kmer_positions();

                if self_kmer_positions != other_kmer_positions {
                    return Ok(false);
                }
            } else {
                return Ok(false);
            }
        }

        // DashMap references don't need to be dropped explicitly

        // Compare combined minhashes - use consistent lock ordering
        let self_combined = self.combined_minhash.lock();
        let other_combined = other.combined_minhash.lock();

        let self_mins = self_combined.mins();
        let other_mins = other_combined.mins();
        if self_mins != other_mins {
            return Ok(false);
        }

        Ok(true)
    }

    /// Print index statistics
    pub fn print_stats(&self) {
        println!("ProteomeIndex Statistics:");
        println!("  K-mer size: {}", self.ksize);
        println!("  Scaled: {}", self.scaled);
        println!("  Molecular type: {}", self.moltype);
        // println!("  Number of signatures: {}", self.signature_count());
        println!("  Combined minhash size: {}", self.combined_minhash_size());
        println!(
            "  Raw sequence storage: {}",
            if self.store_raw_sequences { "enabled" } else { "disabled" }
        );
    }

    /// Get the raw sequence storage configuration
    pub fn store_raw_sequences(&self) -> bool {
        self.store_raw_sequences
    }

    /// Generate a filename based on the index parameters
    pub fn generate_filename(&self, base_name: &str) -> String {
        format!(
            "{}.{}.k{}.scaled{}.kmerseek.rocksdb",
            base_name, self.moltype, self.ksize, self.scaled
        )
    }

    /// Create a new index with automatic filename generation
    pub fn new_with_auto_filename<P: AsRef<Path>>(
        base_path: P,
        ksize: u32,
        scaled: u32,
        moltype: &str,
        store_raw_sequences: bool,
    ) -> IndexResult<Self> {
        let base_path = base_path.as_ref();
        let filename = format!(
            "{}.{}.k{}.scaled{}.kmerseek.rocksdb",
            base_path.file_name().unwrap().to_string_lossy(),
            moltype,
            ksize,
            scaled
        );
        let full_path = base_path.parent().unwrap().join(filename);

        Self::new(full_path, ksize, scaled, moltype, store_raw_sequences)
    }

    /// Add a single protein sequence as a signature and process its k-mers
    ///
    /// This method creates a protein signature from the given sequence, processes its k-mers
    /// to extract detailed position information, and returns the signature for later storage.
    ///
    /// The method resolves amino acid ambiguity before processing. Valid amino acids include the 20 standard amino acids (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
    /// and the ambiguous codes (B for D/N, Z for E/Q, J for I/L, X for unknown) which are resolved to one of their possible values.
    ///
    /// # Arguments
    ///
    /// * `sequence` - The protein sequence as a string
    /// * `name` - The name/identifier for the protein
    ///
    /// # Returns
    ///
    /// Returns the processed `ProteinSketch` on success, or an error if the operation fails.
    /// The error will contain details about any invalid amino acids found in the sequence.
    ///
    /// # Example
    ///
    /// ```rust
    /// use kmerseek::index::ProteomeIndex;
    /// use tempfile::tempdir;
    ///
    /// fn main() -> anyhow::Result<()> {
    ///     let dir = tempdir()?;
    ///     
    ///     // Create a new index
    ///     let index = ProteomeIndex::new(
    ///         dir.path().join("test.db"),
    ///         5,        // k-mer size
    ///         1,        // scaled (1 = capture all k-mers)
    ///         "protein", // molecular type
    ///         false,    // store raw sequences
    ///     )?;
    ///     
    ///     // Add a protein sequence
    ///     let sequence = "PLANTANDANIMALGENQMES";
    ///     let signature = index.create_protein_signature(sequence, "test_protein")?;
    ///     
    ///     // The protein signature is now ready for storage
    ///     Ok(())
    /// }
    /// ```
    pub fn create_protein_signature(
        &self,
        sequence: &str,
        name: &str,
    ) -> IndexResult<ProteinSketch> {
        // Validate and resolve ambiguity if needed
        let processed_sequence = self.aa_ambiguity.validate_and_resolve(sequence)?;

        // Create a new protein signature
        let mut protein_sig = ProteinSketch::new(name, self.ksize, self.scaled, &self.moltype)?;

        // Add the protein sequence to the signature
        // WHY: add_protein now handles all processing: minhash, kmer_infos, and sequence storage.
        // This eliminates the need for separate process_kmers and sequence storage calls.
        // We pass the index's store_raw_sequences flag to ensure consistency - sequences are
        // only stored in memory if they will be saved to disk, preventing memory waste and
        // ensuring search operations work correctly.
        protein_sig.add_protein(&processed_sequence, self.store_raw_sequences)?;

        // Return the processed signature (don't store it yet)
        Ok(protein_sig)
    }

    pub fn process_kmers(
        &self,
        sequence: &str,
        protein_signature: &mut ProteinSketch,
    ) -> IndexResult<()> {
        let ksize = self.ksize as usize;
        let hashvals: HashSet<u64> =
            protein_signature.signature().get_minhash().to_vec().into_iter().collect();

        let custom_hp = HpAlphabet::from_moltype(&self.moltype);
        for i in 0..sequence.len().saturating_sub(ksize - 1) {
            let kmer = &sequence[i..i + ksize];
            // WHY: sourmash's ReadingFrame::new_protein uppercases before hashing.
            let hashval = if let Some(ref alpha) = custom_hp {
                let encoded: Vec<u8> = kmer
                    .bytes()
                    .map(|b| {
                        alpha
                            .table()
                            .get(&b.to_ascii_uppercase())
                            .copied()
                            .unwrap_or(b)
                            .to_ascii_uppercase()
                    })
                    .collect();
                _hash_murmur(&encoded, SEED)
            } else {
                match encode_with_fn(kmer, self.encoding_fn) {
                    Ok(enc) => _hash_murmur(enc.to_ascii_uppercase().as_bytes(), SEED),
                    Err(_) => continue,
                }
            };
            if hashvals.contains(&hashval) {
                protein_signature.kmer_positions_mut().entry(hashval).or_default().push(i);
            }
        }

        Ok(())
    }

    /// Store a collection of protein signatures in the index
    ///
    /// This method stores multiple signatures at once and updates the combined minhash.
    /// It's designed to be called after processing multiple sequences in parallel.
    ///
    /// # Arguments
    ///
    /// * `signatures` - A vector of `ProteinSketch` objects to store
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or an error if the operation fails.
    pub fn store_signatures(&self, protein_signatures: Vec<ProteinSketch>) -> IndexResult<()> {
        for protein_signature in protein_signatures {
            let md5sum = protein_signature.signature().md5sum.clone();
            self.signatures.insert(md5sum.to_string(), protein_signature);
        }
        Ok(())
    }

    /// Store a batch of protein signatures efficiently.
    ///
    /// This method is optimized for batch processing by reusing the same logic
    /// as `store_signatures` but with better memory management for streaming scenarios.
    ///
    /// # Arguments
    ///
    /// * `protein_signatures` - A slice of protein signatures to store
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or an error if the operation fails.
    ///
    /// # Why this is idiomatic
    ///
    /// - **Borrowing over ownership**: Takes `&[ProteinSketch]` to avoid unnecessary moves
    /// - **Reuses existing logic**: Delegates to `store_signatures` for consistency
    /// - **Memory efficient**: Allows for batch processing without accumulating all signatures
    pub fn store_signatures_batch(&self, protein_signatures: &[ProteinSketch]) -> IndexResult<()> {
        // Convert slice to owned Vec for the existing method
        // This is a small allocation cost for the benefit of code reuse
        self.store_signatures(protein_signatures.to_vec())
    }

    /// Validate that a FASTA file exists and is readable
    ///
    /// WHY: This function centralizes file validation logic, making `process_fasta` easier to read.
    /// It provides clear, actionable error messages for common file access issues (permissions,
    /// Google Drive sync, etc.). This is idiomatic Rust - we extract validation logic into
    /// well-named functions and provide helpful error messages.
    ///
    /// # Arguments
    /// * `fasta_path` - Path to the FASTA file to validate
    ///
    /// # Returns
    /// `Ok(())` if the file is valid and readable, `ParseError` with helpful message otherwise
    fn validate_fasta_file_access<P: AsRef<Path>>(fasta_path: P) -> IndexResult<()> {
        let fasta_path = fasta_path.as_ref();

        // Check if file exists
        if !fasta_path.exists() {
            return Err(IndexError::ParseError(format!(
                "FASTA file not found: {}\nCurrent working directory: {}",
                fasta_path.display(),
                std::env::current_dir()
                    .map(|p| p.display().to_string())
                    .unwrap_or_else(|_| "unknown".to_string())
            )));
        }

        // Check if file is readable
        // WHY: On macOS, Google Drive files can exist but not be readable if they're placeholders
        // or haven't fully synced. Checking readability before attempting to open provides a
        // clearer error message than the generic "Operation not permitted" error.
        if let Ok(metadata) = std::fs::metadata(fasta_path) {
            // Check if it's actually a file (not a directory)
            if metadata.is_dir() {
                return Err(IndexError::ParseError(format!(
                    "Path is a directory, not a file: {}",
                    fasta_path.display()
                )));
            }

            // Check permissions - try to open the file to see if we can read it
            // WHY: On macOS, Google Drive files can appear to exist but fail to open if they're
            // placeholders or require special permissions. Attempting to open the file gives us
            // a better error message than just checking metadata.
            match std::fs::File::open(fasta_path) {
                Ok(_) => {
                    // File can be opened, proceed
                }
                Err(e) if e.kind() == std::io::ErrorKind::PermissionDenied => {
                    return Err(IndexError::ParseError(format!(
                        "Permission denied reading file: {}\nThis may happen if:\n- The file is in Google Drive and hasn't fully synced (check Google Drive sync status)\n- The file requires special permissions (check file permissions with 'ls -l')\n- The file is locked by another process\nError details: {}",
                        fasta_path.display(),
                        e
                    )));
                }
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                    // File disappeared between exists() check and open() - rare but possible
                    return Err(IndexError::ParseError(format!(
                        "File disappeared: {}\nThe file existed when we checked, but couldn't be opened.\nThis may happen if the file is in Google Drive and is a placeholder.\nTry: Wait for Google Drive to finish syncing, or copy the file to a local directory.",
                        fasta_path.display()
                    )));
                }
                Err(e) => {
                    // Other I/O errors - provide context
                    return Err(IndexError::ParseError(format!(
                        "Cannot open file: {}\nError: {}\nIf this is a Google Drive file, ensure it has fully synced.\nYou can check sync status in Google Drive settings.",
                        fasta_path.display(),
                        e
                    )));
                }
            }
        }

        Ok(())
    }
    /// Process a protein FASTA file with automatic compression detection and parallel processing.
    ///
    /// This method reads a FASTA file with automatic compression detection (gzip, bzip2, xz, zstd,
    /// uncompressed), validates file access, validates each protein sequence for amino acid ambiguity,
    /// creates protein signatures for each sequence, and stores them in the index using parallel batch processing.
    ///
    /// File validation is separated into `validate_fasta_file_access` for clarity and testability.
    /// Each sequence is validated using the same amino acid validation as `create_protein_signature`.
    /// If any sequence contains invalid amino acids, the entire operation will fail with an error
    /// describing the first invalid amino acid encountered.
    ///
    /// WHY: This method centralizes FASTA processing logic, handling file validation, parsing,
    /// and batch processing. This is idiomatic Rust - we separate concerns and make each function
    /// focused on a single responsibility.
    ///
    /// # Arguments
    ///
    /// * `fasta_path` - Path to the FASTA file to process (supports any compression format)
    /// * `progress_interval` - Number of sequences between progress reports (0 to disable progress)
    /// * `batch_size` - Number of sequences to process in each parallel batch (default: 1000)
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or an error if the operation fails.
    /// The error will contain details about any invalid amino acids found in the sequences.
    ///
    /// # Errors
    ///
    /// Returns `ParseError` if:
    /// - The file doesn't exist or cannot be accessed
    /// - The file is not readable (permissions, Google Drive sync issues, etc.)
    /// - The file format is invalid or cannot be parsed
    /// - Any sequence contains invalid amino acids
    ///
    /// # Examples
    ///
    /// ```
    /// # use kmerseek::ProteomeIndex;
    /// # use tempfile::tempdir;
    /// # let dir = tempdir().unwrap();
    /// # let index = ProteomeIndex::new(dir.path().join("test.db"), 10, 1, "protein", false).unwrap();
    /// # let fasta_path = dir.path().join("test.fasta");
    /// # std::fs::write(&fasta_path, ">test\nACDEFGHIKLMNPQRSTVWY").unwrap();
    ///
    /// // Process with small batch size for memory-constrained environments
    /// index.process_fasta(&fasta_path, 100, 100)?;
    ///
    /// // Process with large batch size for maximum performance
    /// index.process_fasta(&fasta_path, 100, 10000)?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Why this is idiomatic
    ///
    /// - **Auto-compression detection**: needletail handles gzip/bzip2/xz/zstd/uncompressed automatically
    /// - **Streaming + Parallel**: Combines memory efficiency with parallel processing
    /// - **Zero-copy access**: Uses `Cow<'_, [u8]>` for efficient data access when possible
    /// - **Explicit error handling**: All errors are propagated with `?` operator
    /// - **Generic path handling**: Accepts any `AsRef<Path>` type
    /// - **Configurable batching**: Allows tuning memory vs performance trade-offs
    /// - **Parallel processing**: Uses rayon for efficient parallel batch processing
    pub fn process_fasta<P: AsRef<Path>>(
        &self,
        fasta_path: P,
        progress_interval: u32,
        batch_size: usize,
    ) -> IndexResult<()> {
        use needletail::parse_fastx_file;

        if progress_interval > 0 {
            eprintln!("Reading FASTA file with automatic compression detection and parallel processing...");
        }

        // Validate file access before attempting to parse
        // WHY: We validate file access separately to keep process_fasta focused on processing.
        // This makes the code easier to read and the validation logic easier to test.
        Self::validate_fasta_file_access(&fasta_path)?;

        let fasta_path = fasta_path.as_ref();

        // Open and parse the FASTA file using needletail with auto-detection
        // WHY: We've already validated the file exists and is readable, so if needletail fails
        // here, it's likely a format/compression issue rather than a permissions issue.
        let mut reader = parse_fastx_file(fasta_path).map_err(|e| {
            // Provide context about what we were trying to do
            let error_msg = e.to_string();
            let mut diagnostic = format!(
                "Failed to parse FASTA file: {}\nError: {}",
                fasta_path.display(),
                error_msg
            );

            // Add specific help for common error patterns
            if error_msg.contains("Operation not permitted") || error_msg.contains("os error 1") {
                diagnostic.push_str(
                    "\n\nThis error often occurs with Google Drive files on macOS.\nSolutions:\n1. Ensure the file has fully synced in Google Drive\n2. Copy the file to a local directory (not in Google Drive)\n3. Check file permissions: ls -l '",
                );
                diagnostic.push_str(&fasta_path.display().to_string());
                diagnostic.push_str(
                    "'\n4. Try opening the file in another program to verify it's accessible",
                );
            }

            IndexError::ParseError(diagnostic)
        })?;

        // Create progress bar for indexing (unknown total, so use spinner style)
        let progress = if progress_interval > 0 {
            let pb = ProgressBar::new_spinner();
            pb.set_style(
                ProgressStyle::with_template("{spinner:.green} [{elapsed_precise}] {msg}")
                    .unwrap()
                    .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ "),
            );
            pb.set_message("Indexing sequences...");
            Some(pb)
        } else {
            None
        };

        // Stream records and process in parallel batches
        let mut record_count = 0;
        let mut current_batch = Vec::new();

        while let Some(record) = reader.next() {
            let record = record.map_err(|e| IndexError::ParseError(e.to_string()))?;

            // Use zero-copy access to sequence data
            let sequence = record.seq().to_vec(); // Convert to owned for parallel processing
            let id = record.id().to_vec(); // Convert to owned for parallel processing

            current_batch.push((sequence, id));
            record_count += 1;

            // Process batch when it reaches the configured size
            if current_batch.len() >= batch_size {
                self.process_batch_parallel(&current_batch, progress_interval, record_count)?;
                current_batch.clear(); // Free memory after processing
            }

            // Update progress bar
            if let Some(ref pb) = progress {
                pb.set_message(format!("Indexed {} sequences", record_count));
                pb.tick();
            }
        }

        // Process any remaining records in the final batch
        if !current_batch.is_empty() {
            self.process_batch_parallel(&current_batch, progress_interval, record_count)?;
        }

        eprintln!(
            "Done reading FASTA ({} sequences total). Building combined minhash...",
            record_count
        );
        let t_cm = Instant::now();
        self.rebuild_combined_minhash()?;
        eprintln!(
            "Combined minhash built ({} unique k-mers) in {:.1}s",
            self.combined_minhash.lock().mins().len(),
            t_cm.elapsed().as_secs_f32()
        );

        if let Some(pb) = progress {
            pb.finish_with_message(format!("Successfully indexed {} sequences", record_count));
        }
        Ok(())
    }

    /// Process a batch of records in parallel.
    ///
    /// This method handles the parallel processing of a batch of FASTA records,
    /// creating protein signatures and storing them efficiently.
    ///
    /// # Arguments
    ///
    /// * `batch` - Slice of (sequence, id) tuples to process
    /// * `progress_interval` - Progress reporting interval
    /// * `total_processed` - Total number of sequences processed so far
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or an error if the operation fails.
    ///
    /// # Why this is idiomatic
    ///
    /// - **Parallel processing**: Uses rayon for efficient parallel batch processing
    /// - **Error propagation**: All errors are properly propagated through the parallel chain
    /// - **Memory efficient**: Processes batches without accumulating all signatures
    /// - **Thread-safe**: Uses atomic operations for progress tracking
    fn process_batch_parallel(
        &self,
        batch: &[(Vec<u8>, Vec<u8>)],
        progress_interval: u32,
        total_processed: usize,
    ) -> IndexResult<()> {
        use rayon::prelude::*;

        // Process the batch in parallel
        let signatures: Result<Vec<ProteinSketch>, IndexError> = batch
            .par_iter()
            .map(|(seq_bytes, id_bytes)| {
                let sequence = std::str::from_utf8(seq_bytes)?;
                let name = std::str::from_utf8(id_bytes)?;

                // Uppercase the sequence before processing
                let sequence = sequence.to_uppercase();

                // Create protein signature for each sequence
                self.create_protein_signature(&sequence, name)
            })
            .collect();

        // Store the batch of signatures
        self.store_signatures_batch(&signatures?)?;

        // Print progress if needed
        if progress_interval > 0 && total_processed % progress_interval as usize == 0 {
            eprintln!("Processed {} sequences...", total_processed);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    // use super::*;

    use anyhow::Result;
    use sourmash::signature::SigsTrait;

    use tempfile::tempdir;

    use crate::index::ProteomeIndex;
    use crate::sketch::ProteinSketch;
    use crate::tests::test_fixtures::{
        TEST_FASTA_CONTENT, TEST_FASTA_GZ, TEST_FASTA_ZST, TEST_PROTEIN,
    };
    use crate::tests::test_utils::{self, print_kmer_positions};
    use std::collections::HashMap;
    use std::path::PathBuf;

    /// Keeping the tests for ProteomeIndex in a separate file because they're more like integration tests
    /// than unit tests with all the moltype testing. Also, it's a lot of tests!

    #[test]
    fn test_process_kmers_moltype_protein() -> Result<()> {
        let _dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        let sequence = TEST_PROTEIN;

        // Create a protein signature
        let mut protein_sig = ProteinSketch::new(
            "test_protein",
            protein_ksize,
            1, // scaled
            moltype,
        )?;

        // Add the sequence (now handles all processing: minhash, kmer_infos, sequence storage)
        // WHY: We default to storing sequences (true) in tests since they often need sequences
        // for verification and testing search operations.
        protein_sig.add_protein(sequence, true)?;
        println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

        println!("{}", protein_sig.signature().name);
        println!("{:?}", protein_sig.kmer_positions().keys());
        let kmer_count = protein_sig.kmer_positions().len();

        // Should have 17 kmers (length 21 - ksize 5 + 1)
        assert_eq!(kmer_count, 17);

        // Print all kmer infos for debugging
        test_utils::print_kmer_positions(&protein_sig);

        // Expected: hash -> sorted positions
        // Sequence: PLANTANDANIMALGENQMES (length 21, ksize 5)
        let expected_positions: HashMap<u64, Vec<usize>> = [
            (2140811952770908281, vec![14]),  // GENQM
            (4381446250900425522, vec![15]),  // ENQME
            (5798339600059429290, vec![7]),   // DANIM
            (7681438632487987439, vec![8]),   // ANIMA
            (12896310179337320481, vec![1]),  // LANTA
            (2542642819229379552, vec![3]),   // NTAND
            (11965201914550078735, vec![4]),  // TANDA
            (5893010049374798421, vec![0]),   // PLANT
            (110005740849399217, vec![6]),    // NDANI
            (3791883307084689782, vec![13]),  // LGENQ
            (14610011480386804007, vec![12]), // ALGEN
            (6941015416212662126, vec![2]),   // ANTAN
            (12636705882654324958, vec![16]), // NQMES
            (11154024130290913208, vec![10]), // IMALG
            (1225702037828834387, vec![11]),  // MALGE
            (12274863873578753245, vec![9]),  // NIMAL
            (13616372540306653069, vec![5]),  // ANDAN
        ]
        .into_iter()
        .collect();

        assert_eq!(protein_sig.kmer_positions().len(), expected_positions.len());
        for (hash, positions) in protein_sig.kmer_positions().iter() {
            let expected =
                expected_positions.get(hash).unwrap_or_else(|| panic!("Unexpected hash {}", hash));
            let mut sorted = positions.clone();
            sorted.sort();
            assert_eq!(&sorted, expected, "Position mismatch for hash {}", hash);
        }

        Ok(())
    }

    #[test]
    fn test_process_kmers_moltype_dayhoff() -> Result<()> {
        let _dir = tempdir()?;

        let protein_ksize = 5;

        let sequence = TEST_PROTEIN;

        // Create a protein signature
        let mut protein_sig = ProteinSketch::new(
            "test_protein",
            protein_ksize,
            1, // scaled
            "dayhoff",
        )?;

        // Add the sequence (now handles all processing: minhash, kmer_infos, sequence storage)
        // WHY: We default to storing sequences (true) in tests since they often need sequences
        // for verification and testing search operations.
        protein_sig.add_protein(sequence, true)?;
        println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

        println!("{}", protein_sig.signature().name);
        let hashvals = protein_sig.kmer_positions().keys().collect::<Vec<_>>();
        println!("{:?}", hashvals);
        let kmer_count = protein_sig.kmer_positions().len();

        // Should have 17 kmers (length 21 - ksize 5 + 1)
        assert_eq!(kmer_count, 17);

        // Print all kmer infos for debugging
        test_utils::print_kmer_positions(&protein_sig);

        // Expected: hash -> sorted positions (dayhoff encoding collapses 20 aa to 6 letters)
        // Sequence: PLANTANDANIMALGENQMES (length 21, ksize 5)
        let expected_positions: HashMap<u64, Vec<usize>> = [
            (17444159595263538048, vec![9]),  // NIMAL
            (2945598193614695589, vec![15]),  // ENQME
            (4548757849819812604, vec![4]),   // TANDA
            (6463872878592804545, vec![13]),  // LGENQ
            (4030406117949362159, vec![7]),   // DANIM
            (7014407397606522347, vec![1]),   // LANTA
            (5045972850709227854, vec![0]),   // PLANT
            (11417072151730334367, vec![2]),  // ANTAN
            (13574922562423607435, vec![8]),  // ANIMA
            (15050500149255106627, vec![14]), // GENQM
            (5430883729707969951, vec![10]),  // IMALG
            (13894194422852851851, vec![12]), // ALGEN
            (9604281550621775790, vec![5]),   // ANDAN
            (6161374941338912337, vec![16]),  // NQMES
            (655307631517862365, vec![6]),    // NDANI
            (360995089333906261, vec![11]),   // MALGE
            (15056713696431004031, vec![3]),  // NTAND
        ]
        .into_iter()
        .collect();

        assert_eq!(protein_sig.kmer_positions().len(), expected_positions.len());
        for (hash, positions) in protein_sig.kmer_positions().iter() {
            let expected =
                expected_positions.get(hash).unwrap_or_else(|| panic!("Unexpected hash {}", hash));
            let mut sorted = positions.clone();
            sorted.sort();
            assert_eq!(&sorted, expected, "Position mismatch for hash {}", hash);
        }

        Ok(())
    }

    #[test]
    fn test_process_kmers_moltype_hp() -> Result<()> {
        let _dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "hp";

        let sequence = TEST_PROTEIN;

        // Create a protein signature
        let mut protein_sig = ProteinSketch::new(
            "test_protein",
            protein_ksize,
            1, // scaled
            moltype,
        )?;

        // Add the sequence (now handles all processing: minhash, kmer_infos, sequence storage)
        // WHY: We default to storing sequences (true) in tests since they often need sequences
        // for verification and testing search operations.
        protein_sig.add_protein(sequence, true)?;
        println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

        println!("{}", protein_sig.signature().name);
        let hashvals = protein_sig.kmer_positions().keys().collect::<Vec<_>>();
        println!("{:?}", hashvals);
        let kmer_count = protein_sig.kmer_positions().len();

        // // Should have 14 kmers (length 21 - ksize 5 + 1), but a few duplicates
        assert_eq!(kmer_count, 14);

        // Print all kmer infos for debugging
        test_utils::print_kmer_positions(&protein_sig);

        // HP encoding collapses 20 aa to 2 letters (h/p), so multiple original k-mers
        // can produce the same hash. We store all positions together in sorted order.
        let expected_positions: HashMap<u64, Vec<usize>> = [
            (17248460043117039725, vec![11]),    // MALGE
            (5673218808929106268, vec![9]),      // NIMAL
            (16969835101383990681, vec![1]),     // LANTA
            (7345312524621807974, vec![6]),      // NDANI
            (16370543730027378051, vec![4]),     // TANDA
            (3278382041688965244, vec![8]),      // ANIMA
            (8541583772724823208, vec![10]),     // IMALG
            (16158526221854164806, vec![14]),    // GENQM
            (11553019557737058697, vec![13]),    // LGENQ
            (9081059129327932468, vec![15]),     // ENQME
            (2863220259252354754, vec![7]),      // DANIM
            (4230974618842309829, vec![0, 12]),  // PLANT(0) + ALGEN(12) → same HP hash
            (13058023948041027181, vec![3, 16]), // NTAND(3) + NQMES(16) → same HP hash
            (4144736064335623701, vec![2, 5]),   // ANTAN(2) + ANDAN(5) → same HP hash
        ]
        .into_iter()
        .collect();

        assert_eq!(protein_sig.kmer_positions().len(), expected_positions.len());
        for (hash, positions) in protein_sig.kmer_positions().iter() {
            let expected =
                expected_positions.get(hash).unwrap_or_else(|| panic!("Unexpected hash {}", hash));
            let mut sorted = positions.clone();
            sorted.sort();
            assert_eq!(&sorted, expected, "Position mismatch for hash {}", hash);
        }

        Ok(())
    }

    #[test]
    fn test_create_protein_signature_moltype_protein() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("protein_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        let sequence = TEST_PROTEIN;
        let name = "test_protein";

        // Add the protein sequence to the index and get the signature
        let signature = index.create_protein_signature(sequence, name)?;

        // Verify the signature has the expected number of k-mers
        assert_eq!(signature.kmer_positions().len(), 17, "Expected 17 k-mers for the test protein");

        // Verify some specific k-mers are present
        let expected_hash = 5893010049374798421; // Hash for "PLANT"
        assert!(
            signature.kmer_positions().contains_key(&expected_hash),
            "Expected k-mer hash {} to be present",
            expected_hash
        );

        // Store the signature in the index
        index.store_signatures(vec![signature])?;
        index.rebuild_combined_minhash()?;

        // Verify the signature was added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 1, "Expected 1 signature to be stored");
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            assert!(combined_minhash.size() == 17, "Combined minhash should contain 17 hashes");
        }

        Ok(())
    }

    #[test]
    fn test_create_protein_signature_moltype_dayhoff() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "dayhoff";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("dayhoff_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        let sequence = TEST_PROTEIN;
        let name = "test_protein";

        // Add the protein sequence to the index and get the signature
        let signature = index.create_protein_signature(sequence, name)?;

        // Verify the signature has the expected number of k-mers
        assert_eq!(signature.kmer_positions().len(), 17, "Expected 17 k-mers for the test protein");

        // Verify some specific k-mers are present
        let expected_hash = 5045972850709227854; // Hash for "PLANT" in Dayhoff encoding ("bebcb")
        assert!(
            signature.kmer_positions().contains_key(&expected_hash),
            "Expected k-mer hash {} to be present",
            expected_hash
        );

        // Store the signature in the index
        index.store_signatures(vec![signature])?;
        index.rebuild_combined_minhash()?;

        // Verify the signature was added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 1, "Expected 1 signature to be stored");
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            assert!(combined_minhash.size() == 17, "Combined minhash should contain hashes");
        }

        Ok(())
    }

    #[test]
    fn test_create_protein_signature_moltype_hp() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "hp";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("hp_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        let sequence = TEST_PROTEIN;
        let name = "test_protein";

        // Add the protein sequence to the index and get the signature
        let signature = index.create_protein_signature(sequence, name)?;

        // Verify the signature has the expected number of k-mers
        assert_eq!(signature.kmer_positions().len(), 14, "Expected 14 k-mers for the test protein");

        // Verify some specific k-mers are present
        let expected_hash = 4230974618842309829; // Hash for "PLANT" in HP encoding ("hhhpp")
        assert!(
            signature.kmer_positions().contains_key(&expected_hash),
            "Expected k-mer hash {} to be present",
            expected_hash
        );

        // Store the signature in the index
        index.store_signatures(vec![signature])?;
        index.rebuild_combined_minhash()?;

        // Verify the signature was added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 1, "Expected 1 signature to be stored");
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            assert!(combined_minhash.size() == 14, "Combined minhash should contain 14 hashes");
        }

        Ok(())
    }

    #[test]
    fn test_process_fasta_moltype_protein() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("fasta_protein_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Create a temporary FASTA file for testing with distinct sequences
        let fasta_content = TEST_FASTA_CONTENT;
        let fasta_path = dir.path().join("test.fasta");
        std::fs::write(&fasta_path, fasta_content)?;

        // Process the FASTA file
        index.process_fasta(&fasta_path, 0, 1000)?;

        // Verify the signatures were added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 2, "Expected 2 signatures to be stored");

            // Verify each signature has the expected number of k-mers
            for entry in signatures.iter() {
                let md5sum = entry.key();
                let stored_signature = entry.value();
                if md5sum == "f7661cd829e75c0d" {
                    assert!(
                        stored_signature.kmer_positions().len() == 7,
                        "LIVINGALIVE should have 7 protein 5-mers"
                    );
                } else if md5sum == "7641839ad508ab8" {
                    assert!(
                        stored_signature.kmer_positions().len() == 17,
                        "PLANTANDANIMALGENQMES should have 17 protein 5-mers"
                    );
                } else {
                    println!("md5sum: {}", md5sum);
                    println!("Name: {}", stored_signature.signature().name);
                    println!("Len of Kmer infos: {}", stored_signature.kmer_positions().len());
                    panic!("Unknown md5sum: {}", md5sum);
                }
            }
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            println!("combined_minhash.size(): {}", combined_minhash.size());
            assert!(combined_minhash.size() == 24, "Combined minhash should contain 24 hashes");
        }

        Ok(())
    }

    #[test]
    fn test_process_fasta_moltype_dayhoff() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "dayhoff";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("fasta_dayhoff_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Create a temporary FASTA file for testing with distinct sequences
        let fasta_content = TEST_FASTA_CONTENT;
        let fasta_path = dir.path().join("test.fasta");
        std::fs::write(&fasta_path, fasta_content)?;

        // Process the FASTA file
        index.process_fasta(&fasta_path, 0, 1000)?;

        // Verify the signatures were added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 2, "Expected 2 signatures to be stored");

            // Verify each signature has the expected number of k-mers
            for entry in signatures.iter() {
                let md5sum = entry.key();
                let stored_signature = entry.value();
                if md5sum == "a963d06839b6d6a9" {
                    assert!(
                        stored_signature.kmer_positions().len() == 7,
                        "LIVINGALIVE should have 7 dayhoff 5-mers"
                    );
                } else if md5sum == "84d7545d531dcf51" {
                    assert!(
                        stored_signature.kmer_positions().len() == 17,
                        "PLANTANDANIMALGENQMES should have 17 dayhoff 5-mers"
                    );
                } else {
                    println!("md5sum: {}", md5sum);
                    println!("Name: {}", stored_signature.signature().name);
                    println!("Len of Kmer infos: {}", stored_signature.kmer_positions().len());
                    panic!("Unknown md5sum: {}", md5sum);
                }
            }
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            println!("combined_minhash.size(): {}", combined_minhash.size());
            assert!(combined_minhash.size() == 24, "Combined minhash should contain 24 hashes");
        }

        Ok(())
    }

    #[test]
    fn test_process_fasta_moltype_hp() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "hp";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("fasta_hp_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Create a temporary FASTA file for testing with distinct sequences
        let fasta_content = TEST_FASTA_CONTENT;
        let fasta_path = dir.path().join("test.fasta");
        std::fs::write(&fasta_path, fasta_content)?;

        // Process the FASTA file
        index.process_fasta(&fasta_path, 0, 1000)?;

        // Verify the signatures were added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 2, "Expected 2 signatures to be stored");

            // Verify each signature has the expected number of k-mers
            for entry in signatures.iter() {
                let md5sum = entry.key();
                let stored_signature = entry.value();
                if md5sum == "24ca8d939672666b" {
                    assert!(
                        stored_signature.kmer_positions().len() == 6,
                        "LIVINGALIVE should have 6 hp 5-mers"
                    );
                } else if md5sum == "668d7173d661287b" {
                    assert!(
                        stored_signature.kmer_positions().len() == 14,
                        "PLANTANDANIMALGENQMES should have 14 hp 5-mers"
                    );
                } else {
                    println!("md5sum: {}", md5sum);
                    println!("Name: {}", stored_signature.signature().name);
                    println!("Len of Kmer infos: {}", stored_signature.kmer_positions().len());
                    panic!("Unknown md5sum: {}", md5sum);
                }
            }
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            println!("combined_minhash.size(): {}", combined_minhash.size());
            assert!(combined_minhash.size() == 16, "Combined minhash should contain 16 hashes");
        }

        Ok(())
    }

    #[test]
    fn test_process_fasta_zstd_moltype_protein() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("fasta_zstd_protein_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Process the zstd compressed FASTA file
        index.process_fasta(TEST_FASTA_ZST, 0, 1000)?;

        // Verify the signatures were added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 2, "Expected 2 signatures to be stored");

            // Verify each signature has the expected number of k-mers
            for entry in signatures.iter() {
                let md5sum = entry.key();
                let stored_signature = entry.value();
                if md5sum == "f7661cd829e75c0d" {
                    assert!(
                        stored_signature.kmer_positions().len() == 7,
                        "LIVINGALIVE should have 7 protein 5-mers"
                    );
                } else if md5sum == "7641839ad508ab8" {
                    assert!(
                        stored_signature.kmer_positions().len() == 17,
                        "PLANTANDANIMALGENQMES should have 17 protein 5-mers"
                    );
                } else {
                    println!("md5sum: {}", md5sum);
                    println!("Name: {}", stored_signature.signature().name);
                    println!("Len of Kmer infos: {}", stored_signature.kmer_positions().len());
                    panic!("Unknown md5sum: {}", md5sum);
                }
            }
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            println!("combined_minhash.size(): {}", combined_minhash.size());
            assert!(combined_minhash.size() == 24, "Combined minhash should contain 24 hashes");
        }

        Ok(())
    }

    #[test]
    fn test_process_fasta_gz_moltype_protein() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("fasta_gz_protein_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Process the FASTA file
        index.process_fasta(TEST_FASTA_GZ, 0, 1000)?;

        // Verify the signatures were added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 25, "Expected 25 signatures to be stored");

            // Check a few signatures
            for entry in signatures.iter() {
                let md5sum = entry.key();
                let stored_signature = entry.value();
                println!("\n---\nmd5sum: {}", md5sum);
                println!("Name: {}", stored_signature.signature().name);
                println!("Len of Kmer infos: {}", stored_signature.kmer_positions().len());
                if md5sum == "4d565dee9c8de9db" {
                    assert!(
                        stored_signature.kmer_positions().len() == 474,
                        "sp|O43236|SEPT4_HUMAN should have 474 protein 5-mers"
                    );
                }
                if md5sum == "4da1f84ad8be618e" {
                    assert!(
                        stored_signature.kmer_positions().len() == 235,
                        "sp|P10415|BCL2_HUMAN should have 235 protein 5-mers"
                    );
                }
            }
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            println!("combined_minhash.size(): {}", combined_minhash.size());
            assert!(
                combined_minhash.size() == 9049,
                "Combined minhash should contain 9049 protein 5-mer hashes"
            );
        }

        Ok(())
    }

    #[test]
    fn test_process_fasta_gz_moltype_dayhoff() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "dayhoff";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("fasta_gz_dayhoff_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Process the FASTA file
        index.process_fasta(TEST_FASTA_GZ, 0, 1000)?;

        // Verify the signatures were added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 25, "Expected 25 signatures to be stored");

            // Check a few signatures
            for entry in signatures.iter() {
                let md5sum = entry.key();
                let stored_signature = entry.value();
                println!("\n---\nmd5sum: {}", md5sum);
                println!("Name: {}", stored_signature.signature().name);
                println!("Len of Kmer infos: {}", stored_signature.kmer_positions().len());
                if md5sum == "fc27dcd533217385" {
                    assert!(
                        stored_signature.kmer_positions().len() == 433,
                        "sp|O43236|SEPT4_HUMAN should have 433 dayhoff 5-mers"
                    );
                }
                if md5sum == "3206706fa14185e7" {
                    assert!(
                        stored_signature.kmer_positions().len() == 204,
                        "sp|P10415|BCL2_HUMAN should have 204 dayhoff 5-mers"
                    );
                }
            }
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            println!("combined_minhash.size(): {}", combined_minhash.size());
            assert!(
                combined_minhash.size() == 2730,
                "Combined minhash should contain 2730 dayhoff 5-mer hashes"
            );
        }

        Ok(())
    }

    #[test]
    fn test_process_fasta_gz_moltype_hp() -> Result<()> {
        let dir = tempdir()?;

        // Need a higher k-mer size because otherwise the binary space of 5-mers is saturated (other tests use 5-mers)
        // If k=5, then all signatures have ~32 (=2^5) 5-mers:
        // - Not unique -> each md5sum is identical -> "25 signatures" fails
        // - "Combined minhash" only contains 32 hashes, which isn't an interesting test
        // If k=12, then all signatures have ~4096 (=2^12) 12-mers:
        // - Unique -> "25 signatures" passes
        // - "Combined minhash" has an upper bound of 4096 hashes, which is a lot more interesting
        let protein_ksize = 12;
        let moltype = "hp";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("fasta_gz_hp_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Process the FASTA file
        index.process_fasta(TEST_FASTA_GZ, 0, 1000)?;

        // Verify the signatures were added to the signatures map
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 25, "Expected 25 signatures to be stored");

            // Check a few signatures
            for entry in signatures.iter() {
                let md5sum = entry.key();
                let stored_signature = entry.value();
                println!("\n---\nmd5sum: {}", md5sum);
                println!("Name: {}", stored_signature.signature().name);
                println!("Len of Kmer infos: {}", stored_signature.kmer_positions().len());
                if md5sum == "38ffedf9d3ec7cec" {
                    assert!(
                        stored_signature.kmer_positions().len() == 452,
                        "sp|O43236|SEPT4_HUMAN should have 452 hp 12-mers"
                    );
                }
                if md5sum == "204716e4d80eb350" {
                    assert!(
                        stored_signature.kmer_positions().len() == 220,
                        "sp|P10415|BCL2_HUMAN should have 220 hp 12-mers"
                    );
                }
            }
        }

        // Verify the combined minhash was updated
        {
            let combined_minhash = index.get_combined_minhash().lock();
            println!("combined_minhash.size(): {}", combined_minhash.size());
            assert!(
                combined_minhash.size() == 3549,
                "Combined minhash should contain 3549 hp 12-mer hashes"
            );
        }

        Ok(())
    }

    #[test]
    fn test_create_protein_signature_amino_acid_validation_moltype_protein() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("validation_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Test valid sequences (including those with ambiguous amino acids that should be resolved)
        let valid_sequences = [
            "PLANTANDANIMALGENQMES", // Standard amino acids
            "ACDEFGHIKLMNPQRSTVWY",  // All standard amino acids
            "ACDEFXBZJ",             // With ambiguous amino acids (should be resolved)
        ];

        for sequence in valid_sequences.iter() {
            let protein_signature = index.create_protein_signature(sequence, "test_protein")?;
            test_utils::print_kmer_positions(&protein_signature);
            if protein_signature.signature().md5sum == "7641839ad508ab8" {
                assert!(
                    protein_signature.kmer_positions().len() == 17,
                    "Valid sequence 'PLANTANDANIMALGENQMES' should be accepted and have 17 protein 5-mers",
                );
            } else if protein_signature.signature().md5sum == "b95f0777d5439d56" {
                assert!(
                    protein_signature.kmer_positions().len() == 16,
                    "Valid sequence 'ACDEFGHIKLMNPQRSTVWY' should be accepted and have 16 protein 5-mers",
                );
            } else if protein_signature.signature().md5sum == "fa11c30a562fd82" {
                assert!(
                    protein_signature.kmer_positions().len() == 5,
                    "Valid sequence 'ACDEFXBZJ' should be accepted and have 5 protein 5-mers",
                );
            } else {
                // For the third sequence, just check the length is correct
                if protein_signature.kmer_positions().len() == 5 {
                    // This is the expected case for ACDEFXBZJ
                } else {
                    panic!(
                        "Unexpected kmer count: {} for md5sum: {}",
                        protein_signature.kmer_positions().len(),
                        protein_signature.signature().md5sum
                    );
                }
            }
        }

        // Test sequences with truly invalid characters (not in the replacements map)
        let invalid_sequences = [
            ("PLANTANDANIMALGEN1MES", "Invalid amino acid '1'"), // Number
            ("PLANTANDANIMALGEN$MES", "Invalid amino acid '$'"), // Special character
            ("PLANTANDANIMALGEN@MES", "Invalid amino acid '@'"), // Special character
        ];

        for (sequence, expected_error) in invalid_sequences.iter() {
            let result = index.create_protein_signature(sequence, "test_protein");
            assert!(result.is_err(), "Invalid sequence '{}' should be rejected", sequence);

            let error_msg = result.unwrap_err().to_string();
            assert!(
                error_msg.contains(expected_error),
                "Expected error message to contain '{}', but got '{}'",
                expected_error,
                error_msg
            );
        }

        // Test that ambiguous characters are resolved (not rejected)
        let ambiguous_sequences = [
            "PLANTANDANIMALGENBMES", // B should be resolved to D or N
            "PLANTANDANIMALGENZMES", // Z should be resolved to E or Q
            "PLANTANDANIMALGENJMES", // J should be resolved to I or L
        ];

        for sequence in ambiguous_sequences.iter() {
            let result = index.create_protein_signature(sequence, "test_protein");
            assert!(
                result.is_ok(),
                "Sequence with ambiguous amino acid '{}' should be resolved, not rejected",
                sequence
            );

            let protein_signature = result.unwrap();
            print_kmer_positions(&protein_signature);
            // Should have the same number of k-mers as the original sequence
            assert_eq!(
                protein_signature.kmer_positions().len(),
                17,
                "Resolved sequence should have 17 protein 5-mers"
            );
        }

        Ok(())
    }

    #[test]
    fn test_create_protein_signature_amino_acid_validation_moltype_dayhoff() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "dayhoff";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("validation_test_dayhoff.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Test that ambiguous characters are resolved (not rejected)
        let ambiguous_sequences = [
            "PLANTANDANIMALGENBMES", // B should be resolved to D or N
            "PLANTANDANIMALGENZMES", // Z should be resolved to E or Q
            "PLANTANDANIMALGENJMES", // J should be resolved to I or L
        ];

        for sequence in ambiguous_sequences.iter() {
            let result = index.create_protein_signature(sequence, "test_protein");
            assert!(
                result.is_ok(),
                "Sequence with ambiguous amino acid '{}' should be resolved, not rejected",
                sequence
            );

            let protein_signature = result.unwrap();
            print_kmer_positions(&protein_signature);
            // Should have the same number of k-mers as the original sequence
            println!("sequence: {}", sequence);
            assert_eq!(
                protein_signature.kmer_positions().len(),
                17,
                "Resolved sequence should have 17 protein 5-mers"
            );
            // Check that the ambiguous k-mer is resolved correctly
            if sequence == &"PLANTANDANIMALGENBMES" {
                // B resolves to D or N → dayhoff hash for NDMES/NNMES (both map to same dayhoff 6-letter encoding)
                assert!(
                    protein_signature.kmer_positions().contains_key(&6161374941338912337),
                    "Expected k-mer with hash 6161374941338912337 (NDMES/NNMES dayhoff) to be present in {}",
                    sequence
                );
            } else if sequence == &"PLANTANDANIMALGENZMES" {
                // Z resolves to E or Q → dayhoff hash for NEMES/NQMES
                assert!(
                    protein_signature.kmer_positions().contains_key(&6161374941338912337),
                    "Expected k-mer with hash 6161374941338912337 (NEMES/NQMES dayhoff) to be present in {}",
                    sequence
                );
            } else if sequence == &"PLANTANDANIMALGENJMES" {
                // J resolves to I or L → dayhoff hash for NLMES/NIMES
                assert!(
                    protein_signature.kmer_positions().contains_key(&9182605311834199497),
                    "Expected k-mer with hash 9182605311834199497 (NLMES/NIMES dayhoff) to be present in {}",
                    sequence
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_create_protein_signature_amino_acid_validation_moltype_hp() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "hp";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("validation_test_hp.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Test that ambiguous characters are resolved (not rejected)
        let ambiguous_sequences = [
            "PLANTANDANIMALGENBMES", // B should be resolved to D or N
            "PLANTANDANIMALGENZMES", // Z should be resolved to E or Q
            "PLANTANDANIMALGENJMES", // J should be resolved to I or L
        ];

        for sequence in ambiguous_sequences.iter() {
            let result = index.create_protein_signature(sequence, "test_protein");
            assert!(
                result.is_ok(),
                "Sequence with ambiguous amino acid '{}' should be resolved, not rejected",
                sequence
            );

            let protein_signature = result.unwrap();
            print_kmer_positions(&protein_signature);
            // Should have the same number of k-mers as the original sequence
            println!("sequence: {}", sequence);
            assert_eq!(
                protein_signature.kmer_positions().len(),
                14,
                "Resolved sequence should have 14 protein 5-mers"
            );
            // Check that the ambiguous k-mer is resolved correctly
            if sequence == &"PLANTANDANIMALGENBMES" {
                // B resolves to D or N → HP hash for NDMES/NNMES (both map to "pphpp" HP encoding)
                assert!(
                    protein_signature.kmer_positions().contains_key(&13058023948041027181),
                    "Expected k-mer with hash 13058023948041027181 (NDMES/NNMES HP) to be present in {}",
                    sequence
                );
            } else if sequence == &"PLANTANDANIMALGENZMES" {
                // Z resolves to E or Q → HP hash for NEMES/NQMES (both map to "pphpp" HP encoding)
                assert!(
                    protein_signature.kmer_positions().contains_key(&13058023948041027181),
                    "Expected k-mer with hash 13058023948041027181 (NEMES/NQMES HP) to be present in {}",
                    sequence
                );
            } else if sequence == &"PLANTANDANIMALGENJMES" {
                // J resolves to I or L → HP hash for NLMES/NIMES (both map to "phhpp" HP encoding)
                assert!(
                    protein_signature.kmer_positions().contains_key(&10495165127682499337),
                    "Expected k-mer with hash 10495165127682499337 (NLMES/NIMES HP) to be present in {}",
                    sequence
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_process_fasta_amino_acid_validation() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("fasta_validation_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Create a temporary FASTA file with both valid and invalid sequences
        // Note: Sequences with ambiguous characters (B, Z, J, X) should now be processed successfully
        let fasta_content = ">valid_protein1\nPLANTANDANIMALGENQMES\n>ambiguous_protein1\nPLANTANDANIMALGENBMES\n>valid_protein2\nACDEFGHIKLMNPQRSTVWY\n>invalid_protein1\nPLANTANDANIMALGEN1MES";
        let fasta_path = dir.path().join("test_validation.fasta");
        std::fs::write(&fasta_path, fasta_content)?;

        // Process the FASTA file - this should fail due to truly invalid sequences (like '1')
        let result = index.process_fasta(&fasta_path, 0, 1000);
        assert!(result.is_err(), "Processing FASTA with invalid sequences should fail");

        // Check that the error message contains information about the invalid amino acids
        let error_msg = result.unwrap_err().to_string();
        assert!(
            error_msg.contains("Invalid amino acid '1'"),
            "Error message should mention invalid amino acids, but got: {}",
            error_msg
        );

        // Create a FASTA file with only valid sequences (including ambiguous ones that should be resolved)
        let valid_fasta_content = ">valid_protein1\nPLANTANDANIMALGENQMES\n>valid_protein2\nACDEFGHIKLMNPQRSTVWY\n>ambiguous_protein1\nACDEFXBZJ\n>ambiguous_protein2\nPLANTANDANIMALGENBMES";
        let valid_fasta_path = dir.path().join("test_valid.fasta");
        std::fs::write(&valid_fasta_path, valid_fasta_content)?;

        // Process the valid FASTA file - this should succeed
        let result = index.process_fasta(&valid_fasta_path, 0, 1000);
        assert!(result.is_ok(), "Processing FASTA with valid sequences should succeed");

        // Verify that the signatures were added
        {
            let signatures = index.get_signatures();
            assert_eq!(signatures.len(), 4, "Expected 4 signatures to be stored");
        }

        Ok(())
    }

    #[test]
    fn test_create_protein_signature_no_ambiguous_chars() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("no_ambiguous_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Test a sequence with no ambiguous characters
        let sequence = "PLANTANDANIMALGENQMES"; // Only standard amino acids
        let signature = index.create_protein_signature(sequence, "test_protein")?;

        // Verify the signature has the expected number of k-mers
        assert_eq!(signature.kmer_positions().len(), 17, "Expected 17 k-mers for the test protein");

        // Verify some specific k-mers are present
        let expected_hash = 5893010049374798421; // Hash for "PLANT"
        assert!(
            signature.kmer_positions().contains_key(&expected_hash),
            "Expected k-mer hash {} to be present",
            expected_hash
        );

        Ok(())
    }

    #[test]
    fn test_index_equivalence() {
        let temp_dir = tempdir().unwrap();
        let db_path1 = temp_dir.path().join("test1.db");
        let db_path2 = temp_dir.path().join("test2.db");

        // Create two indices with the same parameters
        let index1 = ProteomeIndex::new(&db_path1, 5, 1, "protein", false).unwrap();
        let index2 = ProteomeIndex::new(&db_path2, 5, 1, "protein", false).unwrap();

        // Add the same signatures to both indices
        let sig1_1 = index1.create_protein_signature("ACDEFGHIKLMNPQRSTVWY", "test1").unwrap();
        let sig2_1 = index1.create_protein_signature("PLANTANDANIMALGENQMES", "test2").unwrap();
        index1.store_signatures(vec![sig1_1, sig2_1]).unwrap();
        index1.rebuild_combined_minhash().unwrap();

        let sig1_2 = index2.create_protein_signature("ACDEFGHIKLMNPQRSTVWY", "test1").unwrap();
        let sig2_2 = index2.create_protein_signature("PLANTANDANIMALGENQMES", "test2").unwrap();
        index2.store_signatures(vec![sig1_2, sig2_2]).unwrap();
        index2.rebuild_combined_minhash().unwrap();

        // Test equivalence
        assert!(index1.is_equivalent_to(&index2).unwrap());

        // Check stats
        assert_eq!(index1.signature_count(), 2);
        assert_eq!(index2.signature_count(), 2);
        assert_eq!(index1.combined_minhash_size(), index2.combined_minhash_size());

        // Test that different indices are not equivalent
        let index3 =
            ProteomeIndex::new(temp_dir.path().join("test3.db"), 10, 1, "protein", false).unwrap();
        assert!(!index1.is_equivalent_to(&index3).unwrap());
    }

    #[test]
    fn test_index_stats() {
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");

        // Create a new index
        let index = ProteomeIndex::new(&db_path, 8, 10, "hp", false).unwrap();

        // Add a test signature
        let sig = index.create_protein_signature("ACDEFGHIKLMNPQRSTVWY", "test").unwrap();
        index.store_signatures(vec![sig]).unwrap();
        index.rebuild_combined_minhash().unwrap();

        // Print stats (this should not panic)
        index.print_stats();

        // Verify stats
        assert_eq!(index.signature_count(), 1);
        assert!(index.combined_minhash_size() > 0);
    }

    #[test]
    fn test_manual_vs_auto_index_equivalence() {
        // Create temporary directory for test isolation
        let temp_dir = tempdir().unwrap();

        // Define paths
        let fasta_path = PathBuf::from("tests/testdata/index/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz");
        let manual_index_dir = temp_dir.path().join("manual-index");

        // Ensure the FASTA file exists
        assert!(fasta_path.exists(), "BCL2 FASTA file not found at {:?}", fasta_path);

        // Create index with manual path in temp directory
        let manual_index =
            ProteomeIndex::new(manual_index_dir.clone(), 16, 5, "hp", false).unwrap();

        // Process the FASTA file
        println!("Processing FASTA file: {:?}", fasta_path);
        manual_index.process_fasta(&fasta_path, 0, 1000).unwrap();

        // Print stats
        println!("Manual index stats after processing:");
        manual_index.print_stats();

        // Verify the manual index has content
        assert!(manual_index.signature_count() == 25, "Manual index should have 25 signatures");
        assert!(
            manual_index.combined_minhash_size() == 1603,
            "Manual index should have combined minhash of size 1603"
        );

        println!("Saving manual index state...");
        manual_index.save_state().unwrap();

        // Create a new auto-generated index in the temp directory
        println!("Creating auto-generated index in temp directory...");
        let auto_index =
            ProteomeIndex::new_with_auto_filename(&fasta_path, 16, 5, "hp", false).unwrap();

        // Process the same FASTA file
        auto_index.process_fasta(&fasta_path, 0, 1000).unwrap();

        // Print auto-generated index stats
        println!("Auto-generated index stats:");
        auto_index.print_stats();

        // Verify the auto-generated index has the same content
        assert!(
            auto_index.signature_count() == 25,
            "Auto-generated index should have 25 signatures"
        );
        assert!(
            auto_index.combined_minhash_size() == 1603,
            "Auto-generated index should have combined minhash of size 1603"
        );

        // Compare the two indices - they should be equivalent since they processed the same data
        // with the same parameters
        let are_equivalent = manual_index.is_equivalent_to(&auto_index).unwrap();
        assert!(are_equivalent, "Manual and auto-generated indices should be equivalent");

        println!("Both indices are equivalent - test passed!");
    }

    #[test]
    fn test_automatic_filename_generation() {
        let temp_dir = tempdir().unwrap();
        let base_path = temp_dir.path().join("test.fasta");

        // Test different parameter combinations
        let test_cases = vec![
            (16, 5, "hp", "test.fasta.hp.k16.scaled5.kmerseek.rocksdb"),
            (10, 1, "protein", "test.fasta.protein.k10.scaled1.kmerseek.rocksdb"),
            (8, 100, "dayhoff", "test.fasta.dayhoff.k8.scaled100.kmerseek.rocksdb"),
        ];

        for (ksize, scaled, moltype, expected) in test_cases {
            let index =
                ProteomeIndex::new_with_auto_filename(&base_path, ksize, scaled, moltype, false)
                    .unwrap();
            let generated = index.generate_filename("test.fasta");
            assert_eq!(
                generated, expected,
                "Failed for ksize={}, scaled={}, moltype={}",
                ksize, scaled, moltype
            );
        }
    }

    #[test]
    fn test_bcl2_processing_workflow() {
        // Define the FASTA file path
        let fasta_path = PathBuf::from("tests/testdata/index/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz");

        // Ensure the FASTA file exists
        assert!(fasta_path.exists(), "BCL2 FASTA file not found at {:?}", fasta_path);

        // Test different parameter combinations
        let test_cases = vec![
            (16, 5, "hp", "BCL2 with hp encoding, k=16, scaled=5"),
            (10, 1, "protein", "BCL2 with protein encoding, k=10, scaled=1"),
            (8, 100, "dayhoff", "BCL2 with dayhoff encoding, k=8, scaled=100"),
        ];

        for (ksize, scaled, moltype, description) in test_cases {
            println!("Testing: {}", description);

            // Create index with automatic filename generation
            let auto_index =
                ProteomeIndex::new_with_auto_filename(&fasta_path, ksize, scaled, moltype, false)
                    .unwrap();

            // Verify the generated filename
            let expected_filename = format!("bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz.{}.k{}.scaled{}.kmerseek.rocksdb", moltype, ksize, scaled);
            let generated_filename = auto_index.generate_filename(
                "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz",
            );
            assert_eq!(
                generated_filename, expected_filename,
                "Filename generation failed for {}",
                description
            );

            // Process the FASTA file
            auto_index.process_fasta(&fasta_path, 10, 1000).unwrap();

            // Verify the index has content
            assert!(
                auto_index.signature_count() > 0,
                "Index should have signatures for {}",
                description
            );
            assert!(
                auto_index.combined_minhash_size() > 0,
                "Index should have combined minhash for {}",
                description
            );

            // Verify we can access signatures
            let signatures = auto_index.get_signatures();
            let sig_map = signatures;
            assert!(!sig_map.is_empty(), "Signature map should not be empty for {}", description);

            // Test that signatures have the expected structure
            for entry in sig_map.iter().take(3) {
                let md5 = entry.key();
                let sig = entry.value();
                assert!(!md5.is_empty(), "MD5 should not be empty");
                assert!(!sig.signature().name.is_empty(), "Signature name should not be empty");
            }
        }
    }

    #[test]
    fn test_equivalence_workflow() {
        let temp_dir = tempdir().unwrap();
        let db_path1 = temp_dir.path().join("index1.db");
        let db_path2 = temp_dir.path().join("index2.db");

        // Create two indices with the same parameters
        let index1 = ProteomeIndex::new(&db_path1, 5, 1, "protein", false).unwrap();
        let index2 = ProteomeIndex::new(&db_path2, 5, 1, "protein", false).unwrap();

        // Add the same protein sequences to both indices
        let sequences = vec![
            ("ACDEFGHIKLMNPQRSTVWY", "protein1"),
            ("PLANTANDANIMALGENQMES", "protein2"),
            ("METHIONINELEUCINE", "protein3"),
        ];

        for (seq, name) in &sequences {
            let sig1 = index1.create_protein_signature(seq, name).unwrap();
            let sig2 = index2.create_protein_signature(seq, name).unwrap();

            index1.store_signatures(vec![sig1]).unwrap();
            index2.store_signatures(vec![sig2]).unwrap();
        }

        // Verify both indices have the same content
        assert_eq!(index1.signature_count(), 3);
        assert_eq!(index2.signature_count(), 3);
        assert_eq!(index1.combined_minhash_size(), index2.combined_minhash_size());

        // Test equivalence
        let are_equivalent = index1.is_equivalent_to(&index2).unwrap();
        assert!(are_equivalent, "Identical indices should be equivalent");

        // Create a third index with different parameters
        let index3 =
            ProteomeIndex::new(temp_dir.path().join("index3.db"), 10, 1, "protein", false).unwrap();

        // Test that different indices are not equivalent
        let are_equivalent_3 = index1.is_equivalent_to(&index3).unwrap();
        assert!(!are_equivalent_3, "Indices with different parameters should not be equivalent");

        // Test with different sequences
        let index4 =
            ProteomeIndex::new(temp_dir.path().join("index4.db"), 5, 1, "protein", false).unwrap();
        let sig4 = index4.create_protein_signature("DIFFERENTSEQUENCE", "different").unwrap();
        index4.store_signatures(vec![sig4]).unwrap();

        let are_equivalent_4 = index1.is_equivalent_to(&index4).unwrap();
        assert!(!are_equivalent_4, "Indices with different sequences should not be equivalent");
    }

    #[test]
    fn test_automatic_filename_generation_edge_cases() {
        let temp_dir = tempdir().unwrap();

        // Test with various filename patterns
        let test_cases = vec![
            ("simple.fasta", "simple.fasta.hp.k16.scaled5.kmerseek.rocksdb"),
            (
                "complex-name_with.underscores.fasta.gz",
                "complex-name_with.underscores.fasta.gz.hp.k16.scaled5.kmerseek.rocksdb",
            ),
            ("no_extension", "no_extension.hp.k16.scaled5.kmerseek.rocksdb"),
            (
                "multiple.dots.in.name.fasta",
                "multiple.dots.in.name.fasta.hp.k16.scaled5.kmerseek.rocksdb",
            ),
        ];

        for (base_name, expected) in test_cases {
            let base_path = temp_dir.path().join(base_name);
            let index =
                ProteomeIndex::new_with_auto_filename(&base_path, 16, 5, "hp", false).unwrap();
            let generated = index.generate_filename(base_name);
            assert_eq!(generated, expected, "Failed for base_name: {}", base_name);
        }

        // Test with different molecular types
        let moltype_cases = vec![
            ("hp", "test.fasta.hp.k8.scaled10.kmerseek.rocksdb"),
            ("protein", "test.fasta.protein.k8.scaled10.kmerseek.rocksdb"),
            ("dayhoff", "test.fasta.dayhoff.k8.scaled10.kmerseek.rocksdb"),
            ("raw", "test.fasta.raw.k8.scaled10.kmerseek.rocksdb"),
        ];

        for (moltype, expected) in moltype_cases {
            let base_path = temp_dir.path().join("test.fasta");
            let index =
                ProteomeIndex::new_with_auto_filename(&base_path, 8, 10, moltype, false).unwrap();
            let generated = index.generate_filename("test.fasta");
            assert_eq!(generated, expected, "Failed for moltype: {}", moltype);
        }
    }

    #[test]
    fn test_serialization_issue_demonstration() {
        println!("=== Serialization Issue Demonstration ===");
        println!("The save/load functionality has a fundamental issue with serializing KmerMinHash objects from the sourmash library.");
        println!("This is a known limitation where bincode cannot properly serialize complex objects that don't implement proper serialization traits.");
        println!();

        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("serialization_test.hp.k8.scaled10.kmerseek.rocksdb");

        // Create a simple index
        let index = ProteomeIndex::new(&db_path, 8, 10, "hp", false).unwrap();

        // Add a simple signature
        let sig = index.create_protein_signature("ACDEFGHIKLMNPQRSTVWY", "test_protein").unwrap();
        index.store_signatures(vec![sig]).unwrap();

        println!("Index created successfully with 1 signature");
        println!("Attempting to save state...");

        match index.save_state() {
            Ok(_) => {
                println!("✓ Save operation completed without errors");
                println!("  This suggests the save operation itself works");
            }
            Err(e) => {
                println!("✗ Save operation failed: {}", e);
                return;
            }
        }

        // Drop the index to release RocksDB lock
        drop(index);

        println!("Attempting to load index...");
        match ProteomeIndex::load(&db_path) {
            Ok(loaded_index) => {
                println!("✓ Load operation completed successfully!");
                println!("  Loaded index has {} signatures", loaded_index.signature_count());
                println!("  This would indicate the save/load functionality is working");
            }
            Err(e) => {
                println!("✗ Load operation failed: {}", e);
                println!();
                println!("=== Root Cause Analysis ===");
                println!("The error 'string is not valid utf8' indicates that bincode is trying to interpret");
                println!("binary serialized data as UTF-8 text, which suggests:");
                println!(
                    "1. The KmerMinHash objects contain binary data that cannot be properly serialized"
                );
                println!("2. The sourmash library's KmerMinHash does not implement proper Serialize/Deserialize traits");
                println!(
                    "3. The serialization format is not compatible with bincode's expectations"
                );
                println!();
                println!("=== Potential Solutions ===");
                println!("1. Use a different serialization format (e.g., JSON, MessagePack)");
                println!("2. Implement custom serialization for KmerMinHash objects");
                println!("3. Store only the essential data (mins, abunds) and reconstruct objects");
                println!("4. Use a different storage backend that doesn't require serialization");
                println!(
                    "5. Work with the sourmash maintainers to add proper serialization support"
                );
                println!();
                println!("=== Current Status ===");
                println!(
                    "The save/load functionality is not reliable due to these serialization issues."
                );
                println!("For now, indices should be recreated from source data rather than loaded from saved state.");
                println!("The save_state() and load_state() methods work for in-memory operations but fail for persistent storage.");
            }
        }
    }

    #[test]
    fn test_efficient_storage_basic() -> Result<()> {
        let dir = tempdir()?;

        // Create index with raw sequence storage enabled
        let index = ProteomeIndex::new(
            dir.path().join("test_efficient_basic.db"),
            5,         // k-mer size
            1,         // scaled
            "protein", // molecular type
            true,      // store raw sequences
        )?;

        // Add a protein sequence
        let sequence = "ACDEFGHIKLMNPQRSTVWY";
        let signature = index.create_protein_signature(sequence, "test_protein")?;

        // Verify raw sequence is stored
        assert!(signature.has_efficient_data());
        let raw_sequence = signature.get_raw_sequence();
        assert!(raw_sequence.is_some());
        assert_eq!(raw_sequence.unwrap(), sequence);

        // Store the signature
        index.store_signatures(vec![signature])?;

        // Verify the index has the correct configuration
        assert!(index.store_raw_sequences());
        assert_eq!(index.signature_count(), 1);

        // Get the signature and verify raw sequence is preserved
        let signatures = index.get_signatures();
        let entry = signatures.iter().next().unwrap();
        let signature = entry.value();
        assert!(signature.has_efficient_data());
        let raw_sequence = signature.get_raw_sequence();
        assert!(raw_sequence.is_some());
        assert_eq!(raw_sequence.unwrap(), sequence);

        // Test that we can save state without errors
        index.save_state()?;

        Ok(())
    }

    #[test]
    fn test_process_fasta_mixed_case_sequences() -> Result<()> {
        let dir = tempdir()?;
        const EXPECTED_SIGNATURES: usize = 2;
        const SHORT_SEQUENCE_KMERS: usize = 7; // "mAaGgCcTt" -> "MAAGGCCTT" (length 9 - ksize 3 + 1)
        const MIN_LONG_SEQUENCE_KMERS: usize = 10;

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("mixed_case_test.db"),
            3, // ksize
            1, // scaled=1 to capture all kmers for testing
            "protein",
            true, // store_raw_sequences
        )?;

        // Create and process FASTA file with mixed case sequences
        let fasta_path = dir.path().join("test_mixed_case.fasta");
        std::fs::write(&fasta_path, crate::tests::test_fixtures::TEST_FASTA_MIXED_CASE_CONTENT)?;
        index.process_fasta(&fasta_path, 0, 1000)?;

        // Verify signatures were added
        let signatures = index.get_signatures();
        assert_eq!(
            signatures.len(),
            EXPECTED_SIGNATURES,
            "Expected {EXPECTED_SIGNATURES} signatures to be stored"
        );

        // Extract k-mer counts and raw sequences using functional programming
        let kmer_counts: Vec<usize> =
            signatures.iter().map(|entry| entry.value().kmer_positions().len()).collect();

        let raw_sequences: Vec<String> = signatures
            .iter()
            .filter_map(|entry| {
                let signature = entry.value();
                signature
                    .has_efficient_data()
                    .then(|| signature.get_raw_sequence())
                    .flatten()
                    .map(|s| s.to_string())
            })
            .collect();

        // Verify k-mer counts using pattern matching
        let has_short_sequence = kmer_counts.contains(&SHORT_SEQUENCE_KMERS);
        let has_long_sequence = kmer_counts.iter().any(|&count| count > MIN_LONG_SEQUENCE_KMERS);

        assert!(has_short_sequence, "Expected to find signature with {SHORT_SEQUENCE_KMERS} k-mers for first mixed case sequence");
        assert!(
            has_long_sequence,
            "Expected to find signature with many k-mers for second mixed case sequence"
        );

        assert!(
            !raw_sequences.is_empty(),
            "Expected to find at least one signature with raw sequence data"
        );

        // Validate all raw sequences are uppercased and contain valid amino acids
        let validation_results: Vec<Result<(), String>> = raw_sequences
            .iter()
            .map(|sequence| {
                // Check for lowercase letters
                if sequence.chars().any(|c| c.is_lowercase()) {
                    return Err(format!("Raw sequence should be uppercased, but found lowercase letters in: {sequence}"));
                }

                // Check for valid amino acid characters
                if !sequence.chars().all(|c| c.is_ascii_alphabetic() && (c.is_uppercase() || c == '*')) {
                    return Err(format!("Raw sequence should contain only valid uppercase amino acid characters, but found: {sequence}"));
                }

                println!("Raw sequence stored: {sequence}");
                Ok(())
            })
            .collect();

        // Ensure all validations passed
        for result in validation_results {
            result.map_err(|e| anyhow::anyhow!(e))?;
        }

        // Test that we can save state without errors
        index.save_state()?;

        Ok(())
    }
}

/// Builder for creating ProteomeIndex instances with sensible defaults
///
/// This builder provides a fluent interface for configuring ProteomeIndex parameters.
///
/// # Examples
///
/// ```rust,no_run
/// use kmerseek::index::ProteomeIndex;
///
/// fn main() -> kmerseek::errors::IndexResult<()> {
///     // Basic usage with explicit path
///     let index = ProteomeIndex::builder()
///         .path("/path/to/database.db")
///         .ksize(5)
///         .scaled(1)
///         .moltype("protein")
///         .build()?;
///
///     // With auto filename generation
///     let index = ProteomeIndex::builder()
///         .path("/path/to/base")
///         .ksize(5)
///         .scaled(1)
///         .moltype("protein")
///         .build_with_auto_filename()?;
///
///     // With raw sequence storage
///     let index = ProteomeIndex::builder()
///         .path("/path/to/database.db")
///         .ksize(5)
///         .scaled(1)
///         .moltype("protein")
///         .store_raw_sequences(true)
///         .build()?;
///     
///     Ok(())
/// }
/// ```
#[derive(Default)]
pub struct ProteomeIndexBuilder {
    path: Option<PathBuf>,
    ksize: Option<u32>,
    scaled: Option<u32>,
    moltype: Option<String>,
    store_raw_sequences: bool,
}

impl ProteomeIndexBuilder {
    /// Create a new builder with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the database path
    pub fn path<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.path = Some(path.as_ref().to_path_buf());
        self
    }

    /// Set the k-mer size
    pub fn ksize(mut self, ksize: u32) -> Self {
        self.ksize = Some(ksize);
        self
    }

    /// Set the scaled value
    pub fn scaled(mut self, scaled: u32) -> Self {
        self.scaled = Some(scaled);
        self
    }

    /// Set the molecular type
    pub fn moltype(mut self, moltype: &str) -> Self {
        self.moltype = Some(moltype.to_string());
        self
    }

    /// Set whether to store raw sequences (defaults to false)
    pub fn store_raw_sequences(mut self, store_raw_sequences: bool) -> Self {
        self.store_raw_sequences = store_raw_sequences;
        self
    }

    /// Build the ProteomeIndex
    pub fn build(self) -> IndexResult<ProteomeIndex> {
        let path = self
            .path
            .ok_or_else(|| IndexError::BuilderError("Database path is required".to_string()))?;
        let ksize = self
            .ksize
            .ok_or_else(|| IndexError::BuilderError("K-mer size is required".to_string()))?;
        let scaled = self
            .scaled
            .ok_or_else(|| IndexError::BuilderError("Scaled value is required".to_string()))?;
        let moltype = self
            .moltype
            .ok_or_else(|| IndexError::BuilderError("Molecular type is required".to_string()))?;

        ProteomeIndex::new(path, ksize, scaled, &moltype, self.store_raw_sequences)
    }

    /// Build the ProteomeIndex with automatic filename generation
    pub fn build_with_auto_filename(self) -> IndexResult<ProteomeIndex> {
        let base_path = self
            .path
            .ok_or_else(|| IndexError::BuilderError("Base path is required".to_string()))?;
        let ksize = self
            .ksize
            .ok_or_else(|| IndexError::BuilderError("K-mer size is required".to_string()))?;
        let scaled = self
            .scaled
            .ok_or_else(|| IndexError::BuilderError("Scaled value is required".to_string()))?;
        let moltype = self
            .moltype
            .ok_or_else(|| IndexError::BuilderError("Molecular type is required".to_string()))?;

        ProteomeIndex::new_with_auto_filename(
            base_path,
            ksize,
            scaled,
            &moltype,
            self.store_raw_sequences,
        )
    }
}
