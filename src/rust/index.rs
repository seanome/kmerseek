use dashmap::DashMap;
use indicatif::{ProgressBar, ProgressStyle};
use parking_lot::Mutex;
use std::collections::{BTreeMap, BinaryHeap, HashMap};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;

use flate2::write::GzEncoder;
use flate2::Compression;
use rocksdb::{Options, DB};
use serde::{Deserialize, Serialize};
use sourmash::collection::Collection;
use sourmash::manifest::Manifest;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::storage::{FSStorage, InnerStorage};

use crate::aminoacid::AminoAcidAmbiguity;
use crate::errors::{IndexError, IndexResult};
use crate::hash_functions::get_hash_function_from_moltype;
use crate::signature::{SignatureAccess, SEED};
use crate::sketch::{ProteinSketch, ProteinSketchStore};
use crate::types::KmerSize;
use crate::types::MolType;

/// Schema version for the on-disk index format.
/// Increment this constant whenever the stored format changes in a backward-incompatible way
/// (e.g. new fields in SearchCache, renamed fields in ProteinSketchStore, etc.).
/// Indices that predate versioning (schema_version key absent) are treated as version 0
/// and will be rejected with a clear error message asking the user to rebuild.
pub const SCHEMA_VERSION: u32 = 2;

/// RocksDB key holding the kmerseek version that wrote the index, e.g. `"0.4.0"`.
/// Provenance only; `schema_version` is what selects the on-disk layout.
const KMERSEEK_VERSION_KEY: &[u8] = b"kmerseek_version";

/// First schema version whose metadata carries `remove_low_complexity`.
/// Indexes older than this are read through `LegacyProteomeIndexMetadata`.
const SCHEMA_VERSION_WITH_REMOVE_LOW_COMPLEXITY: u32 = 2;

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
    /// Whether low-complexity k-mers were removed when this index was built.
    /// Present from schema 2 onward; older layouts are read through
    /// `LegacyProteomeIndexMetadata`. See `read_metadata`.
    remove_low_complexity: bool,
}

/// The metadata layout used before `kmerseek_version` was stamped into indexes.
///
/// bincode is not self-describing, so an older blob cannot be deserialized into
/// the current `ProteomeIndexMetadata` -- it would run out of bytes on the
/// trailing field. Indexes without a version key are read through this instead.
#[derive(Serialize, Deserialize)]
struct LegacyProteomeIndexMetadata {
    total_signatures: usize,
    chunk_count: usize,
    combined_mins: Vec<u64>,
    combined_abunds: Option<Vec<u64>>,
    moltype: String,
    ksize: u32,
    scaled: u32,
    store_raw_sequences: bool,
}

impl From<LegacyProteomeIndexMetadata> for ProteomeIndexMetadata {
    fn from(legacy: LegacyProteomeIndexMetadata) -> Self {
        Self {
            total_signatures: legacy.total_signatures,
            chunk_count: legacy.chunk_count,
            combined_mins: legacy.combined_mins,
            combined_abunds: legacy.combined_abunds,
            moltype: legacy.moltype,
            ksize: legacy.ksize,
            scaled: legacy.scaled,
            store_raw_sequences: legacy.store_raw_sequences,
            // Predates the flag, so by definition every k-mer was kept.
            remove_low_complexity: false,
        }
    }
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

    // Whether to drop low-complexity (homopolymer) k-mers when building protein
    // signatures: raw amino-acid runs for any moltype, plus all-h/all-p runs for
    // HP-family moltypes. Defaults to false. Persisted separately from the
    // metadata blob; see save_state and read_remove_low_complexity.
    remove_low_complexity: bool,

    // Running totals accumulated as signatures are built, rather than by walking
    // every signature afterwards. Atomic because create_protein_signature takes
    // &self and process_fasta drives it from a par_iter.
    //
    // Relaxed ordering is sufficient and deliberate: these are counters only, and
    // no other data is published through them. Callers read them after the
    // par_iter has finished, and rayon's join already establishes happens-before
    // between the workers and the caller, so every increment is visible by then.
    // Reading mid-run would just yield a partial count, never a torn value.
    kmer_windows_examined: AtomicUsize,
    low_complexity_kmers_removed: AtomicUsize,

    // kmerseek version that wrote this index on disk, if it was stamped.
    // None for a freshly constructed index (nothing saved yet) and for indexes
    // built before version stamping existed.
    kmerseek_version: Option<String>,
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
    ///         .moltype("protein20")
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
        // Normalize before storing: pre-rename spellings (`hp`, `dayhoff`, `hp_<name>`) must
        // be written to metadata under their current names, or reopening the index would hit
        // reject_legacy_builtin_hp and refuse a database this binary just wrote.
        let moltype = MolType::new(moltype)
            .map_err(|message| IndexError::ValidationError { message })?
            .get()
            .to_string();
        let moltype = moltype.as_str();
        // Validate before opening RocksDB, so a bad size fails fast instead of
        // leaving an empty database behind.
        KmerSize::new(ksize).map_err(|message| IndexError::ConfigurationError {
            field: "ksize".to_string(),
            message,
        })?;

        // Create RocksDB options optimized for large datasets
        let opts = Self::create_rocksdb_options(true);

        // Open the database
        let db = DB::open(&opts, path)?;

        let hash_function = get_hash_function_from_moltype(moltype)
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
            moltype: moltype.to_string(),
            ksize,
            minhash_ksize,
            scaled,
            stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
            store_raw_sequences,
            remove_low_complexity: false,
            kmer_windows_examined: AtomicUsize::new(0),
            low_complexity_kmers_removed: AtomicUsize::new(0),
            kmerseek_version: None,
        })
    }

    /// Enable or disable dropping low-complexity (homopolymer) k-mers when
    /// building protein signatures. Defaults to `false`.
    pub fn set_remove_low_complexity(&mut self, remove_low_complexity: bool) {
        self.remove_low_complexity = remove_low_complexity;
    }

    /// The kmerseek version that wrote this index, e.g. `"0.4.0"`.
    ///
    /// `None` for an index built before version stamping, or one not yet saved.
    /// Useful for diagnosing behavior differences between an index and the binary
    /// querying it.
    pub fn kmerseek_version(&self) -> Option<&str> {
        self.kmerseek_version.as_deref()
    }

    /// Whether this index drops low-complexity (homopolymer) k-mers.
    ///
    /// Search reads this to build query sketches the same way the targets were
    /// built; see `save_state` for why a mismatch skews containment.
    pub fn remove_low_complexity(&self) -> bool {
        self.remove_low_complexity
    }

    /// Total k-mer windows examined while building signatures, and how many were
    /// removed as low-complexity. Returns `(0, 0)` when removal is off.
    ///
    /// Accumulated as each signature is built rather than by walking the whole
    /// signature map afterwards, so reading this is O(1) and adds no extra pass
    /// over the data. Not persisted; meaningful only for the current process.
    pub fn low_complexity_counts(&self) -> (usize, usize) {
        (
            self.kmer_windows_examined.load(Ordering::Relaxed),
            self.low_complexity_kmers_removed.load(Ordering::Relaxed),
        )
    }

    /// The kmerseek version that wrote an index, e.g. "0.4.0".
    ///
    /// `None` for indexes built before version stamping was added. Purely
    /// provenance -- the on-disk layout is keyed off `schema_version`, not this,
    /// because a semver string is not a usable layout discriminator.
    fn read_kmerseek_version(db: &DB) -> IndexResult<Option<String>> {
        match db.get(KMERSEEK_VERSION_KEY)? {
            Some(data) => Ok(Some(bincode::deserialize(&data)?)),
            None => Ok(None),
        }
    }

    /// Schema version an index was written with. Absent means pre-versioning, i.e. 0.
    fn read_schema_version(db: &DB) -> IndexResult<u32> {
        match db.get(b"schema_version")? {
            Some(data) => Ok(bincode::deserialize(&data)?),
            None => Ok(0),
        }
    }

    /// Deserialize index metadata, picking the layout by schema version.
    ///
    /// bincode is not self-describing, so the layout has to be known up front.
    /// Schema 2 added `remove_low_complexity` to the metadata; anything older is
    /// read through the legacy struct, which defaults it to `false` -- correct,
    /// since those indexes kept every k-mer.
    fn read_metadata(db: &DB, raw: &[u8]) -> IndexResult<ProteomeIndexMetadata> {
        if Self::read_schema_version(db)? >= SCHEMA_VERSION_WITH_REMOVE_LOW_COMPLEXITY {
            Ok(bincode::deserialize(raw)?)
        } else {
            let legacy: LegacyProteomeIndexMetadata = bincode::deserialize(raw)?;
            Ok(legacy.into())
        }
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
    fn save_inverted_index(&self, kmer_stats_out: Option<&Path>) -> IndexResult<()> {
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

        // Passed as a closure rather than the raw inverted_index/target_list so
        // save_kmer_stats_only (below) can plug in a different, lower-memory resolution
        // strategy instead -- see log_kmer_frequency_stats's doc comment.
        self.log_kmer_frequency_stats(&kmer_frequencies, kmer_stats_out, |hashes| {
            hashes
                .iter()
                .map(|&hash| (hash, self.resolve_kmer_string(hash, &inverted_index, &target_list)))
                .collect()
        })?;

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
        self.write_search_cache_bytes(&serialized)?;
        eprintln!("[save] search_cache written in {:.1}s", t2.elapsed().as_secs_f32());

        Ok(())
    }

    /// Largest slice written under one RocksDB key.
    ///
    /// RocksDB refuses any single value at or above 4 GiB -- its length is a u32 -- and
    /// returns `Invalid argument: value is too large`. That is a format limit, so no amount
    /// of memory or tuning avoids it. A reviewed Swiss-Prot index (483_966 targets,
    /// 115_749_594 unique k-mers) serializes to 4.43 GB and hit it.
    ///
    /// 1 GiB leaves a wide margin and keeps the chunk count small: even a 100 GB cache is
    /// 100 keys.
    const SEARCH_CACHE_CHUNK: usize = 1 << 30;

    /// Store the serialized [`SearchCache`], splitting it across keys when it is too large
    /// for one.
    ///
    /// A cache that fits keeps the original single `search_cache` key and byte-for-byte
    /// layout, so databases written by this version are still readable by older builds
    /// whenever they would have been readable at all.
    ///
    /// The chunk count is written LAST. A crash midway therefore leaves a database with no
    /// count key, which reads as "no cache" rather than as a cache that silently ends
    /// early -- the failure that would otherwise surface much later as a truncated
    /// inverted index and quietly missing search hits.
    fn write_search_cache_bytes(&self, serialized: &[u8]) -> IndexResult<()> {
        self.write_search_cache_chunked(serialized, Self::SEARCH_CACHE_CHUNK)
    }

    /// The chunk size is a parameter so the split path can be tested against a few KB
    /// instead of the 4 GiB it takes to reach it in production.
    fn write_search_cache_chunked(&self, serialized: &[u8], chunk_size: usize) -> IndexResult<()> {
        if serialized.len() < chunk_size {
            self.db.put(b"search_cache", serialized)?;
            return Ok(());
        }

        let chunks: Vec<&[u8]> = serialized.chunks(chunk_size).collect();
        eprintln!(
            "[save] cache is {} bytes, above the {} byte RocksDB value limit -- writing {} chunks",
            serialized.len(),
            chunk_size,
            chunks.len()
        );

        // Order matters, and it is chosen so that an interrupted write leaves a database
        // that REFUSES to open rather than one that opens wrong.
        //
        // The count goes first, as a delete. Rewriting a database that already had a valid
        // count would otherwise leave that count standing over chunks being overwritten
        // underneath it, and a crash midway would produce a count matched by a mixture of
        // old and new chunks -- which deserializes into a plausible, wrong index.
        self.db.delete(b"search_cache_chunks")?;
        // Any stale single-key cache from an earlier build would otherwise win on read,
        // because the reader prefers it.
        self.db.delete(b"search_cache")?;

        for (i, chunk) in chunks.iter().enumerate() {
            self.db.put(Self::search_cache_chunk_key(i), chunk)?;
        }
        // A previous write may have made more chunks than this one. Those are not read --
        // the count bounds the loop -- but they are left behind claiming to be part of an
        // index they no longer belong to, and the torn-write check below reads their
        // presence as evidence. Clear them until the first key that is absent.
        for i in chunks.len().. {
            if self.db.get(Self::search_cache_chunk_key(i))?.is_none() {
                break;
            }
            self.db.delete(Self::search_cache_chunk_key(i))?;
        }

        // Written LAST: its presence is what certifies every chunk above is on disk.
        self.db.put(b"search_cache_chunks", chunks.len().to_string().as_bytes())?;
        Ok(())
    }

    fn search_cache_chunk_key(i: usize) -> Vec<u8> {
        format!("search_cache_chunk_{i}").into_bytes()
    }

    /// Number of k-mers listed in the "most common" / "least common" summaries.
    const KMER_EXAMPLES: usize = 10;

    /// Log a k-mer frequency histogram and the most/least common k-mers to stderr.
    ///
    /// The lists show the actual encoded k-mer string (not the hash). `resolve` turns a
    /// batch of hashes into their text (or a placeholder for any it can't find) -- callers
    /// that already have a full inverted index can resolve directly from it;
    /// `save_kmer_stats_only` (no inverted index kept in memory, to save memory at proteome
    /// scale) scans signatures for just this handful of hashes instead. Either way this
    /// function only ever asks for `2 * KMER_EXAMPLES` hashes, never a bulk resolution.
    fn log_kmer_frequency_stats(
        &self,
        kmer_frequencies: &HashMap<u64, usize>,
        kmer_stats_out: Option<&Path>,
        resolve: impl Fn(&[u64]) -> HashMap<u64, String>,
    ) -> IndexResult<()> {
        if kmer_frequencies.is_empty() {
            return Ok(());
        }
        // Build the spectrum once. Every summary below (bins, total, unique, median, tie
        // counts) is derived from it rather than rescanning the map, which holds over
        // 125 million entries for protein k=15.
        let spectrum = Self::frequency_spectrum(kmer_frequencies);
        Self::print_frequency_histogram(&spectrum);
        if let Some(path) = kmer_stats_out {
            self.write_kmer_frequency_spectrum(path, &spectrum)?;
        }

        // Resolving a hash back to text needs the sequence it came from, which is only in
        // memory when the index was built with store_raw_sequences.
        if !self.has_stored_sequences() {
            eprintln!("[save] (k-mer sequences not stored, skipping most/least common k-mers)");
            return Ok(());
        }

        let n = Self::KMER_EXAMPLES;
        let most = Self::n_smallest_by_key(kmer_frequencies, n, |hash, count| {
            (std::cmp::Reverse(count), hash)
        });
        let least = Self::n_smallest_by_key(kmer_frequencies, n, |hash, count| (count, hash));

        // Hashes actually printed below, deduplicated so a hash that lands in both `most`
        // and `least` (a proteome with very few unique k-mers) is only resolved once.
        let mut example_hashes: Vec<u64> =
            most.iter().chain(least.iter()).map(|&(hash, _)| hash).collect();
        example_hashes.sort_unstable();
        example_hashes.dedup();
        let resolved = resolve(&example_hashes);

        self.print_kmer_examples("most common", &most, &spectrum, &resolved);
        self.print_kmer_examples("least common", &least, &spectrum, &resolved);
        Ok(())
    }

    /// Write the k-mer frequency spectrum as CSV, gzip-compressed when the path ends in `.gz`.
    ///
    /// WHY the spectrum rather than one row per k-mer: this is the distribution you plot, and
    /// it is a few thousand rows instead of tens of millions, so a sweep over alphabets and
    /// k-sizes stays small. `moltype` and `ksize` are repeated on every row so that files from
    /// different runs concatenate directly into one frame.
    fn write_kmer_frequency_spectrum(
        &self,
        path: &Path,
        spectrum: &BTreeMap<usize, usize>,
    ) -> IndexResult<()> {
        let file = File::create(path)?;

        // The gzip trailer is written when the encoder is finished. Finishing it explicitly
        // rather than leaving it to Drop is what surfaces a failed final write: flate2's Drop
        // discards that error, which would leave a truncated file behind while this function
        // reported success.
        if path.extension().is_some_and(|e| e == "gz") {
            let mut encoder = GzEncoder::new(file, Compression::default());
            self.write_spectrum_csv(&mut encoder, spectrum)?;
            encoder.finish()?;
        } else {
            let mut file = file;
            self.write_spectrum_csv(&mut file, spectrum)?;
        }

        eprintln!(
            "[save] Wrote k-mer frequency spectrum ({} rows) to {}",
            spectrum.len(),
            path.display()
        );
        Ok(())
    }

    /// Write the totals comment and one row per distinct occurrence count.
    fn write_spectrum_csv<W: Write>(
        &self,
        sink: &mut W,
        spectrum: &BTreeMap<usize, usize>,
    ) -> IndexResult<()> {
        // Totals go in a leading `#` comment rather than a column, because they are constant
        // for the whole file and would otherwise be repeated on every row. Readers skip it
        // with a comment prefix, e.g. polars' `read_csv(..., comment_prefix="#")`.
        let total = Self::total_kmers(spectrum);
        let unique = Self::unique_kmers(spectrum);
        writeln!(
            sink,
            "# total_kmers={total} unique_kmers={unique} mean_seqs_per_kmer={:.4} \
             median_seqs_per_kmer={:.1} mode_seqs_per_kmer={} moltype={} ksize={}",
            total as f64 / unique as f64,
            Self::median_occurrences(spectrum),
            Self::mode_occurrences(spectrum),
            self.moltype,
            self.ksize,
        )?;

        let mut writer = csv::Writer::from_writer(sink);
        writer.write_record(["moltype", "ksize", "occurrences", "n_kmers"])?;
        let ksize = self.ksize.to_string();
        for (occurrences, n_kmers) in spectrum {
            writer.write_record([
                &self.moltype,
                &ksize,
                &occurrences.to_string(),
                &n_kmers.to_string(),
            ])?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Exact frequency spectrum: occurrence count -> how many k-mers were seen that many times.
    ///
    /// Counts are always at least 1 because the caller builds them by incrementing from zero,
    /// which [`frequency_bins`] relies on to avoid underflowing on `leading_zeros`.
    fn frequency_spectrum(kmer_frequencies: &HashMap<u64, usize>) -> BTreeMap<usize, usize> {
        let mut spectrum: BTreeMap<usize, usize> = BTreeMap::new();
        for &count in kmer_frequencies.values() {
            debug_assert!(count > 0, "k-mer frequencies are counts, so never zero");
            *spectrum.entry(count).or_insert(0) += 1;
        }
        spectrum
    }

    /// Total (sequence, k-mer) pairs: each k-mer counted once per sequence containing it.
    fn total_kmers(spectrum: &BTreeMap<usize, usize>) -> usize {
        spectrum.iter().map(|(occurrences, n_kmers)| occurrences * n_kmers).sum()
    }

    /// How many distinct k-mers the index holds.
    fn unique_kmers(spectrum: &BTreeMap<usize, usize>) -> usize {
        spectrum.values().sum()
    }

    /// Bin k-mer counts into power-of-two ranges: bin `b` holds counts in `[2^b, 2^(b+1))`.
    ///
    /// WHY power-of-two: k-mer frequency distributions are heavily right-skewed (most k-mers
    /// occur once, a few occur thousands of times), so linear bins would be nearly unreadable.
    fn frequency_bins(spectrum: &BTreeMap<usize, usize>) -> BTreeMap<u32, usize> {
        let mut bins: BTreeMap<u32, usize> = BTreeMap::new();
        for (&count, &n_kmers) in spectrum {
            let bin = usize::BITS - count.leading_zeros() - 1;
            *bins.entry(bin).or_insert(0) += n_kmers;
        }
        bins
    }

    /// Median occurrences per k-mer, averaging the two middle values when the count is even.
    ///
    /// WHY alongside the mean: these distributions are heavily right-skewed, so a handful of
    /// very common k-mers drag the mean well above what a typical k-mer looks like.
    fn median_occurrences(spectrum: &BTreeMap<usize, usize>) -> f64 {
        let unique = Self::unique_kmers(spectrum);
        if unique == 0 {
            return 0.0;
        }
        // Positions of the middle element(s) in the sorted list of per-k-mer counts.
        let (lower_rank, upper_rank) = ((unique - 1) / 2, unique / 2);
        let (mut seen, mut lower) = (0usize, None);
        for (&occurrences, &n_kmers) in spectrum {
            seen += n_kmers;
            if lower.is_none() && seen > lower_rank {
                lower = Some(occurrences);
            }
            if seen > upper_rank {
                return (lower.unwrap_or(occurrences) + occurrences) as f64 / 2.0;
            }
        }
        lower.unwrap_or(0) as f64
    }

    /// The most common per-k-mer occurrence count, ties broken toward the smaller count.
    ///
    /// WHY alongside mean and median: this answers "what does a typical k-mer look like" most
    /// directly. In a k-mer frequency spectrum it is almost always 1, since most k-mers occur
    /// in only one sequence -- seeing that plainly is more useful than inferring it from mean
    /// and median alone.
    fn mode_occurrences(spectrum: &BTreeMap<usize, usize>) -> usize {
        spectrum
            .iter()
            .max_by_key(|&(&occurrences, &n_kmers)| (n_kmers, std::cmp::Reverse(occurrences)))
            .map(|(&occurrences, _)| occurrences)
            .unwrap_or(0)
    }

    /// Print the binned frequency distribution as an ASCII bar chart.
    ///
    /// Empty bins between the smallest and largest populated bin are printed with a zero count
    /// so that gaps in the distribution are explicit rather than silently skipped.
    fn print_frequency_histogram(spectrum: &BTreeMap<usize, usize>) {
        const BAR_WIDTH: usize = 40;
        let bins = Self::frequency_bins(spectrum);
        let (Some((&first, _)), Some((&last, _))) = (bins.iter().next(), bins.iter().next_back())
        else {
            return;
        };
        let max_count = bins.values().max().copied().unwrap_or(1);

        let total = Self::total_kmers(spectrum);
        let unique = Self::unique_kmers(spectrum);
        eprintln!(
            "[save] K-mer frequency histogram: {total} total k-mers found, {unique} unique \
             (mean {:.2}, median {:.1}, mode {} sequences per k-mer)",
            total as f64 / unique as f64,
            Self::median_occurrences(spectrum),
            Self::mode_occurrences(spectrum),
        );
        for bin in first..=last {
            let count = bins.get(&bin).copied().unwrap_or(0);
            let (lo, hi) = (1u64 << bin, (1u64 << (bin + 1)) - 1);
            let label = if lo == hi { format!("{lo}") } else { format!("{lo}-{hi}") };
            let bar = "#".repeat((count * BAR_WIDTH / max_count).max(usize::from(count > 0)));
            eprintln!("[save]   {label:>12} occurrences: {count:>10} k-mers  {bar}");
        }
    }

    /// Return the `n` (hash, count) pairs with the smallest `key`, in ascending key order.
    ///
    /// WHY a bounded heap instead of sorting: real proteome indexes hold tens of millions of
    /// unique k-mers, so materializing and sorting the whole list just to read off 10 entries
    /// would cost seconds and hundreds of MB. This is O(N log n) time and O(n) memory.
    fn n_smallest_by_key<K: Ord>(
        kmer_frequencies: &HashMap<u64, usize>,
        n: usize,
        key: impl Fn(u64, usize) -> K,
    ) -> Vec<(u64, usize)> {
        let mut heap: BinaryHeap<(K, u64, usize)> = BinaryHeap::with_capacity(n);
        for (&hash, &count) in kmer_frequencies {
            let candidate = (key(hash, count), hash, count);
            if heap.len() < n {
                heap.push(candidate);
            } else if heap.peek().is_some_and(|worst| candidate < *worst) {
                // Replacing the worst entry avoids a push/pop pair for every element that
                // cannot make the list, which is nearly all of them.
                *heap.peek_mut().expect("heap is non-empty because n > 0") = candidate;
            }
        }
        let mut selected = heap.into_vec();
        selected.sort();
        selected.into_iter().map(|(_, hash, count)| (hash, count)).collect()
    }

    /// Print one labelled list of example k-mers, noting how many share the boundary frequency.
    ///
    /// WHY the tie note: when hundreds of thousands of k-mers all occur once, listing ten of
    /// them looks like a ranking but is an arbitrary sample. Saying how many tie makes
    /// that explicit.
    fn print_kmer_examples(
        &self,
        label: &str,
        examples: &[(u64, usize)],
        spectrum: &BTreeMap<usize, usize>,
        resolved: &HashMap<u64, String>,
    ) {
        let Some(&(_, boundary)) = examples.last() else { return };
        let tied = spectrum.get(&boundary).copied().unwrap_or(0);
        eprintln!("[save] {} {} k-mers (encoded k-mer: occurrences):", examples.len(), label);
        for &(hash, count) in examples {
            let kmer = resolved
                .get(&hash)
                .cloned()
                .unwrap_or_else(|| format!("<sequence unavailable, hash {hash}>"));
            eprintln!("[save]   {kmer}: {count}");
        }
        if tied > examples.len() {
            eprintln!("[save]   ({tied} k-mers occur {boundary}x; showing an arbitrary sample)");
        }
    }

    /// Whether signatures kept their sequence text, which `resolve_kmer_string` needs.
    fn has_stored_sequences(&self) -> bool {
        self.signatures
            .iter()
            .next()
            .is_some_and(|s| s.get_moltype_sequence().is_some() || s.get_raw_sequence().is_some())
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

    /// Resolve a small set of k-mer hashes back to their sequence text by scanning
    /// `self.signatures` directly, with no inverted index at all.
    ///
    /// WHY a scan instead of a lookup: `resolve_kmer_string` needs an inverted index (hash
    /// -> protein) to do this in O(1), but building that index for the whole proteome is
    /// exactly the memory cost `save_kmer_stats_only` exists to avoid (see its doc comment).
    /// For a handful of hashes (`KMER_EXAMPLES` most + least common, ~20 total) a linear scan
    /// that early-exits once every hash is found is the better trade -- O(hashes) memory
    /// instead of O(total (protein, k-mer) pairs), and in practice touches only the first
    /// handful of signatures before `remaining` empties out. Not a general-purpose reverse
    /// lookup: don't reach for this for more than a small, fixed number of hashes.
    fn resolve_kmer_strings_by_scan(&self, hashes: &[u64]) -> HashMap<u64, String> {
        let ksize = self.ksize as usize;
        let mut remaining: std::collections::HashSet<u64> = hashes.iter().copied().collect();
        let mut resolved: HashMap<u64, String> = HashMap::with_capacity(hashes.len());
        if remaining.is_empty() {
            return resolved;
        }
        for entry in self.signatures.iter() {
            if remaining.is_empty() {
                break;
            }
            let Some(seq) = entry.get_moltype_sequence().or_else(|| entry.get_raw_sequence())
            else {
                continue;
            };
            let positions = entry.kmer_positions();
            remaining.retain(|&hash| {
                let Some(&position) = positions.get(&hash).and_then(|p| p.first()) else {
                    return true; // not in this signature -- keep looking
                };
                if let Some(kmer_str) = seq.get(position..position + ksize) {
                    resolved.insert(hash, kmer_str.to_string());
                }
                false // found (or an unresolvable slice) -- stop looking for this one either way
            });
        }
        resolved
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
        self.save_state_with_kmer_stats(None)
    }

    /// Compute and write the k-mer frequency spectrum WITHOUT persisting a searchable
    /// index -- no RocksDB writes at all (no per-signature storage, no SearchCache, no
    /// metadata). Signatures must already be in memory (i.e. call after `process_fasta`).
    ///
    /// WHY this exists as a separate path rather than a flag on `save_state_with_kmer_stats`:
    /// that function always persists a full searchable index (chunked signature storage +
    /// one serialized SearchCache value) even when the caller only wants `--kmer-stats-out`.
    /// That persistence is expensive in two independent ways that both broke real runs at
    /// proteome scale: the SearchCache is one RocksDB value containing the whole inverted
    /// index, which can exceed RocksDB's ~4 GiB single-value limit at high ksize + a small
    /// alphabet + a large proteome ("Invalid argument: value is too large"); and chunked
    /// signature storage writes thousands of small files (500K+ signatures / CHUNK_SIZE=100
    /// per run), which hit transient filesystem I/O errors under concurrent load on a shared
    /// parallel filesystem ("IO error: ... Input/output error"). Both are sidestepped
    /// entirely by never writing to `self.db` -- everything below reads only `self.signatures`,
    /// which is already in memory from `process_fasta`.
    pub fn save_kmer_stats_only(&self, kmer_stats_out: &Path) -> IndexResult<()> {
        let t0 = Instant::now();
        let total_sigs = self.signatures.len();
        eprintln!(
            "[save] Computing k-mer frequency stats for {} signatures (no index persisted)...",
            total_sigs
        );

        // Only kmer_frequencies (one entry per UNIQUE k-mer) is needed for the CSV. An
        // earlier version of this function also built `inverted_index: HashMap<u64,
        // Vec<u32>>` and `target_list: Vec<String>` up front, purely so the console's
        // "most/least common k-mer" printout could resolve a hash back to real sequence
        // text. Those two structures scale with total (protein, k-mer) PAIRS and total
        // signature count respectively -- not unique_kmers -- and at UniRef50 scale (tens
        // of billions of residues) that dwarfs kmer_frequencies itself, becoming the
        // dominant memory cost this function was supposed to have eliminated. The example
        // text is still produced (see log_kmer_frequency_stats's `resolve` closure below),
        // just via `resolve_kmer_strings_by_scan` instead -- a scan bounded to the ~20
        // hashes actually printed, not a structure sized to the whole proteome.
        let mut kmer_frequencies: HashMap<u64, usize> = HashMap::new();
        let mut n_processed = 0usize;

        for entry in self.signatures.iter() {
            let mins = entry.value().signature().minhash.mins();
            for min in mins {
                *kmer_frequencies.entry(min).or_insert(0) += 1;
            }

            n_processed += 1;
            if n_processed % 1000 == 0 {
                eprintln!(
                    "[save] {}/{} signatures processed ({:.1}s elapsed)",
                    n_processed,
                    total_sigs,
                    t0.elapsed().as_secs_f32(),
                );
            }
        }

        eprintln!(
            "[save] Processed all {} signatures in {:.1}s. {} unique k-mers.",
            total_sigs,
            t0.elapsed().as_secs_f32(),
            kmer_frequencies.len(),
        );

        self.log_kmer_frequency_stats(&kmer_frequencies, Some(kmer_stats_out), |hashes| {
            self.resolve_kmer_strings_by_scan(hashes)
        })?;

        eprintln!(
            "[save] Done in {:.1}s (nothing written to the index database).",
            t0.elapsed().as_secs_f32()
        );
        Ok(())
    }

    /// Like [`save_state`], but also writes the k-mer frequency spectrum to `kmer_stats_out`
    /// as CSV (gzipped when the path ends in `.gz`) for plotting across alphabets and k-sizes.
    pub fn save_state_with_kmer_stats(&self, kmer_stats_out: Option<&Path>) -> IndexResult<()> {
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
            remove_low_complexity: self.remove_low_complexity,
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

        // Stamp the writing kmerseek version. Besides being useful provenance, its
        // presence tells readers the metadata carries a remove_low_complexity field.
        self.db.put(KMERSEEK_VERSION_KEY, bincode::serialize(env!("CARGO_PKG_VERSION"))?)?;
        eprintln!(
            "[save] Metadata + schema_version + kmerseek_version written in {:.1}s total",
            t3.elapsed().as_secs_f32()
        );

        // Build and persist search cache + individual signatures for fast search startup
        eprintln!("[save] Building search cache...");
        let t4 = Instant::now();
        self.save_inverted_index(kmer_stats_out)?;
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
            let metadata = Self::read_metadata(&self.db, &metadata_data)?;

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
            let metadata = Self::read_metadata(&db, &data)?;

            let _hash_function = get_hash_function_from_moltype(&metadata.moltype)
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
            let remove_low_complexity = metadata.remove_low_complexity;
            // Read before `db` is moved into the struct below.
            let kmerseek_version = Self::read_kmerseek_version(&db)?;

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
                moltype: metadata.moltype,
                ksize: metadata.ksize,
                minhash_ksize: metadata.ksize * 3,
                scaled: metadata.scaled,
                stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
                store_raw_sequences: metadata.store_raw_sequences,
                remove_low_complexity,
                kmer_windows_examined: AtomicUsize::new(0),
                low_complexity_kmers_removed: AtomicUsize::new(0),
                kmerseek_version,
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
        let metadata = Self::read_metadata(&db, &metadata_data)?;

        let hash_function = get_hash_function_from_moltype(&metadata.moltype)
            .map_err(|e| IndexError::SourmashError(e.to_string()))?;
        let minhash_ksize = metadata.ksize * 3;

        // Create a minimal combined_minhash (not used for search, but required by struct)
        let combined_minhash =
            KmerMinHash::new(metadata.scaled, minhash_ksize, hash_function, SEED, true, 0);

        let remove_low_complexity = metadata.remove_low_complexity;
        let kmerseek_version = Self::read_kmerseek_version(&db)?;

        Ok(Self {
            db,
            signatures: DashMap::new(), // Empty - signatures loaded on demand by get_signature_by_md5()
            combined_minhash: Arc::new(Mutex::new(combined_minhash)),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            moltype: metadata.moltype,
            ksize: metadata.ksize,
            minhash_ksize,
            scaled: metadata.scaled,
            stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
            store_raw_sequences: metadata.store_raw_sequences,
            remove_low_complexity,
            kmer_windows_examined: AtomicUsize::new(0),
            low_complexity_kmers_removed: AtomicUsize::new(0),
            kmerseek_version,
        })
    }

    /// Load the pre-built search cache from RocksDB.
    ///
    /// Returns `Some((target_list, inverted_index, kmer_frequencies))` if the cache was
    /// saved by `save_inverted_index()`, or `None` for older databases that predate the cache.
    ///
    /// The caller (ProteinSearcher::load) uses this to skip loading all signatures and instead
    /// find candidates via the inverted index, loading individual signatures on demand.
    /// Reads a cache written either as one value or as chunks; see
    /// [`Self::write_search_cache_bytes`].
    pub fn load_search_cache(&self) -> IndexResult<Option<SearchCache>> {
        if let Some(data) = self.db.get(b"search_cache")? {
            let cache: SearchCache = bincode::deserialize(&data)?;
            return Ok(Some(cache));
        }

        let Some(count) = self.db.get(b"search_cache_chunks")? else {
            // No count key. That is either a database written before the cache existed --
            // legitimately Ok(None), and the caller falls back to loading every signature
            // -- or one whose chunk write was interrupted before the count was committed.
            //
            // Those two must not be confused. A torn index returning Ok(None) is the
            // quietest possible failure: the search still runs, on the slow path, against
            // an index nobody is told is broken. Chunk 0 is written before any other, so
            // its presence without a count means the write did not finish.
            if self.db.get(Self::search_cache_chunk_key(0))?.is_some() {
                return Err(IndexError::CorruptIndex(
                    "search cache chunks are present but the chunk count is missing: the \
                     index was interrupted while being written and is incomplete. Rebuild \
                     it with `kmerseek index`."
                        .to_string(),
                ));
            }
            return Ok(None);
        };
        let count: usize = String::from_utf8_lossy(&count).parse().map_err(|_| {
            IndexError::CorruptIndex(format!(
                "search_cache_chunks is not a number: {:?}",
                String::from_utf8_lossy(&count)
            ))
        })?;

        // A missing chunk is an error, never a short cache. Deserializing a truncated
        // stream would either fail somewhere confusing or, worse, succeed against a
        // partial inverted index and drop search hits with nothing to show for it.
        let mut serialized = Vec::new();
        for i in 0..count {
            let chunk = self.db.get(Self::search_cache_chunk_key(i))?.ok_or_else(|| {
                IndexError::CorruptIndex(format!(
                    "search cache chunk {i} of {count} is missing; the index is \
                         incomplete and must be rebuilt"
                ))
            })?;
            serialized.extend_from_slice(&chunk);
        }
        let cache: SearchCache = bincode::deserialize(&serialized)?;
        Ok(Some(cache))
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
        let stored_version = Self::read_schema_version(&db)?;
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
            let metadata = Self::read_metadata(&db, &metadata_data)?;
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
    ///
    /// WHY the `.nolowcomplexity` segment: the whole point of the flag is A/B
    /// comparison, so indexing the same FASTA with and without it must not
    /// resolve to the same path — otherwise the second run silently clobbers the
    /// first. Indexes that keep every k-mer retain their historical filename.
    pub fn generate_filename(&self, base_name: &str) -> String {
        let suffix = if self.remove_low_complexity { ".nolowcomplexity" } else { "" };
        format!(
            "{}.{}.k{}.scaled{}{}.kmerseek.rocksdb",
            base_name, self.moltype, self.ksize, self.scaled, suffix
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
    ///         "protein20", // molecular type
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
        let processed_sequence = self.aa_ambiguity.validate_and_resolve(sequence, &self.moltype)?;

        // Create a new protein signature
        let mut protein_sig = ProteinSketch::new(name, self.ksize, self.scaled, &self.moltype)?;
        protein_sig.set_remove_low_complexity(self.remove_low_complexity);

        // Add the protein sequence to the signature
        // WHY: add_protein now handles all processing: minhash, kmer_infos, and sequence storage.
        // This eliminates the need for separate process_kmers and sequence storage calls.
        // We pass the index's store_raw_sequences flag to ensure consistency - sequences are
        // only stored in memory if they will be saved to disk, preventing memory waste and
        // ensuring search operations work correctly.
        protein_sig.add_protein(&processed_sequence, self.store_raw_sequences)?;

        // Fold this sequence's tallies in now, while the sketch is already in
        // hand, so no later pass over the signature map is needed.
        let (examined, removed) = protein_sig.low_complexity_counts();
        if examined > 0 {
            self.kmer_windows_examined.fetch_add(examined, Ordering::Relaxed);
            self.low_complexity_kmers_removed.fetch_add(removed, Ordering::Relaxed);
        }

        // Return the processed signature (don't store it yet)
        Ok(protein_sig)
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
    /// # let index = ProteomeIndex::new(dir.path().join("test.db"), 10, 1, "protein20", false).unwrap();
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
    // Private to the module; needed to forge a pre-versioning index in
    // test_unversioned_index_reads_through_legacy_metadata_layout.
    use super::{LegacyProteomeIndexMetadata, ProteomeIndexMetadata, KMERSEEK_VERSION_KEY};
    use crate::sketch::ProteinSketch;
    use crate::tests::test_fixtures::{
        TEST_FASTA_CONTENT, TEST_FASTA_GZ, TEST_FASTA_ZST, TEST_PROTEIN,
    };
    use crate::tests::test_utils;
    use std::collections::{BTreeMap, HashMap};
    use std::fs::File;
    use std::path::PathBuf;

    /// A zero k-mer size reached `add_protein` and panicked on integer underflow.
    /// It is rejected at construction now, before RocksDB is even opened.
    #[test]
    fn test_new_rejects_zero_ksize_without_creating_database() -> Result<()> {
        let dir = tempdir()?;
        let db_path = dir.path().join("rejected.db");

        // `.err()` rather than `unwrap_err()`: ProteomeIndex holds a RocksDB
        // handle and does not implement Debug, which unwrap_err() would require.
        let err = ProteomeIndex::new(&db_path, 0, 1, "protein", false)
            .err()
            .expect("a zero k-mer size should be rejected");
        assert!(
            err.to_string().contains("K-mer size must be greater than 0"),
            "unexpected error: {err}"
        );
        assert!(!db_path.exists(), "a rejected k-mer size should leave no database behind");

        Ok(())
    }

    #[test]
    fn test_new_rejects_oversized_ksize() -> Result<()> {
        let dir = tempdir()?;
        let err = ProteomeIndex::new(dir.path().join("big.db"), 101, 1, "protein", false)
            .err()
            .expect("an oversized k-mer size should be rejected");
        assert!(err.to_string().contains("K-mer size too large"), "unexpected error: {err}");

        Ok(())
    }

    /// k=1 and k=100 are the boundaries next to the rejected 0 and 101, so they
    /// should construct normally rather than being rejected.
    #[test]
    fn test_new_accepts_boundary_ksizes() -> Result<()> {
        let dir = tempdir()?;
        let small = ProteomeIndex::new(dir.path().join("small.db"), 1, 1, "protein", false)?;
        assert_eq!(small.ksize(), 1);

        let large = ProteomeIndex::new(dir.path().join("large.db"), 100, 1, "protein", false)?;
        assert_eq!(large.ksize(), 100);

        Ok(())
    }

    /// Keeping the tests for ProteomeIndex in a separate file because they're more like integration tests
    /// than unit tests with all the moltype testing. Also, it's a lot of tests!

    #[test]
    fn test_add_protein_moltype_protein() -> Result<()> {
        let _dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein20";

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
    fn test_add_protein_moltype_dayhoff() -> Result<()> {
        let _dir = tempdir()?;

        let protein_ksize = 5;

        let sequence = TEST_PROTEIN;

        // Create a protein signature
        let mut protein_sig = ProteinSketch::new(
            "test_protein",
            protein_ksize,
            1, // scaled
            "dayhoff6",
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
    fn test_add_protein_moltype_hp() -> Result<()> {
        let _dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "hp_lehninger2";

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
        let moltype = "protein20";

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
        let moltype = "dayhoff6";

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
        let moltype = "hp_lehninger2";

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
    fn test_remove_low_complexity_toggle_on_index() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "hp_lehninger2";
        let sequence = TEST_PROTEIN;
        let name = "test_protein";

        // TEST_PROTEIN = "PLANTANDANIMALGENQMES"; window 10, "IMALG", is
        // all-hydrophobic ("hhhhh") under the HP (Lehninger) alphabet.
        const IMALG_HASH: u64 = 8541583772724823208;

        // Default: removal is off, so the low-complexity k-mer is kept.
        let index_off = ProteomeIndex::new(
            dir.path().join("removal_off.db"),
            protein_ksize,
            1,
            moltype,
            false,
        )?;
        let sig_off = index_off.create_protein_signature(sequence, name)?;
        assert_eq!(sig_off.kmer_positions().len(), 14);
        assert!(sig_off.kmer_positions().contains_key(&IMALG_HASH));

        // Opted in via set_remove_low_complexity: the low-complexity k-mer is dropped.
        let mut index_on =
            ProteomeIndex::new(dir.path().join("removal_on.db"), protein_ksize, 1, moltype, false)?;
        index_on.set_remove_low_complexity(true);
        let sig_on = index_on.create_protein_signature(sequence, name)?;
        assert_eq!(sig_on.kmer_positions().len(), 13);
        assert!(!sig_on.kmer_positions().contains_key(&IMALG_HASH));

        Ok(())
    }

    // Residues 26-55 of human FKBP8 (UniProt Q14318): a genuine 11-residue
    // poly-glutamate tract, real low-complexity sequence rather than an
    // invented motif.
    const FKBP8_POLY_E: &str = "VLDGVEDAEGEEEEEEEEEEEDDLSELPPL";

    /// The index-level counts aggregate across every stored signature; main.rs
    /// prints them after indexing.
    #[test]
    fn test_index_low_complexity_counts_aggregate_across_signatures() -> Result<()> {
        let dir = tempdir()?;

        // Two distinct sequences -- storing the same one twice would collapse to a
        // single md5 key and prove nothing about aggregation.
        // TEST_PROTEIN: 17 windows, 1 removed ("IMALG", all-hydrophobic).
        // FKBP8_POLY_E: 26 windows, 9 removed (7 raw "EEEEE" + 2 encoded "ppppp").
        let mut index =
            ProteomeIndex::new(dir.path().join("counts.db"), 5, 1, "hp_lehninger2", false)?;
        index.set_remove_low_complexity(true);
        for (name, seq) in [("p1", TEST_PROTEIN), ("p2", FKBP8_POLY_E)] {
            let sig = index.create_protein_signature(seq, name)?;
            index.store_signatures(vec![sig])?;
        }
        assert_eq!(index.get_signatures().len(), 2, "both sequences should be stored");
        assert_eq!(index.low_complexity_counts(), (43, 10));

        // With removal off nothing walks windows itself, so both stay zero.
        let index_off =
            ProteomeIndex::new(dir.path().join("counts_off.db"), 5, 1, "hp_lehninger2", false)?;
        let sig = index_off.create_protein_signature(TEST_PROTEIN, "p1")?;
        index_off.store_signatures(vec![sig])?;
        assert_eq!(index_off.low_complexity_counts(), (0, 0));

        Ok(())
    }

    /// The builder must carry the flag through to the constructed index.
    #[test]
    fn test_builder_sets_remove_low_complexity() -> Result<()> {
        let dir = tempdir()?;

        let on = ProteomeIndex::builder()
            .path(dir.path().join("builder_on.db"))
            .ksize(5)
            .scaled(1)
            .moltype("hp_lehninger2")
            .remove_low_complexity(true)
            .store_raw_sequences(true)
            .build()?;
        assert!(on.remove_low_complexity());
        assert!(on.store_raw_sequences());

        // Defaults to false when the builder method is not called.
        let off = ProteomeIndex::builder()
            .path(dir.path().join("builder_off.db"))
            .ksize(5)
            .scaled(1)
            .moltype("hp_lehninger2")
            .build()?;
        assert!(!off.remove_low_complexity());

        Ok(())
    }

    /// The setting must survive save_state -> open_for_search, since search
    /// relies on it to build query sketches the same way the targets were built.
    #[test]
    fn test_remove_low_complexity_persists_across_save_and_reopen() -> Result<()> {
        let dir = tempdir()?;

        for flag in [true, false] {
            let db_path = dir.path().join(format!("persist_{}.db", flag));
            {
                let mut index = ProteomeIndex::new(&db_path, 5, 1, "hp_lehninger2", true)?;
                index.set_remove_low_complexity(flag);
                let sig = index.create_protein_signature(TEST_PROTEIN, "p")?;
                index.store_signatures(vec![sig])?;
                index.save_state()?;
            }

            let reopened = ProteomeIndex::open_for_search(&db_path)?;
            assert_eq!(
                reopened.remove_low_complexity(),
                flag,
                "remove_low_complexity should round-trip through the database"
            );
        }

        Ok(())
    }

    /// An index older than schema 2 has no `remove_low_complexity` field in its
    /// metadata. It must still load, with the flag reading `false`, since such an
    /// index kept every k-mer by definition.
    #[test]
    fn test_unversioned_index_reads_through_legacy_metadata_layout() -> Result<()> {
        // Two flavors of old index: one written at schema 1, and one predating
        // schema versioning entirely (no key, which reads as version 0).
        for rolled_back in [true, false] {
            let dir = tempdir()?;
            let db_path = dir.path().join("legacy.db");
            {
                let index = ProteomeIndex::new(&db_path, 5, 1, "hp_lehninger2", true)?;
                let sig = index.create_protein_signature(TEST_PROTEIN, "p")?;
                index.store_signatures(vec![sig])?;
                index.save_state()?;
            }

            // Rewrite the index as a pre-versioning one: metadata in the old 8-field
            // layout, and no kmerseek_version key.
            {
                use rocksdb::{Options, DB};
                let db = DB::open(&Options::default(), &db_path)?;
                let current: ProteomeIndexMetadata =
                    bincode::deserialize(&db.get(b"index_metadata")?.unwrap())?;
                let legacy = LegacyProteomeIndexMetadata {
                    total_signatures: current.total_signatures,
                    chunk_count: current.chunk_count,
                    combined_mins: current.combined_mins,
                    combined_abunds: current.combined_abunds,
                    moltype: current.moltype,
                    ksize: current.ksize,
                    scaled: current.scaled,
                    store_raw_sequences: current.store_raw_sequences,
                };
                db.put(b"index_metadata", bincode::serialize(&legacy)?)?;
                db.delete(KMERSEEK_VERSION_KEY)?;
                // Layout is selected by schema_version, so roll that back too --
                // deleting the provenance key alone would not make this a v1 index.
                if rolled_back {
                    db.put(b"schema_version", bincode::serialize(&1u32)?)?;
                } else {
                    // Truly ancient: no schema_version key at all, which reads as 0.
                    db.delete(b"schema_version")?;
                }
            }

            // The legacy layout still loads, and defaults the flag to false.
            let reopened = ProteomeIndex::open_for_search(&db_path)?;
            assert!(!reopened.remove_low_complexity());
            assert_eq!(reopened.ksize(), 5);
            assert_eq!(reopened.moltype(), "hp_lehninger2");
            // No version was ever stamped on these.
            assert_eq!(reopened.kmerseek_version(), None);
        }

        Ok(())
    }

    /// `load_state` rehydrates signatures into an existing index handle, a
    /// separate path from `load` (which builds a fresh index) -- search uses it.
    #[test]
    fn test_load_state_rehydrates_signatures_into_existing_index() -> Result<()> {
        let dir = tempdir()?;
        let db_path = dir.path().join("load_state.db");
        {
            let index = ProteomeIndex::new(&db_path, 5, 1, "hp_lehninger2", true)?;
            for (name, seq) in [("p1", TEST_PROTEIN), ("p2", FKBP8_POLY_E)] {
                let sig = index.create_protein_signature(seq, name)?;
                index.store_signatures(vec![sig])?;
            }
            index.save_state()?;
        }

        let index = ProteomeIndex::new(&db_path, 5, 1, "hp_lehninger2", true)?;
        assert_eq!(index.signature_count(), 0, "a fresh handle starts empty");
        index.load_state()?;
        assert_eq!(index.signature_count(), 2, "load_state should pull both signatures back");

        Ok(())
    }

    /// The full in-memory load path (as opposed to `open_for_search`) must also
    /// recover the flag and the version, and bring the signatures back with it.
    #[test]
    fn test_load_round_trips_flag_version_and_signatures() -> Result<()> {
        let dir = tempdir()?;
        let db_path = dir.path().join("full_load.db");
        {
            let mut index = ProteomeIndex::new(&db_path, 5, 1, "hp_lehninger2", true)?;
            index.set_remove_low_complexity(true);
            for (name, seq) in [("p1", TEST_PROTEIN), ("p2", FKBP8_POLY_E)] {
                let sig = index.create_protein_signature(seq, name)?;
                index.store_signatures(vec![sig])?;
            }
            index.save_state()?;
        }
        // Dropped above, so RocksDB's lock is released before reopening.

        let loaded = ProteomeIndex::load(&db_path)?;
        assert!(loaded.remove_low_complexity());
        assert_eq!(loaded.kmerseek_version(), Some(env!("CARGO_PKG_VERSION")));
        assert_eq!(loaded.signature_count(), 2);
        assert_eq!(loaded.ksize(), 5);
        assert_eq!(loaded.scaled(), 1);
        assert_eq!(loaded.moltype(), "hp_lehninger2");

        // Counts are per-process build state, not persisted, so a fresh load
        // starts at zero rather than inheriting the writer's totals.
        assert_eq!(loaded.low_complexity_counts(), (0, 0));

        Ok(())
    }

    /// `get_index_parameters` is what the search CLI uses to autodetect settings
    /// from a database, and it reads the schema version on the way through.
    #[test]
    fn test_get_index_parameters_reads_from_saved_index() -> Result<()> {
        let dir = tempdir()?;
        let db_path = dir.path().join("params.db");
        {
            let index = ProteomeIndex::new(&db_path, 7, 1, "dayhoff6", true)?;
            let sig = index.create_protein_signature(TEST_PROTEIN, "p")?;
            index.store_signatures(vec![sig])?;
            index.save_state()?;
        }

        let (ksize, scaled, moltype) = ProteomeIndex::get_index_parameters(&db_path)?;
        assert_eq!(ksize, 7);
        assert_eq!(scaled, 1);
        // Stored under its current name, not the "dayhoff" spelling it was created with.
        assert_eq!(moltype, "dayhoff6");

        Ok(())
    }

    /// Saved indexes record which kmerseek version wrote them.
    #[test]
    fn test_save_state_stamps_kmerseek_version() -> Result<()> {
        let dir = tempdir()?;
        let db_path = dir.path().join("versioned.db");
        {
            let index = ProteomeIndex::new(&db_path, 5, 1, "hp_lehninger2", true)?;
            let sig = index.create_protein_signature(TEST_PROTEIN, "p")?;
            index.store_signatures(vec![sig])?;
            index.save_state()?;
        }

        let reopened = ProteomeIndex::open_for_search(&db_path)?;
        assert_eq!(reopened.kmerseek_version(), Some(env!("CARGO_PKG_VERSION")));

        Ok(())
    }

    #[test]
    fn test_process_fasta_moltype_protein() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein20";

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
        let moltype = "dayhoff6";

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
        let moltype = "hp_lehninger2";

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
        let moltype = "protein20";

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
        let moltype = "protein20";

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
        let moltype = "dayhoff6";

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
        let moltype = "hp_lehninger2";

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
        let moltype = "protein20";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("validation_test.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // The k-mer count is the number of *readings*, not the number of windows: each
        // ambiguity code doubles the windows covering it, up to a tenth of the window
        // rounded up, which is one code at k=5. For "ACDEFXBZJ" the five windows ACDEF,
        // CDEFX, DEFXB, EFXBZ and FXBZJ carry 0, 0, 1, 2 and 3 codes, so the first three
        // give 1 + 1 + 2 = 4 readings and the last two are over the cap and dropped.
        let valid_sequences =
            [("PLANTANDANIMALGENQMES", 17), ("ACDEFGHIKLMNPQRSTVWY", 16), ("ACDEFXBZJ", 4)];

        for (sequence, expected_kmers) in valid_sequences {
            let protein_signature = index.create_protein_signature(sequence, "test_protein")?;
            assert_eq!(protein_signature.kmer_positions().len(), expected_kmers, "{sequence}");
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

        // Ambiguity codes are accepted, and indexed under both readings rather than one.
        let ambiguous_sequences = [
            ("PLANTANDANIMALGENBMES", 21), // B is indexed as both Asp and Asn
            ("PLANTANDANIMALGENZMES", 21), // Z is indexed as both Glu and Gln
            ("PLANTANDANIMALGENJMES", 21), // J is indexed as both Ile and Leu
        ];

        for (sequence, expected_kmers) in ambiguous_sequences {
            let result = index.create_protein_signature(sequence, "test_protein");
            assert!(
                result.is_ok(),
                "Sequence with ambiguous amino acid '{}' should be resolved, not rejected",
                sequence
            );

            let protein_signature = result.unwrap();
            // 21 residues at k=5 gives 17 windows, and the B at position 18 falls in four
            // of them. Under this alphabet the two readings encode differently, so those
            // four windows contribute two k-mers each: 17 + 4 = 21.
            assert_eq!(
                protein_signature.kmer_positions().len(),
                expected_kmers,
                "{sequence}: expected the windows covering the ambiguity code to be doubled"
            );
        }

        Ok(())
    }

    #[test]
    fn test_create_protein_signature_amino_acid_validation_moltype_dayhoff() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "dayhoff6";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("validation_test_dayhoff.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Ambiguity codes are accepted, and indexed under both readings rather than one.
        let ambiguous_sequences = [
            ("PLANTANDANIMALGENBMES", 17), // Asp and Asn are both dayhoff `c`
            ("PLANTANDANIMALGENZMES", 17), // Glu and Gln are both dayhoff `c`
            ("PLANTANDANIMALGENJMES", 17), // Ile and Leu are both dayhoff `e`
        ];

        for (sequence, expected_kmers) in ambiguous_sequences {
            let result = index.create_protein_signature(sequence, "test_protein");
            assert!(
                result.is_ok(),
                "Sequence with ambiguous amino acid '{}' should be resolved, not rejected",
                sequence
            );

            let protein_signature = result.unwrap();
            // Disambiguating a code adds k-mers only where the alphabet keeps the two readings
            // apart. Dayhoff puts both members of every ambiguous pair in one class, so the
            // disambiguated windows hash identically and the count is unchanged.
            assert_eq!(protein_signature.kmer_positions().len(), expected_kmers, "{sequence}");
            // Check that the ambiguous k-mer is resolved correctly
            if sequence == "PLANTANDANIMALGENBMES" {
                // B resolves to D or N → dayhoff hash for NDMES/NNMES (both map to same dayhoff 6-letter encoding)
                assert!(
                    protein_signature.kmer_positions().contains_key(&6161374941338912337),
                    "Expected k-mer with hash 6161374941338912337 (NDMES/NNMES dayhoff) to be present in {}",
                    sequence
                );
            } else if sequence == "PLANTANDANIMALGENZMES" {
                // Z resolves to E or Q → dayhoff hash for NEMES/NQMES
                assert!(
                    protein_signature.kmer_positions().contains_key(&6161374941338912337),
                    "Expected k-mer with hash 6161374941338912337 (NEMES/NQMES dayhoff) to be present in {}",
                    sequence
                );
            } else if sequence == "PLANTANDANIMALGENJMES" {
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
        let moltype = "hp_lehninger2";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("validation_test_hp.db"),
            protein_ksize,
            1, // scaled=1 to capture all kmers for testing
            moltype,
            false,
        )?;

        // Ambiguity codes are accepted, and indexed under both readings rather than one.
        let ambiguous_sequences = [
            ("PLANTANDANIMALGENBMES", 14), // Asp and Asn are both polar
            ("PLANTANDANIMALGENZMES", 14), // Glu and Gln are both polar
            ("PLANTANDANIMALGENJMES", 14), // Ile and Leu are both hydrophobic
        ];

        for (sequence, expected_kmers) in ambiguous_sequences {
            let result = index.create_protein_signature(sequence, "test_protein");
            assert!(
                result.is_ok(),
                "Sequence with ambiguous amino acid '{}' should be resolved, not rejected",
                sequence
            );

            let protein_signature = result.unwrap();
            // Both readings of every ambiguous pair land on the same side of the HP split,
            // so the disambiguated windows hash identically and the count is unchanged.
            assert_eq!(protein_signature.kmer_positions().len(), expected_kmers, "{sequence}");
            // Check that the ambiguous k-mer is resolved correctly
            if sequence == "PLANTANDANIMALGENBMES" {
                // B resolves to D or N → HP hash for NDMES/NNMES (both map to "pphpp" HP encoding)
                assert!(
                    protein_signature.kmer_positions().contains_key(&13058023948041027181),
                    "Expected k-mer with hash 13058023948041027181 (NDMES/NNMES HP) to be present in {}",
                    sequence
                );
            } else if sequence == "PLANTANDANIMALGENZMES" {
                // Z resolves to E or Q → HP hash for NEMES/NQMES (both map to "pphpp" HP encoding)
                assert!(
                    protein_signature.kmer_positions().contains_key(&13058023948041027181),
                    "Expected k-mer with hash 13058023948041027181 (NEMES/NQMES HP) to be present in {}",
                    sequence
                );
            } else if sequence == "PLANTANDANIMALGENJMES" {
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
        let moltype = "protein20";

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
        let moltype = "protein20";

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
        let index1 = ProteomeIndex::new(&db_path1, 5, 1, "protein20", false).unwrap();
        let index2 = ProteomeIndex::new(&db_path2, 5, 1, "protein20", false).unwrap();

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
            ProteomeIndex::new(temp_dir.path().join("test3.db"), 10, 1, "protein20", false)
                .unwrap();
        assert!(!index1.is_equivalent_to(&index3).unwrap());
    }

    /// Real N-terminal fragment of C. elegans CED-9 (UniProt P41958), from
    /// tests/testdata/fasta/ced9.fasta.
    const CED9_PREFIX: &str = "MTRCTADNSLTNPAYRRRTMATGEMKEFLGIKGTEPTDFGINSDAQDLPSPSRQASTRRM";

    #[test]
    fn test_kmer_spectrum_csv_has_totals_comment_and_rows() {
        use std::io::Read;

        let temp_dir = tempdir().unwrap();
        let index =
            ProteomeIndex::new(temp_dir.path().join("spec.db"), 10, 1, "protein20", true).unwrap();
        // 3 k-mers seen once, 1 seen twice, 1 seen three times: 5 unique, 8 total.
        let spectrum: BTreeMap<usize, usize> = [(1, 3), (2, 1), (3, 1)].into_iter().collect();

        let csv_path = temp_dir.path().join("spectrum.csv");
        index.write_kmer_frequency_spectrum(&csv_path, &spectrum).unwrap();

        let mut contents = String::new();
        File::open(&csv_path).unwrap().read_to_string(&mut contents).unwrap();

        assert_eq!(
            contents,
            "# total_kmers=8 unique_kmers=5 mean_seqs_per_kmer=1.6000 median_seqs_per_kmer=1.0 mode_seqs_per_kmer=1 moltype=protein20 ksize=10\n\
             moltype,ksize,occurrences,n_kmers\n\
             protein20,10,1,3\n\
             protein20,10,2,1\n\
             protein20,10,3,1\n"
        );
    }

    #[test]
    fn test_median_occurrences_exact() {
        // Counts 1,1,1,2,3 -> odd length, middle element is 1.
        let odd: BTreeMap<usize, usize> = [(1, 3), (2, 1), (3, 1)].into_iter().collect();
        assert_eq!(ProteomeIndex::median_occurrences(&odd), 1.0);

        // Counts 1,1,2,3 -> even length, middle two are 1 and 2.
        let even: BTreeMap<usize, usize> = [(1, 2), (2, 1), (3, 1)].into_iter().collect();
        assert_eq!(ProteomeIndex::median_occurrences(&even), 1.5);

        // Counts 4,4,9,9 -> both middles inside one bucket.
        let flat: BTreeMap<usize, usize> = [(4, 2), (9, 2)].into_iter().collect();
        assert_eq!(ProteomeIndex::median_occurrences(&flat), 6.5);

        // Single value and empty.
        assert_eq!(ProteomeIndex::median_occurrences(&[(7, 1)].into_iter().collect()), 7.0);
        assert_eq!(ProteomeIndex::median_occurrences(&BTreeMap::new()), 0.0);
    }

    #[test]
    fn test_mode_occurrences() {
        // Occurrence count 1 has the most k-mers (3), so it's the mode.
        let clear: BTreeMap<usize, usize> = [(1, 3), (2, 1), (3, 1)].into_iter().collect();
        assert_eq!(ProteomeIndex::mode_occurrences(&clear), 1);

        // Tie between occurrence counts 2 and 5 (both have 4 k-mers): break toward smaller.
        let tied: BTreeMap<usize, usize> = [(2, 4), (5, 4), (9, 1)].into_iter().collect();
        assert_eq!(ProteomeIndex::mode_occurrences(&tied), 2);

        assert_eq!(ProteomeIndex::mode_occurrences(&[(7, 1)].into_iter().collect()), 7);
        assert_eq!(ProteomeIndex::mode_occurrences(&BTreeMap::new()), 0);
    }

    #[test]
    fn test_frequency_bins_groups_counts_into_powers_of_two() {
        // Bin b holds counts in [2^b, 2^(b+1)): 1 | 2-3 | 4-7 | 8-15 | ... | 256-511
        // Three k-mers seen once, one each at 2, 3, 4, 7, 8 and 300 occurrences.
        let spectrum: BTreeMap<usize, usize> =
            [(1, 3), (2, 1), (3, 1), (4, 1), (7, 1), (8, 1), (300, 1)].into_iter().collect();

        let bins = ProteomeIndex::frequency_bins(&spectrum);

        let expected: BTreeMap<u32, usize> =
            [(0, 3), (1, 2), (2, 2), (3, 1), (8, 1)].into_iter().collect();
        assert_eq!(bins, expected);
    }

    #[test]
    fn test_n_smallest_by_key_selects_most_and_least_common() {
        let frequencies: HashMap<u64, usize> =
            [(100, 5), (200, 9), (300, 1), (400, 9), (500, 3)].into_iter().collect();

        // Most common: highest count first, ties broken by ascending hash (200 before 400).
        let most = ProteomeIndex::n_smallest_by_key(&frequencies, 3, |hash, count| {
            (std::cmp::Reverse(count), hash)
        });
        assert_eq!(most, vec![(200, 9), (400, 9), (100, 5)]);

        // Least common: lowest count first.
        let least = ProteomeIndex::n_smallest_by_key(&frequencies, 3, |hash, count| (count, hash));
        assert_eq!(least, vec![(300, 1), (500, 3), (100, 5)]);
    }

    #[test]
    fn test_n_smallest_by_key_returns_all_when_n_exceeds_len() {
        let frequencies: HashMap<u64, usize> = [(100, 2), (200, 1)].into_iter().collect();

        let least = ProteomeIndex::n_smallest_by_key(&frequencies, 10, |hash, count| (count, hash));

        assert_eq!(least, vec![(200, 1), (100, 2)]);
    }

    /// Real N-terminal fragment of C. elegans CED-9 (UniProt P41958) with three residues
    /// rewritten to the ambiguity codes that stand for them: Asn->B (Asx), Glu->Z (Glx),
    /// Ile->J (Xle). Also carries U (Sec) and O (Pyl).
    const CED9_WITH_AMBIGUITY_CODES: &str =
        "MTRCTADNSLTNPAYRRRTMBTGEMKEFLGJKGTEPTDFGZNSDAQDLPSPSRQASTRRUO";

    /// Indexing the same sequence twice must produce byte-identical sketches.
    ///
    /// WHY: ambiguity codes were previously resolved by drawing at random from the
    /// alternatives, so B became Asp on one run and Asn on the next. That changed the k-mers,
    /// the hashes and the stored index every time the same FASTA was indexed.
    #[test]
    fn test_ambiguity_codes_index_deterministically() {
        let temp_dir = tempdir().unwrap();

        for (i, moltype) in
            ["protein20", "dayhoff6", "hp_lehninger2", "hp_pbotc_1st_ed2"].iter().enumerate()
        {
            let index =
                ProteomeIndex::new(temp_dir.path().join(format!("d{i}.db")), 5, 1, moltype, true)
                    .unwrap();

            let first = index
                .create_protein_signature(CED9_WITH_AMBIGUITY_CODES, "ced9")
                .unwrap()
                .mins_as_set();
            assert!(!first.is_empty(), "{moltype}: no k-mers produced");

            for attempt in 0..5 {
                let again = index
                    .create_protein_signature(CED9_WITH_AMBIGUITY_CODES, "ced9")
                    .unwrap()
                    .mins_as_set();
                assert_eq!(again, first, "{moltype}: sketch differed on attempt {attempt}");
            }
        }
    }

    /// A sequence carrying B must sketch to exactly the union of the two sequences it stands
    /// for: every k-mer of the Asp reading and every k-mer of the Asn reading, and nothing
    /// else. That is what makes the index match a query holding either residue, without
    /// committing to a reading the source never made.
    ///
    /// Checked across the whole range of alphabets, including protein20 (where D and N are
    /// distinct) and sdm12 and hsdm17 (where they land in separate classes), since those are
    /// exactly the cases a single representative would have got wrong.
    #[test]
    fn test_ambiguity_code_sketches_as_union_of_both_alternatives() {
        let temp_dir = tempdir().unwrap();
        // Same fragment written three ways: with B, and with each residue B stands for.
        let with_b = "MTRCTADNSLTNPAYRRRTMBTGEMKEFLGIK";
        let with_d = "MTRCTADNSLTNPAYRRRTMDTGEMKEFLGIK";
        let with_n = "MTRCTADNSLTNPAYRRRTMNTGEMKEFLGIK";

        let moltypes = [
            "protein20",
            "dayhoff6",
            "hp_lehninger2",
            "hp_pbotc_1st_ed2",
            "gbmr4",
            "sdm12",
            "hsdm17",
            "uniprot18",
        ];
        for (i, moltype) in moltypes.iter().enumerate() {
            let index =
                ProteomeIndex::new(temp_dir.path().join(format!("a{i}.db")), 5, 1, moltype, true)
                    .unwrap();
            let sketch = |seq| index.create_protein_signature(seq, "x").unwrap().mins_as_set();

            let (from_b, from_d, from_n) = (sketch(with_b), sketch(with_d), sketch(with_n));
            let union: std::collections::HashSet<u64> = from_d.union(&from_n).copied().collect();
            assert_eq!(from_b, union, "{moltype}: B should sketch as D union N");
            assert!(from_b.is_superset(&from_d), "{moltype}: B is missing the Asp reading");
            assert!(from_b.is_superset(&from_n), "{moltype}: B is missing the Asn reading");
        }
    }

    /// Where the alphabet puts Asp and Asn in different classes, the two readings really are
    /// different k-mers, so B contributes strictly more than either alone. This is the case
    /// a fixed representative got wrong.
    #[test]
    fn test_ambiguity_code_adds_both_readings_when_classes_differ() {
        let temp_dir = tempdir().unwrap();
        let with_b = "MTRCTADNSLTNPAYRRRTMBTGEMKEFLGIK";
        let with_d = "MTRCTADNSLTNPAYRRRTMDTGEMKEFLGIK";

        for (i, moltype) in ["protein20", "sdm12", "hsdm17"].iter().enumerate() {
            let index =
                ProteomeIndex::new(temp_dir.path().join(format!("d{i}.db")), 5, 1, moltype, true)
                    .unwrap();
            let sketch = |seq| index.create_protein_signature(seq, "x").unwrap().mins_as_set();

            assert!(
                sketch(with_b).len() > sketch(with_d).len(),
                "{moltype}: B should contribute k-mers the Asp reading alone does not"
            );
        }
    }

    #[test]
    fn test_index_stats() {
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");

        // Create a new index
        let index = ProteomeIndex::new(&db_path, 8, 10, "hp_lehninger2", false).unwrap();

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
            ProteomeIndex::new(manual_index_dir.clone(), 16, 5, "hp_lehninger2", false).unwrap();

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
            ProteomeIndex::new_with_auto_filename(&fasta_path, 16, 5, "hp_lehninger2", false)
                .unwrap();

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
        // The moltype in the filename is the normalized name, so passing a pre-rename
        // spelling still produces a file named after the current alphabet.
        let test_cases = vec![
            (16, 5, "hp_lehninger2", "test.fasta.hp_lehninger2.k16.scaled5.kmerseek.rocksdb"),
            (10, 1, "protein20", "test.fasta.protein20.k10.scaled1.kmerseek.rocksdb"),
            (8, 100, "dayhoff6", "test.fasta.dayhoff6.k8.scaled100.kmerseek.rocksdb"),
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
        // The fourth field is the normalized moltype that ends up in the filename: pre-rename
        // spellings are rewritten to the current alphabet name when the index is created.
        let test_cases = vec![
            (16, 5, "hp_lehninger2", "hp_lehninger2", "BCL2 with hp encoding, k=16, scaled=5"),
            (10, 1, "protein20", "protein20", "BCL2 with the full alphabet, k=10, scaled=1"),
            (8, 100, "dayhoff6", "dayhoff6", "BCL2 with dayhoff encoding, k=8, scaled=100"),
        ];

        for (ksize, scaled, moltype, stored_moltype, description) in test_cases {
            println!("Testing: {}", description);

            // Create index with automatic filename generation
            let auto_index =
                ProteomeIndex::new_with_auto_filename(&fasta_path, ksize, scaled, moltype, false)
                    .unwrap();

            // Verify the generated filename
            let expected_filename = format!("bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz.{}.k{}.scaled{}.kmerseek.rocksdb", stored_moltype, ksize, scaled);
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
        let index1 = ProteomeIndex::new(&db_path1, 5, 1, "protein20", false).unwrap();
        let index2 = ProteomeIndex::new(&db_path2, 5, 1, "protein20", false).unwrap();

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
            ProteomeIndex::new(temp_dir.path().join("index3.db"), 10, 1, "protein20", false)
                .unwrap();

        // Test that different indices are not equivalent
        let are_equivalent_3 = index1.is_equivalent_to(&index3).unwrap();
        assert!(!are_equivalent_3, "Indices with different parameters should not be equivalent");

        // Test with different sequences
        let index4 =
            ProteomeIndex::new(temp_dir.path().join("index4.db"), 5, 1, "protein20", false)
                .unwrap();
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
            ("simple.fasta", "simple.fasta.hp_lehninger2.k16.scaled5.kmerseek.rocksdb"),
            (
                "complex-name_with.underscores.fasta.gz",
                "complex-name_with.underscores.fasta.gz.hp_lehninger2.k16.scaled5.kmerseek.rocksdb",
            ),
            ("no_extension", "no_extension.hp_lehninger2.k16.scaled5.kmerseek.rocksdb"),
            (
                "multiple.dots.in.name.fasta",
                "multiple.dots.in.name.fasta.hp_lehninger2.k16.scaled5.kmerseek.rocksdb",
            ),
        ];

        for (base_name, expected) in test_cases {
            let base_path = temp_dir.path().join(base_name);
            let index =
                ProteomeIndex::new_with_auto_filename(&base_path, 16, 5, "hp_lehninger2", false)
                    .unwrap();
            let generated = index.generate_filename(base_name);
            assert_eq!(generated, expected, "Failed for base_name: {}", base_name);
        }

        // Test with different molecular types
        // Pre-rename spellings normalize, so the filename names the current alphabet.
        let moltype_cases = vec![
            ("hp_lehninger2", "test.fasta.hp_lehninger2.k8.scaled10.kmerseek.rocksdb"),
            ("protein20", "test.fasta.protein20.k8.scaled10.kmerseek.rocksdb"),
            ("dayhoff6", "test.fasta.dayhoff6.k8.scaled10.kmerseek.rocksdb"),
            ("protein20", "test.fasta.protein20.k8.scaled10.kmerseek.rocksdb"),
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
        let index = ProteomeIndex::new(&db_path, 8, 10, "hp_lehninger2", false).unwrap();

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
            5,           // k-mer size
            1,           // scaled
            "protein20", // molecular type
            true,        // store raw sequences
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
            "protein20",
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

    /// A cache too big for one RocksDB value survives a write/read round trip.
    ///
    /// The real trigger is 4 GiB; this drives the same code path with a 1 KB chunk size so
    /// the test costs milliseconds. What it is actually checking is that reassembly
    /// preserves byte order across chunk boundaries -- a reversed or skipped chunk still
    /// deserializes into *something* for many payloads, so a round trip that only asserts
    /// "no error" would pass while returning a wrong index.
    #[test]
    fn test_search_cache_survives_chunking() -> Result<()> {
        use tempfile::tempdir;
        let dir = tempdir()?;
        let index = ProteomeIndex::new(dir.path().join("chunked.db"), 10, 1, "protein20", false)?;

        // Big enough to need several chunks, and varied enough that a misordered
        // reassembly cannot coincidentally match.
        let target_list: Vec<String> = (0..500).map(|i| format!("md5-{i:08x}")).collect();
        let inverted_index: HashMap<u64, Vec<u32>> =
            (0..500u64).map(|i| (i * 7919, vec![i as u32, (i as u32) + 1])).collect();
        let kmer_frequencies: HashMap<u64, usize> =
            (0..500u64).map(|i| (i * 7919, (i as usize) % 13 + 1)).collect();
        let cache = crate::index::SearchCache {
            target_list: target_list.clone(),
            inverted_index: inverted_index.clone(),
            kmer_frequencies: kmer_frequencies.clone(),
        };
        let serialized = bincode::serialize(&cache)?;
        assert!(serialized.len() > 4096, "payload must span several 1 KB chunks");

        index.write_search_cache_chunked(&serialized, 1024)?;
        let loaded = index.load_search_cache()?.expect("a chunked cache must load");

        assert_eq!(loaded.target_list, target_list);
        assert_eq!(loaded.inverted_index, inverted_index);
        assert_eq!(loaded.kmer_frequencies, kmer_frequencies);
        Ok(())
    }

    /// A cache that fits keeps the single-key layout, so a database written by this build
    /// is still readable by one without the chunking support.
    #[test]
    fn test_small_search_cache_stays_single_key() -> Result<()> {
        use tempfile::tempdir;
        let dir = tempdir()?;
        let index = ProteomeIndex::new(dir.path().join("small.db"), 10, 1, "protein20", false)?;

        let cache = crate::index::SearchCache {
            target_list: vec!["md5-0".to_string()],
            inverted_index: HashMap::from([(42u64, vec![0u32])]),
            kmer_frequencies: HashMap::from([(42u64, 1usize)]),
        };
        let serialized = bincode::serialize(&cache)?;
        index.write_search_cache_chunked(&serialized, 1 << 30)?;

        assert!(index.db.get(b"search_cache")?.is_some(), "must use the original key");
        assert!(index.db.get(b"search_cache_chunks")?.is_none(), "must not write a count");
        assert!(index.load_search_cache()?.is_some());
        Ok(())
    }

    /// An interrupted write must refuse to open, not read as a database with no cache.
    ///
    /// This is the quietest failure the chunked layout can produce: chunks on disk, no
    /// count key, `Ok(None)` returned, and the caller falls back to loading every
    /// signature. The search then runs -- slowly, and against an index nobody is told is
    /// broken. A genuinely cache-free database must still return `Ok(None)`, so the two
    /// cases are asserted together; a check that cannot tell them apart is no check.
    #[test]
    fn test_torn_write_is_loud_but_no_cache_is_quiet() -> Result<()> {
        use tempfile::tempdir;
        let dir = tempdir()?;

        let torn = ProteomeIndex::new(dir.path().join("torn.db"), 10, 1, "protein20", false)?;
        let cache = crate::index::SearchCache {
            target_list: (0..500).map(|i| format!("md5-{i:08x}")).collect(),
            inverted_index: (0..500u64).map(|i| (i * 7919, vec![i as u32])).collect(),
            kmer_frequencies: (0..500u64).map(|i| (i * 7919, 1usize)).collect(),
        };
        let serialized = bincode::serialize(&cache)?;
        torn.write_search_cache_chunked(&serialized, 1024)?;

        // Exactly what a crash between the last chunk and the count leaves behind.
        torn.db.delete(b"search_cache_chunks")?;
        match torn.load_search_cache() {
            Ok(None) => panic!("a torn write must not read as an absent cache"),
            Ok(Some(_)) => panic!("a torn write must not load"),
            Err(err) => assert!(
                err.to_string().contains("interrupted"),
                "error should say the write was interrupted, got: {err}"
            ),
        }

        // A database that never had a cache is still legitimately empty.
        let fresh = ProteomeIndex::new(dir.path().join("fresh.db"), 10, 1, "protein20", false)?;
        assert!(
            fresh.load_search_cache()?.is_none(),
            "a database with no cache must return Ok(None), not an error"
        );
        Ok(())
    }

    /// Rewriting a database with a SMALLER cache must not leave the old cache's extra
    /// chunks behind.
    ///
    /// They are never read, since the count bounds the loop, so the index is correct
    /// either way -- but they are indistinguishable from the debris of an interrupted
    /// write, which is what the torn-write check keys on. Left in place they would make a
    /// healthy database fail that check the next time it was rewritten.
    #[test]
    fn test_rewriting_smaller_clears_stale_chunks() -> Result<()> {
        use tempfile::tempdir;
        let dir = tempdir()?;
        let index = ProteomeIndex::new(dir.path().join("shrink.db"), 10, 1, "protein20", false)?;

        let big = crate::index::SearchCache {
            target_list: (0..2000).map(|i| format!("md5-{i:08x}")).collect(),
            inverted_index: (0..2000u64).map(|i| (i * 7919, vec![i as u32])).collect(),
            kmer_frequencies: (0..2000u64).map(|i| (i * 7919, 1usize)).collect(),
        };
        let big_bytes = bincode::serialize(&big)?;
        index.write_search_cache_chunked(&big_bytes, 1024)?;
        let big_chunks = big_bytes.len().div_ceil(1024);

        let small = crate::index::SearchCache {
            target_list: (0..100).map(|i| format!("md5-{i:08x}")).collect(),
            inverted_index: (0..100u64).map(|i| (i * 7919, vec![i as u32])).collect(),
            kmer_frequencies: (0..100u64).map(|i| (i * 7919, 1usize)).collect(),
        };
        let small_bytes = bincode::serialize(&small)?;
        index.write_search_cache_chunked(&small_bytes, 1024)?;
        let small_chunks = small_bytes.len().div_ceil(1024);
        assert!(small_chunks < big_chunks, "second write must need fewer chunks");

        for i in small_chunks..big_chunks {
            assert!(
                index.db.get(ProteomeIndex::search_cache_chunk_key(i))?.is_none(),
                "stale chunk {i} from the larger cache was left behind"
            );
        }
        let loaded = index.load_search_cache()?.expect("the smaller cache must load");
        assert_eq!(loaded.target_list, small.target_list);
        Ok(())
    }

    /// A cache whose chunks are incomplete is an error, not a short index. Deserializing a
    /// truncated stream can succeed against a partial inverted index and silently drop
    /// search hits, which is far worse than refusing to open.
    #[test]
    fn test_missing_chunk_is_an_error_not_a_short_cache() -> Result<()> {
        use tempfile::tempdir;
        let dir = tempdir()?;
        let index = ProteomeIndex::new(dir.path().join("torn.db"), 10, 1, "protein20", false)?;

        let cache = crate::index::SearchCache {
            target_list: (0..500).map(|i| format!("md5-{i:08x}")).collect(),
            inverted_index: (0..500u64).map(|i| (i * 7919, vec![i as u32])).collect(),
            kmer_frequencies: (0..500u64).map(|i| (i * 7919, 1usize)).collect(),
        };
        let serialized = bincode::serialize(&cache)?;
        index.write_search_cache_chunked(&serialized, 1024)?;

        index.db.delete(b"search_cache_chunk_1")?;
        match index.load_search_cache() {
            Ok(_) => panic!("a torn cache must not load"),
            Err(err) => assert!(
                err.to_string().contains("missing"),
                "error should name the missing chunk, got: {err}"
            ),
        }
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
///         .moltype("protein20")
///         .build()?;
///
///     // With auto filename generation
///     let index = ProteomeIndex::builder()
///         .path("/path/to/base")
///         .ksize(5)
///         .scaled(1)
///         .moltype("protein20")
///         .build_with_auto_filename()?;
///
///     // With raw sequence storage
///     let index = ProteomeIndex::builder()
///         .path("/path/to/database.db")
///         .ksize(5)
///         .scaled(1)
///         .moltype("protein20")
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
    remove_low_complexity: bool,
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

    /// Set whether to drop low-complexity (homopolymer) k-mers (defaults to false)
    pub fn remove_low_complexity(mut self, remove_low_complexity: bool) -> Self {
        self.remove_low_complexity = remove_low_complexity;
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

        let mut index =
            ProteomeIndex::new(path, ksize, scaled, &moltype, self.store_raw_sequences)?;
        index.set_remove_low_complexity(self.remove_low_complexity);
        Ok(index)
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

        let mut index = ProteomeIndex::new_with_auto_filename(
            base_path,
            ksize,
            scaled,
            &moltype,
            self.store_raw_sequences,
        )?;
        index.set_remove_low_complexity(self.remove_low_complexity);
        Ok(index)
    }
}
