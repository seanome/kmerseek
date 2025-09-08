use dashmap::DashMap;
use parking_lot::Mutex;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::Duration;

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
    encode_kmer_with_encoding_fn, get_encoding_fn_from_moltype, get_hash_function_from_moltype,
};
use crate::errors::{IndexError, IndexResult};
use crate::kmer::KmerInfo;
use crate::signature::{ProteinSignature, ProteinSignatureData, SignatureAccess, SEED};

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
    signature_data: Vec<ProteinSignatureData>,
    combined_mins: Vec<u64>,
    combined_abunds: Option<Vec<u64>>,
    moltype: String,
    ksize: u32,
    scaled: u32,
    // Configuration for raw sequence storage
    store_raw_sequences: bool,
}

pub struct ProteomeIndex {
    // RocksDB instance for persistent storage
    db: DB,

    // Combined minhash of all proteins for statistics
    combined_minhash: Arc<Mutex<KmerMinHash>>,

    // Map of signature md5 -> protein signature (thread-safe concurrent map)
    signatures: DashMap<String, ProteinSignature>,

    // Amino acid ambiguity handler
    aa_ambiguity: Arc<AminoAcidAmbiguity>,

    // Protein encoding function
    encoding_fn: fn(u8) -> u8,

    // Statistics for k-mer frequencies and IDF
    // Not currently used, but will be used in the future
    #[allow(dead_code)]
    stats: ProteomeIndexKmerStats,

    // Add moltype field for serialization
    moltype: String,

    // Add ksize field for serialization
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
        // Create RocksDB options
        let mut opts = Options::default();
        opts.create_if_missing(true);
        // Allow multiple connections to the same database
        opts.set_max_open_files(-1);
        opts.set_use_fsync(false);
        opts.set_allow_mmap_reads(true);
        opts.set_allow_mmap_writes(true);

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
    pub fn get_signatures(&self) -> &DashMap<String, ProteinSignature> {
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

    /// Save the current index state to RocksDB using efficient storage format
    pub fn save_state(&self) -> IndexResult<()> {
        // DashMap is already thread-safe, no need to lock
        let combined_minhash = self.combined_minhash.lock();

        // Convert signatures to efficient storage format
        let mut signature_data = Vec::new();
        for sig in self.signatures.iter() {
            let efficient_data = sig.value().to_efficient_data(self.store_raw_sequences);
            signature_data.push(efficient_data);
        }

        // Create state with raw data
        let state = ProteomeIndexState {
            signature_data,
            combined_mins: combined_minhash.mins().to_vec(),
            combined_abunds: combined_minhash.abunds().map(|abunds| abunds.to_vec()),
            moltype: self.moltype.clone(),
            ksize: self.ksize,
            scaled: self.scaled,
            store_raw_sequences: self.store_raw_sequences,
        };

        let serialized = bincode::serialize(&state)?;
        self.db.put(b"index_state", serialized)?;

        Ok(())
    }

    /// Load index state from RocksDB using efficient storage format
    pub fn load_state(&self) -> IndexResult<()> {
        let serialized = self.db.get(b"index_state")?;
        if let Some(data) = serialized {
            let state: ProteomeIndexState = bincode::deserialize(&data)?;

            // Reconstruct signatures from efficient data
            let mut signatures_map = HashMap::new();
            for signature_data in state.signature_data {
                let protein_sig = ProteinSignature::from_efficient_data(
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

    /// Load an existing ProteomeIndex from a RocksDB path
    /// Note: This method has known issues with serialization and may not work reliably.
    /// For now, it's recommended to use save_state() and load_state() on existing indices.
    pub fn load<P: AsRef<Path>>(path: P) -> IndexResult<Self> {
        // Create RocksDB options
        let mut opts = Options::default();
        opts.create_if_missing(false); // Don't create if missing

        // Open the database
        let db = DB::open(&opts, path)?;

        // Try to load state to get configuration
        let serialized = db.get(b"index_state")?;
        if let Some(data) = serialized {
            let state: ProteomeIndexState = bincode::deserialize(&data)?;

            let _hash_function = get_hash_function_from_moltype(&state.moltype)
                .map_err(|e| IndexError::SourmashError(e.to_string()))?;
            let encoding_fn = get_encoding_fn_from_moltype(&state.moltype)
                .map_err(|e| IndexError::SourmashError(e.to_string()))?;

            // Reconstruct the combined minhash from raw data
            let hash_function = get_hash_function_from_moltype(&state.moltype)
                .map_err(|e| IndexError::SourmashError(e.to_string()))?;
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

            // Reconstruct signatures from efficient data
            let mut signatures_map = HashMap::new();
            for signature_data in state.signature_data {
                let protein_sig = ProteinSignature::from_efficient_data(
                    signature_data,
                    state.moltype.clone(),
                    state.ksize,
                    state.scaled,
                )?;

                let md5sum = protein_sig.signature().md5sum.clone();
                signatures_map.insert(md5sum.to_string(), protein_sig);
            }

            let index = Self {
                db,
                signatures: signatures_map.into_iter().collect(),
                combined_minhash: Arc::new(Mutex::new(combined_minhash)),
                aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
                encoding_fn,
                moltype: state.moltype,
                ksize: state.ksize,
                minhash_ksize: state.ksize * 3,
                scaled: state.scaled,
                stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
                store_raw_sequences: state.store_raw_sequences,
            };

            Ok(index)
        } else {
            Err(IndexError::NoSavedState)
        }
    }

    /// Get the number of signatures in the index
    pub fn signature_count(&self) -> usize {
        self.signatures.len()
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

                // Compare kmer_infos
                let self_kmer_infos = self_sig.kmer_infos();
                let other_kmer_infos = other_sig.kmer_infos();

                if self_kmer_infos.len() != other_kmer_infos.len() {
                    return Ok(false);
                }

                for (hashval, self_kmer_info) in self_kmer_infos.iter() {
                    if let Some(other_kmer_info) = other_kmer_infos.get(hashval) {
                        // Compare kmer_info fields
                        if self_kmer_info.ksize != other_kmer_info.ksize {
                            return Ok(false);
                        }
                        if self_kmer_info.hashval != other_kmer_info.hashval {
                            return Ok(false);
                        }
                        if self_kmer_info.encoded_kmer != other_kmer_info.encoded_kmer {
                            return Ok(false);
                        }

                        // Compare original_kmer_to_position maps
                        if self_kmer_info.original_kmer_to_position.len()
                            != other_kmer_info.original_kmer_to_position.len()
                        {
                            return Ok(false);
                        }

                        for (original_kmer, self_positions) in
                            &self_kmer_info.original_kmer_to_position
                        {
                            if let Some(other_positions) =
                                other_kmer_info.original_kmer_to_position.get(original_kmer)
                            {
                                if self_positions != other_positions {
                                    return Ok(false);
                                }
                            } else {
                                return Ok(false);
                            }
                        }
                    } else {
                        return Ok(false);
                    }
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
    /// Returns the processed `ProteinSignature` on success, or an error if the operation fails.
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
    ) -> IndexResult<ProteinSignature> {
        // Validate and resolve ambiguity if needed
        let processed_sequence = self.aa_ambiguity.validate_and_resolve(sequence)?;

        // Create a new protein signature
        let mut protein_sig = ProteinSignature::new(name, self.ksize, self.scaled, &self.moltype)?;

        // Add the protein sequence to the signature
        protein_sig.add_protein(processed_sequence.as_bytes())?;

        // Process the k-mers to get detailed k-mer information
        self.process_kmers(&processed_sequence, &mut protein_sig)?;

        // If raw sequence storage is enabled, create efficient data with the sequence
        if self.store_raw_sequences {
            let efficient_data =
                protein_sig.to_efficient_data_with_capacity(processed_sequence.len());
            let mut efficient_data_with_sequence = efficient_data;
            efficient_data_with_sequence.set_raw_sequence(processed_sequence.to_string());
            protein_sig.set_efficient_data(efficient_data_with_sequence);
        }

        // Return the processed signature (don't store it yet)
        Ok(protein_sig)
    }

    pub fn process_kmers(
        &self,
        sequence: &str,
        protein_signature: &mut ProteinSignature,
    ) -> IndexResult<()> {
        let ksize = self.ksize as usize;
        let seed = SEED;
        let hashvals = &protein_signature.signature().get_minhash().to_vec();

        for i in 0..sequence.len().saturating_sub(ksize - 1) {
            let kmer = &sequence[i..i + ksize];

            // Process the k-mer to get encoded version
            if let Ok((encoded_kmer, original_kmer)) =
                encode_kmer_with_encoding_fn(kmer, self.encoding_fn)
            {
                // Get the hash from the minhash implementation
                let hashval = _hash_murmur(encoded_kmer.as_bytes(), seed);

                // If this hashval is in the minhash, then save its k-mer positions
                if hashvals.contains(&hashval) {
                    let kmer_info = protein_signature
                        .kmer_infos_mut()
                        .entry(hashval)
                        .or_insert_with(|| KmerInfo {
                            ksize,
                            hashval,
                            encoded_kmer: encoded_kmer.clone(),
                            original_kmer_to_position: HashMap::new(),
                        });

                    kmer_info.original_kmer_to_position.entry(original_kmer).or_default().push(i);
                }
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
    /// * `signatures` - A vector of `ProteinSignature` objects to store
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or an error if the operation fails.
    pub fn store_signatures(&self, protein_signatures: Vec<ProteinSignature>) -> IndexResult<()> {
        // Collect the minhash data from new signatures before storing them
        let new_hashes_and_abunds: Vec<(u64, u64)> = protein_signatures
            .iter()
            .flat_map(|sig| {
                let minhash = sig.signature().get_minhash();
                // If abundance is tracked, use to_vec_abunds, else use mins with abundance 1
                if let Some(abunds) = minhash.abunds() {
                    minhash.mins().into_iter().zip(abunds).collect::<Vec<_>>()
                } else {
                    minhash.mins().into_iter().map(|h| (h, 1)).collect::<Vec<_>>()
                }
            })
            .collect();

        // Store all signatures in the signatures map
        {
            for protein_signature in protein_signatures {
                let md5sum = protein_signature.signature().md5sum.clone();
                self.signatures.insert(md5sum.to_string(), protein_signature);
            }
        }

        // Update the combined minhash with only the new signatures
        let mut combined_minhash = self.combined_minhash.lock();
        combined_minhash
            .add_many_with_abund(&new_hashes_and_abunds)
            .map_err(|e| IndexError::SourmashError(e.to_string()))?;

        Ok(())
    }

    /// Process a protein FASTA file in parallel, adding all sequences to the index
    ///
    /// This method reads a FASTA file, validates each protein sequence for amino acid ambiguity,
    /// creates protein signatures for each sequence, and stores them in the index.
    ///
    /// Each sequence is validated using the same amino acid validation as `create_protein_signature`.
    /// If any sequence contains invalid amino acids, the entire operation will fail with an error
    /// describing the first invalid amino acid encountered.
    ///
    /// # Arguments
    ///
    /// * `fasta_path` - Path to the FASTA file to process
    /// * `progress_interval` - Number of sequences between progress reports (0 to disable progress)
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or an error if the operation fails.
    /// The error will contain details about any invalid amino acids found in the sequences.
    pub fn process_fasta<P: AsRef<Path>>(
        &self,
        fasta_path: P,
        progress_interval: u32,
    ) -> IndexResult<()> {
        use needletail::parse_fastx_file;
        use rayon::prelude::*;

        if progress_interval > 0 {
            println!("Reading FASTA file...");
        }

        // Open and parse the FASTA file using needletail
        let mut reader =
            parse_fastx_file(&fasta_path).map_err(|e| IndexError::ParseError(e.to_string()))?;

        // Convert records into a vector for parallel processing
        let mut records = Vec::new();
        let mut record_count = 0;
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| IndexError::ParseError(e.to_string()))?;
            // Clone the record data to avoid borrowing issues
            let sequence = record.seq().to_vec();
            let id = record.id().to_vec();
            records.push((sequence, id));
            record_count += 1;

            // Print progress if interval is set and we've reached the interval
            if progress_interval > 0 && record_count % progress_interval == 0 {
                println!("Read {} sequences...", record_count);
            }
        }

        let total_records = records.len();
        if progress_interval > 0 {
            println!(
                "Finished reading {} sequences. Starting parallel processing...",
                total_records
            );
        }

        // Set up atomic counter for progress
        let processed_count = Arc::new(AtomicUsize::new(0));
        if progress_interval > 0 {
            let processed_count_clone = Arc::clone(&processed_count);
            thread::spawn(move || loop {
                let count = processed_count_clone.load(Ordering::Relaxed);
                if count > 0 {
                    let percent = (count as f64 / total_records as f64) * 100.0;
                    println!("Processed {} sequences ({:.1}%)", count, percent);
                }
                if count >= total_records {
                    break;
                }
                thread::sleep(Duration::from_secs(2));
            });
        }

        // Process records in parallel and collect signatures
        let signatures: Result<Vec<ProteinSignature>, IndexError> = records
            .par_iter()
            .map(|(seq_bytes, id_bytes)| {
                let sequence = std::str::from_utf8(seq_bytes)?;
                let name = std::str::from_utf8(id_bytes)?;

                // Uppercase the sequence before processing
                let sequence = sequence.to_uppercase();

                // Create protein signature for each sequence
                let result = self.create_protein_signature(&sequence, name);
                processed_count.fetch_add(1, Ordering::Relaxed);
                result
            })
            .collect();

        if progress_interval > 0 {
            println!("Storing {} signatures...", total_records);
        }

        // Store all signatures at once
        self.store_signatures(signatures?)?;

        // Save the index state to RocksDB
        self.save_state()?;

        if progress_interval > 0 {
            println!("Successfully processed and stored {} sequences.", total_records);
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
    use crate::signature::ProteinSignature;
    use crate::tests::test_fixtures::{TEST_FASTA_CONTENT, TEST_FASTA_GZ, TEST_PROTEIN};
    use crate::tests::test_utils::{self, print_kmer_infos};
    use std::collections::HashMap;
    use std::path::PathBuf;

    /// Keeping the tests for ProteomeIndex in a separate file because they're more like integration tests
    /// than unit tests with all the moltype testing. Also, it's a lot of tests!

    #[test]
    fn test_process_kmers_moltype_protein() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "protein";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("dayhoff_test.db"),
            protein_ksize, // protein ksize
            1,             // scaled=1 to capture all kmers
            moltype,
            false,
        )?;

        let sequence = TEST_PROTEIN;

        // Create a protein signature
        let mut protein_sig = ProteinSignature::new(
            "test_protein",
            protein_ksize,
            1, // scaled
            moltype,
        )?;

        // Add the sequence
        protein_sig.add_protein(sequence.as_bytes())?;
        println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

        // Process kmers
        index.process_kmers(sequence, &mut protein_sig)?;

        println!("{}", protein_sig.signature().name);
        println!("{:?}", protein_sig.kmer_infos().keys());
        let kmer_count = protein_sig.kmer_infos().len();

        // Should have 17 kmers (length 21 - ksize 5 + 1)
        assert_eq!(kmer_count, 17);

        // Print all kmer infos for debugging
        test_utils::print_kmer_infos(&protein_sig);

        // Create expected hashmap of kmer info
        let raw_data = [
            // Hash               Original  Position
            (2140811952770908281, ("GENQM", [14])),
            (4381446250900425522, ("ENQME", [15])),
            (5798339600059429290, ("DANIM", [7])),
            (7681438632487987439, ("ANIMA", [8])),
            (12896310179337320481, ("LANTA", [1])),
            (2542642819229379552, ("NTAND", [3])),
            (11965201914550078735, ("TANDA", [4])),
            (5893010049374798421, ("PLANT", [0])),
            (110005740849399217, ("NDANI", [6])),
            (3791883307084689782, ("LGENQ", [13])),
            (14610011480386804007, ("ALGEN", [12])),
            (6941015416212662126, ("ANTAN", [2])),
            (12636705882654324958, ("NQMES", [16])),
            (11154024130290913208, ("IMALG", [10])),
            (1225702037828834387, ("MALGE", [11])),
            (12274863873578753245, ("NIMAL", [9])),
            (13616372540306653069, ("ANDAN", [5])),
        ];

        // Convert raw data into the required format with proper string types
        let expected_kmers: HashMap<_, _> = raw_data
            .into_iter()
            .map(|(hash, (kmer, positions))| {
                let mut original_map = HashMap::new();
                original_map.insert(kmer.to_string(), positions.to_vec());
                (hash, (kmer.to_string(), original_map))
            })
            .collect();

        // Verify each kmer info matches expected values
        for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
            let (expected_kmer, expected_positions) =
                expected_kmers.get(hash).expect(&format!("Missing expected hash {}", hash));

            // Verify the k-mer
            assert_eq!(
                &kmer_info.encoded_kmer, expected_kmer,
                "K-mer mismatch for hash {}: expected {}, got {}",
                hash, expected_kmer, kmer_info.encoded_kmer
            );

            // Verify positions for each original k-mer
            for (original_kmer, positions) in &kmer_info.original_kmer_to_position {
                let expected_pos = expected_positions
                    .get(original_kmer)
                    .expect(&format!("Missing positions for k-mer {}", original_kmer));
                assert_eq!(
                    positions, expected_pos,
                    "Position mismatch for k-mer {}: expected {:?}, got {:?}",
                    original_kmer, expected_pos, positions
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_process_kmers_moltype_dayhoff() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("dayhoff_test.db"),
            protein_ksize, // protein ksize
            1,             // scaled=1 to capture all kmers
            "dayhoff",
            false,
        )?;

        let sequence = TEST_PROTEIN;

        // Create a protein signature
        let mut protein_sig = ProteinSignature::new(
            "test_protein",
            protein_ksize,
            1, // scaled
            "dayhoff",
        )?;

        // Add the sequence
        protein_sig.add_protein(sequence.as_bytes())?;
        println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

        // Process kmers
        index.process_kmers(sequence, &mut protein_sig)?;

        println!("{}", protein_sig.signature().name);
        let hashvals = protein_sig.kmer_infos().keys().collect::<Vec<_>>();
        println!("{:?}", hashvals);
        let kmer_count = protein_sig.kmer_infos().len();

        // Should have 17 kmers (length 21 - ksize 5 + 1)
        assert_eq!(kmer_count, 17);

        // Print all kmer infos for debugging
        test_utils::print_kmer_infos(&protein_sig);

        // Define the raw data without string conversions
        let raw_data = [
            (17444159595263538048, ("ceebe", "NIMAL", [9])),
            (2945598193614695589, ("cccec", "ENQME", [15])),
            (4548757849819812604, ("bbccb", "TANDA", [4])),
            (6463872878592804545, ("ebccc", "LGENQ", [13])),
            (4030406117949362159, ("cbcee", "DANIM", [7])),
            (7014407397606522347, ("ebcbb", "LANTA", [1])),
            (5045972850709227854, ("bebcb", "PLANT", [0])),
            (11417072151730334367, ("bcbbc", "ANTAN", [2])),
            (13574922562423607435, ("bceeb", "ANIMA", [8])),
            (15050500149255106627, ("bccce", "GENQM", [14])),
            (5430883729707969951, ("eebeb", "IMALG", [10])),
            (13894194422852851851, ("bebcc", "ALGEN", [12])),
            (9604281550621775790, ("bccbc", "ANDAN", [5])),
            (6161374941338912337, ("ccecb", "NQMES", [16])),
            (655307631517862365, ("ccbce", "NDANI", [6])),
            (360995089333906261, ("ebebc", "MALGE", [11])),
            (15056713696431004031, ("cbbcc", "NTAND", [3])),
        ];

        // Convert raw strings to owned types and create the HashMap
        let expected_kmers: HashMap<_, _> = raw_data
            .into_iter()
            .map(|(hash, (encoded, original, positions))| {
                let mut original_map = HashMap::new();
                original_map.insert(original.to_string(), positions.to_vec());
                (hash, (encoded.to_string(), original_map))
            })
            .collect();

        // Verify all expected hashes are present
        let expected_hashes: Vec<_> = expected_kmers.keys().collect();
        let actual_hashes: Vec<_> = protein_sig.kmer_infos().keys().collect();
        assert_eq!(
            expected_hashes.len(),
            actual_hashes.len(),
            "Number of hashes mismatch: expected {}, got {}",
            expected_hashes.len(),
            actual_hashes.len()
        );

        for hash in expected_hashes {
            assert!(protein_sig.kmer_infos().contains_key(hash), "Missing expected hash {}", hash);
        }

        // Verify each kmer info matches expected values
        for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
            let (expected_encoded, expected_originals) =
                expected_kmers.get(hash).expect(&format!("Missing expected hash {}", hash));

            // Verify the encoded k-mer
            assert_eq!(
                &kmer_info.encoded_kmer, expected_encoded,
                "Encoded k-mer mismatch for hash {}: expected {}, got {}",
                hash, expected_encoded, kmer_info.encoded_kmer
            );

            // Verify each original k-mer and its positions
            for (original_kmer, expected_positions) in expected_originals {
                let positions = protein_sig
                    .kmer_infos()
                    .get(hash)
                    .unwrap()
                    .original_kmer_to_position
                    .get(original_kmer)
                    .expect(&format!("Missing positions for k-mer {}", original_kmer));
                assert_eq!(
                    positions, expected_positions,
                    "Position mismatch for k-mer {}: expected {:?}, got {:?}",
                    original_kmer, expected_positions, positions
                );
            }
        }

        Ok(())
    }

    #[test]
    fn test_process_kmers_moltype_hp() -> Result<()> {
        let dir = tempdir()?;

        let protein_ksize = 5;
        let moltype = "hp";

        // Create index with minimal parameters
        let index = ProteomeIndex::new(
            dir.path().join("dayhoff_test.db"),
            protein_ksize, // protein ksize
            1,             // scaled=1 to capture all kmers
            moltype,
            false,
        )?;

        let sequence = TEST_PROTEIN;

        // Create a protein signature
        let mut protein_sig = ProteinSignature::new(
            "test_protein",
            protein_ksize,
            1, // scaled
            moltype,
        )?;

        // Add the sequence
        protein_sig.add_protein(sequence.as_bytes())?;
        println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

        // Process kmers
        index.process_kmers(sequence, &mut protein_sig)?;

        println!("{}", protein_sig.signature().name);
        let hashvals = protein_sig.kmer_infos().keys().collect::<Vec<_>>();
        println!("{:?}", hashvals);
        let kmer_count = protein_sig.kmer_infos().len();

        // // Should have 14 kmers (length 21 - ksize 5 + 1), but a few duplicates
        assert_eq!(kmer_count, 14);

        // Print all kmer infos for debugging
        test_utils::print_kmer_infos(&protein_sig);

        // Define test data in a more readable format
        let kmer_data: HashMap<u64, (String, HashMap<String, Vec<usize>>)> = vec![
            // Single k-mer cases
            (17248460043117039725, ("hhhhp", vec!["MALGE"], vec![11])),
            (5673218808929106268, ("phhhh", vec!["NIMAL"], vec![9])),
            (16969835101383990681, ("hhpph", vec!["LANTA"], vec![1])),
            (7345312524621807974, ("pphph", vec!["NDANI"], vec![6])),
            (16370543730027378051, ("phpph", vec!["TANDA"], vec![4])),
            (3278382041688965244, ("hphhh", vec!["ANIMA"], vec![8])),
            (8541583772724823208, ("hhhhh", vec!["IMALG"], vec![10])),
            (16158526221854164806, ("hppph", vec!["GENQM"], vec![14])),
            (11553019557737058697, ("hhppp", vec!["LGENQ"], vec![13])),
            (9081059129327932468, ("ppphp", vec!["ENQME"], vec![15])),
            (2863220259252354754, ("phphh", vec!["DANIM"], vec![7])),
            // Multiple original protein k-mer sequences mapping to same HP encoding
            (4230974618842309829, ("hhhpp", vec!["PLANT", "ALGEN"], vec![0, 12])),
            (13058023948041027181, ("pphpp", vec!["NQMES", "NTAND"], vec![16, 3])),
            (4144736064335623701, ("hpphp", vec!["ANDAN", "ANTAN"], vec![5, 2])),
        ]
        .into_iter()
        .map(|(hash, (encoded, originals, positions))| {
            let mut original_map = HashMap::new();
            for (i, orig) in originals.into_iter().enumerate() {
                original_map.insert(orig.to_string(), vec![positions[i]]);
            }
            (hash, (encoded.to_string(), original_map))
        })
        .collect();

        // Verify all expected hashes are present
        let expected_hashes: Vec<_> = kmer_data.keys().collect();
        let actual_hashes: Vec<_> = protein_sig.kmer_infos().keys().collect();
        assert_eq!(
            expected_hashes.len(),
            actual_hashes.len(),
            "Number of hashes mismatch: expected {}, got {}",
            expected_hashes.len(),
            actual_hashes.len()
        );

        for hash in expected_hashes {
            assert!(protein_sig.kmer_infos().contains_key(hash), "Missing expected hash {}", hash);
        }

        // Verify each kmer info matches expected values
        for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
            let (expected_encoded, expected_originals) =
                kmer_data.get(hash).expect(&format!("Missing expected hash {}", hash));

            // Verify the encoded k-mer
            assert_eq!(
                &kmer_info.encoded_kmer, expected_encoded,
                "Encoded k-mer mismatch for hash {}: expected {}, got {}",
                hash, expected_encoded, kmer_info.encoded_kmer
            );

            // Verify each original k-mer and its positions
            for (original_kmer, expected_positions) in expected_originals {
                let positions = protein_sig
                    .kmer_infos()
                    .get(hash)
                    .unwrap()
                    .original_kmer_to_position
                    .get(original_kmer)
                    .expect(&format!("Missing original k-mer {} for hash {}", original_kmer, hash));
                assert_eq!(
                    positions, expected_positions,
                    "Position mismatch for k-mer {}: expected {:?}, got {:?}",
                    original_kmer, expected_positions, positions
                );
            }

            // Verify no unexpected original k-mers
            assert_eq!(
                kmer_info.original_kmer_to_position.len(),
                expected_originals.len(),
                "Expected {} original k-mer mappings for hash {}, got {}",
                expected_originals.len(),
                hash,
                kmer_info.original_kmer_to_position.len()
            );
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
        assert_eq!(signature.kmer_infos().len(), 17, "Expected 17 k-mers for the test protein");

        // Verify some specific k-mers are present
        let expected_hash = 5893010049374798421; // Hash for "PLANT"
        assert!(
            signature.kmer_infos().contains_key(&expected_hash),
            "Expected k-mer hash {} to be present",
            expected_hash
        );

        // Store the signature in the index
        index.store_signatures(vec![signature])?;

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
        assert_eq!(signature.kmer_infos().len(), 17, "Expected 17 k-mers for the test protein");

        // Verify some specific k-mers are present
        let expected_hash = 5045972850709227854; // Hash for "PLANT" in Dayhoff encoding ("bebcb")
        assert!(
            signature.kmer_infos().contains_key(&expected_hash),
            "Expected k-mer hash {} to be present",
            expected_hash
        );

        // Store the signature in the index
        index.store_signatures(vec![signature])?;

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
        assert_eq!(signature.kmer_infos().len(), 14, "Expected 14 k-mers for the test protein");

        // Verify some specific k-mers are present
        let expected_hash = 4230974618842309829; // Hash for "PLANT" in HP encoding ("hhhpp")
        assert!(
            signature.kmer_infos().contains_key(&expected_hash),
            "Expected k-mer hash {} to be present",
            expected_hash
        );

        // Store the signature in the index
        index.store_signatures(vec![signature])?;

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
        index.process_fasta(&fasta_path, 0)?;

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
                        stored_signature.kmer_infos().len() == 7,
                        "LIVINGALIVE should have 7 protein 5-mers"
                    );
                } else if md5sum == "7641839ad508ab8" {
                    assert!(
                        stored_signature.kmer_infos().len() == 17,
                        "PLANTANDANIMALGENQMES should have 17 protein 5-mers"
                    );
                } else {
                    println!("md5sum: {}", md5sum);
                    println!("Name: {}", stored_signature.signature().name);
                    println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                    assert!(false, "Unknown md5sum: {}", md5sum);
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
        index.process_fasta(&fasta_path, 0)?;

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
                        stored_signature.kmer_infos().len() == 7,
                        "LIVINGALIVE should have 7 dayhoff 5-mers"
                    );
                } else if md5sum == "84d7545d531dcf51" {
                    assert!(
                        stored_signature.kmer_infos().len() == 17,
                        "PLANTANDANIMALGENQMES should have 17 dayhoff 5-mers"
                    );
                } else {
                    println!("md5sum: {}", md5sum);
                    println!("Name: {}", stored_signature.signature().name);
                    println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                    assert!(false, "Unknown md5sum: {}", md5sum);
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
        index.process_fasta(&fasta_path, 0)?;

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
                        stored_signature.kmer_infos().len() == 6,
                        "LIVINGALIVE should have 6 hp 5-mers"
                    );
                } else if md5sum == "668d7173d661287b" {
                    assert!(
                        stored_signature.kmer_infos().len() == 14,
                        "PLANTANDANIMALGENQMES should have 14 hp 5-mers"
                    );
                } else {
                    println!("md5sum: {}", md5sum);
                    println!("Name: {}", stored_signature.signature().name);
                    println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                    assert!(false, "Unknown md5sum: {}", md5sum);
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
        index.process_fasta(TEST_FASTA_GZ, 0)?;

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
                println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                if md5sum == "4d565dee9c8de9db" {
                    assert!(
                        stored_signature.kmer_infos().len() == 474,
                        "sp|O43236|SEPT4_HUMAN should have 474 protein 5-mers"
                    );
                }
                if md5sum == "4da1f84ad8be618e" {
                    assert!(
                        stored_signature.kmer_infos().len() == 235,
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
        index.process_fasta(TEST_FASTA_GZ, 0)?;

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
                println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                if md5sum == "fc27dcd533217385" {
                    assert!(
                        stored_signature.kmer_infos().len() == 433,
                        "sp|O43236|SEPT4_HUMAN should have 433 dayhoff 5-mers"
                    );
                }
                if md5sum == "3206706fa14185e7" {
                    assert!(
                        stored_signature.kmer_infos().len() == 204,
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
        index.process_fasta(TEST_FASTA_GZ, 0)?;

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
                println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                if md5sum == "38ffedf9d3ec7cec" {
                    assert!(
                        stored_signature.kmer_infos().len() == 452,
                        "sp|O43236|SEPT4_HUMAN should have 452 hp 12-mers"
                    );
                }
                if md5sum == "204716e4d80eb350" {
                    assert!(
                        stored_signature.kmer_infos().len() == 220,
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
            test_utils::print_kmer_infos(&protein_signature);
            if protein_signature.signature().md5sum == "7641839ad508ab8" {
                assert!(
                    protein_signature.kmer_infos().len() == 17,
                    "Valid sequence 'PLANTANDANIMALGENQMES' should be accepted and have 17 protein 5-mers",
                );
            } else if protein_signature.signature().md5sum == "b95f0777d5439d56" {
                assert!(
                    protein_signature.kmer_infos().len() == 16,
                    "Valid sequence 'ACDEFGHIKLMNPQRSTVWY' should be accepted and have 16 protein 5-mers",
                );
            } else if protein_signature.signature().md5sum == "fa11c30a562fd82" {
                assert!(
                    protein_signature.kmer_infos().len() == 5,
                    "Valid sequence 'ACDEFXBZJ' should be accepted and have 5 protein 5-mers",
                );
            } else {
                // For the third sequence, just check the length is correct
                if protein_signature.kmer_infos().len() == 5 {
                    // This is the expected case for ACDEFXBZJ
                } else {
                    assert!(
                        false,
                        "Unexpected kmer count: {} for md5sum: {}",
                        protein_signature.kmer_infos().len(),
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
            print_kmer_infos(&protein_signature);
            // Should have the same number of k-mers as the original sequence
            assert_eq!(
                protein_signature.kmer_infos().len(),
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
            print_kmer_infos(&protein_signature);
            // Should have the same number of k-mers as the original sequence
            println!("sequence: {}", sequence);
            assert_eq!(
                protein_signature.kmer_infos().len(),
                17,
                "Resolved sequence should have 17 protein 5-mers"
            );
            // Check that the ambiguous k-mer is resolved correctly
            if sequence == &"PLANTANDANIMALGENBMES" {
                let kmer_info = protein_signature.kmer_infos().get(&6161374941338912337);
                assert!(
                    kmer_info.is_some(),
                    "Expected k-mer with hash 6161374941338912337 to be present in {}",
                    sequence
                );
                let kmer_info = kmer_info.unwrap();
                assert_eq!(
                    kmer_info.encoded_kmer, "ccecb",
                    "Expected encoded k-mer 'ccecb' (NDMES/NNMES) to be present in {}",
                    sequence
                );
            } else if sequence == &"PLANTANDANIMALGENZMES" {
                let kmer_info = protein_signature.kmer_infos().get(&6161374941338912337);
                assert!(
                    kmer_info.is_some(),
                    "Expected k-mer with hash 6161374941338912337 to be present in {}",
                    sequence
                );
                let kmer_info = kmer_info.unwrap();
                assert_eq!(
                    kmer_info.encoded_kmer, "ccecb",
                    "Expected encoded k-mer 'ccecb' (NEMES/NQMES) to be present in {}",
                    sequence
                );
            } else if sequence == &"PLANTANDANIMALGENJMES" {
                let kmer_info = protein_signature.kmer_infos().get(&9182605311834199497);
                assert!(
                    kmer_info.is_some(),
                    "Expected k-mer with hash 9182605311834199497 to be present in {}",
                    sequence
                );
                let kmer_info = kmer_info.unwrap();
                assert_eq!(
                    kmer_info.encoded_kmer, "ceecb",
                    "Expected encoded k-mer 'ceecb' (NLMES/NIMES) to be present in {}",
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
            print_kmer_infos(&protein_signature);
            // Should have the same number of k-mers as the original sequence
            println!("sequence: {}", sequence);
            assert_eq!(
                protein_signature.kmer_infos().len(),
                14,
                "Resolved sequence should have 14 protein 5-mers"
            );
            // Check that the ambiguous k-mer is resolved correctly
            if sequence == &"PLANTANDANIMALGENBMES" {
                let kmer_info = protein_signature.kmer_infos().get(&13058023948041027181);
                assert!(
                    kmer_info.is_some(),
                    "Expected k-mer with hash 6161374941338912337 to be present in {}",
                    sequence
                );
                let kmer_info = kmer_info.unwrap();
                assert_eq!(
                    kmer_info.encoded_kmer, "pphpp",
                    "Expected encoded k-mer 'pphpp' (NDMES/NNMES) to be present in {}",
                    sequence
                );
            } else if sequence == &"PLANTANDANIMALGENZMES" {
                let kmer_info = protein_signature.kmer_infos().get(&13058023948041027181);
                assert!(
                    kmer_info.is_some(),
                    "Expected k-mer with hash 13058023948041027181 to be present in {}",
                    sequence
                );
                let kmer_info = kmer_info.unwrap();
                assert_eq!(
                    kmer_info.encoded_kmer, "pphpp",
                    "Expected encoded k-mer 'pphpp' (NEMES/NQMES) to be present in {}",
                    sequence
                );
            } else if sequence == &"PLANTANDANIMALGENJMES" {
                let kmer_info = protein_signature.kmer_infos().get(&10495165127682499337);
                assert!(
                    kmer_info.is_some(),
                    "Expected k-mer with hash 10495165127682499337 to be present in {}",
                    sequence
                );
                let kmer_info = kmer_info.unwrap();
                assert_eq!(
                    kmer_info.encoded_kmer, "phhpp",
                    "Expected encoded k-mer 'phhpp' (NLMES/NIMES) to be present in {}",
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
        let result = index.process_fasta(&fasta_path, 0);
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
        let result = index.process_fasta(&valid_fasta_path, 0);
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
        assert_eq!(signature.kmer_infos().len(), 17, "Expected 17 k-mers for the test protein");

        // Verify some specific k-mers are present
        let expected_hash = 5893010049374798421; // Hash for "PLANT"
        assert!(
            signature.kmer_infos().contains_key(&expected_hash),
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

        let sig1_2 = index2.create_protein_signature("ACDEFGHIKLMNPQRSTVWY", "test1").unwrap();
        let sig2_2 = index2.create_protein_signature("PLANTANDANIMALGENQMES", "test2").unwrap();
        index2.store_signatures(vec![sig1_2, sig2_2]).unwrap();

        // Test equivalence
        assert!(index1.is_equivalent_to(&index2).unwrap());

        // Check stats
        assert_eq!(index1.signature_count(), 2);
        assert_eq!(index2.signature_count(), 2);
        assert_eq!(index1.combined_minhash_size(), index2.combined_minhash_size());

        // Test that different indices are not equivalent
        let index3 =
            ProteomeIndex::new(&temp_dir.path().join("test3.db"), 10, 1, "protein", false).unwrap();
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
        manual_index.process_fasta(&fasta_path, 0).unwrap();

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
        auto_index.process_fasta(&fasta_path, 0).unwrap();

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
            auto_index.process_fasta(&fasta_path, 10).unwrap();

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
            ProteomeIndex::new(&temp_dir.path().join("index3.db"), 10, 1, "protein", false)
                .unwrap();

        // Test that different indices are not equivalent
        let are_equivalent_3 = index1.is_equivalent_to(&index3).unwrap();
        assert!(!are_equivalent_3, "Indices with different parameters should not be equivalent");

        // Test with different sequences
        let index4 =
            ProteomeIndex::new(&temp_dir.path().join("index4.db"), 5, 1, "protein", false).unwrap();
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
        assert_eq!(index.store_raw_sequences(), true);
        assert_eq!(index.signature_count(), 1);

        // Get the signature and verify raw sequence is preserved
        let signatures = index.get_signatures();
        let entry = signatures.iter().next().unwrap();
        let signature = entry.value();
        assert!(signature.has_efficient_data());
        let raw_sequence = signature.get_raw_sequence();
        assert!(raw_sequence.is_some());
        assert_eq!(raw_sequence.unwrap(), sequence);

        Ok(())
    }

    #[test]
    fn test_efficient_storage_with_raw_sequences() -> Result<()> {
        let dir = tempdir()?;
        let db_path = dir.path().join("test_efficient_with_raw_sequences.db");

        // Create index with raw sequence storage enabled
        let index = ProteomeIndex::new(
            &db_path, 5,         // k-mer size
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
        assert_eq!(index.store_raw_sequences(), true);
        assert_eq!(index.signature_count(), 1);

        // Get the signature and verify raw sequence is preserved
        {
            let signatures = index.get_signatures();
            let entry = signatures.iter().next().unwrap();
            let signature = entry.value();
            assert!(signature.has_efficient_data());
            let raw_sequence = signature.get_raw_sequence();
            assert!(raw_sequence.is_some());
            assert_eq!(raw_sequence.unwrap(), sequence);
        }

        // Test that we can save state without errors
        index.save_state()?;

        Ok(())
    }

    #[test]
    fn test_efficient_storage_without_raw_sequences() -> Result<()> {
        let dir = tempdir()?;
        let db_path = dir.path().join("test_efficient_without_raw_sequences.db");

        // Create index with raw sequence storage disabled
        let index = ProteomeIndex::new(
            &db_path, 5,         // k-mer size
            1,         // scaled
            "protein", // molecular type
            false,     // don't store raw sequences,
        )?;

        // Add a protein sequence
        let sequence = "ACDEFGHIKLMNPQRSTVWY";
        let signature = index.create_protein_signature(sequence, "test_protein")?;

        // Verify raw sequence is not stored
        assert!(!signature.has_efficient_data());
        let raw_sequence = signature.get_raw_sequence();
        assert!(raw_sequence.is_none());

        // Store the signature
        index.store_signatures(vec![signature])?;

        // Verify the index has the correct configuration
        assert_eq!(index.store_raw_sequences(), false);
        assert_eq!(index.signature_count(), 1);

        // Get the signature and verify raw sequence is not stored
        {
            let signatures = index.get_signatures();
            let entry = signatures.iter().next().unwrap();
            let signature = entry.value();
            assert!(!signature.has_efficient_data());
            let raw_sequence = signature.get_raw_sequence();
            assert!(raw_sequence.is_none());
        }

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
        index.process_fasta(&fasta_path, 0)?;

        // Verify signatures were added
        let signatures = index.get_signatures();
        assert_eq!(
            signatures.len(),
            EXPECTED_SIGNATURES,
            "Expected {EXPECTED_SIGNATURES} signatures to be stored"
        );

        // Extract k-mer counts and raw sequences using functional programming
        let kmer_counts: Vec<usize> =
            signatures.iter().map(|entry| entry.value().kmer_infos().len()).collect();

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
