use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Duration;

use anyhow::{bail, Result};
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
use crate::kmer::KmerInfo;
use crate::signature::{ProteinSignature, SignatureAccess};

/// Statistics for k-mer frequency analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteomeIndexKmerStats {
    pub idf: HashMap<u64, f64>, // Inverse document frequency for each k-mer hashvalue
    pub frequency: HashMap<u64, f64>, // Raw frequency for each k-mer hashvalue
}

// Represents the serializable state of ProteomeIndex
#[derive(Serialize, Deserialize)]
struct ProteomeIndexState {
    // Store signatures instead of Collection
    signatures: Vec<ProteinSignature>,
    combined_minhash: KmerMinHash,
    moltype: String,
    ksize: u32,
    scaled: u32,
    seed: u64,
}

pub struct ProteomeIndex {
    // RocksDB instance for persistent storage
    db: DB,

    // Combined minhash of all proteins for statistics
    combined_minhash: Arc<Mutex<KmerMinHash>>,

    // Map of signature md5 -> protein signature
    signatures: Arc<Mutex<HashMap<String, ProteinSignature>>>,

    // Amino acid ambiguity handler
    aa_ambiguity: Arc<AminoAcidAmbiguity>,

    // Protein encoding function
    encoding_fn: fn(u8) -> u8,

    // Statistics for k-mer frequencies and IDF
    stats: ProteomeIndexKmerStats,

    // Add moltype field for serialization
    moltype: String,

    // Add ksize field for serialization
    ksize: u32,

    // Add minhash_ksize field for serialization
    // MinHash k-mer size is protein_ksize * 3, as a legacy from Sourmash which was originally designed for DNA
    minhash_ksize: u32,

    // Add scaled field for serialization
    scaled: u32,

    // Add seed field for serialization
    seed: u64,
}

impl ProteomeIndex {
    pub fn new<P: AsRef<Path>>(
        path: P,
        ksize: u32,
        scaled: u32,
        moltype: &str,
        seed: u64,
    ) -> Result<Self> {
        // Create RocksDB options
        let mut opts = Options::default();
        opts.create_if_missing(true);

        // Open the database
        let db = DB::open(&opts, path)?;

        let hash_function = get_hash_function_from_moltype(moltype)?;

        let encoding_fn = get_encoding_fn_from_moltype(moltype)?;

        let minhash_ksize = ksize * 3;
        // Create the minhash sketch
        let minhash = KmerMinHash::new(
            scaled,
            minhash_ksize,
            hash_function,
            seed, // seed
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
            signatures: Arc::new(Mutex::new(HashMap::new())),
            combined_minhash: Arc::new(Mutex::new(minhash)),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            encoding_fn,
            moltype: moltype.to_string(),
            ksize: ksize,
            minhash_ksize: minhash_ksize,
            scaled,
            stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
            seed,
        })
    }

    /// Get a reference to the signatures map (for testing)
    pub fn get_signatures(&self) -> &Arc<Mutex<HashMap<String, ProteinSignature>>> {
        &self.signatures
    }

    /// Get a reference to the combined minhash (for testing)
    pub fn get_combined_minhash(&self) -> &Arc<Mutex<KmerMinHash>> {
        &self.combined_minhash
    }

    /// Save the current index state to RocksDB
    pub fn save_state(&self) -> Result<()> {
        let signatures_map = self.signatures.lock().unwrap();
        let combined_minhash = self.combined_minhash.lock().unwrap();

        let state = ProteomeIndexState {
            signatures: signatures_map.values().cloned().collect(),
            combined_minhash: combined_minhash.clone(),
            moltype: self.moltype.clone(),
            ksize: self.ksize,
            scaled: self.scaled,
            seed: self.seed,
        };

        let serialized = bincode::serialize(&state)?;
        self.db.put(b"index_state", &serialized)?;

        Ok(())
    }

    /// Load index state from RocksDB
    pub fn load_state(&self) -> Result<()> {
        let serialized = self.db.get(b"index_state")?;
        if let Some(data) = serialized {
            let state: ProteomeIndexState = bincode::deserialize(&data)?;

            // Update signatures
            {
                let mut signatures_map = self.signatures.lock().unwrap();
                signatures_map.clear();
                for sig in state.signatures {
                    let md5sum = sig.signature().md5sum.clone();
                    signatures_map.insert(md5sum.to_string(), sig);
                }
            }

            // Update combined minhash
            {
                let mut combined_minhash = self.combined_minhash.lock().unwrap();
                *combined_minhash = state.combined_minhash;
            }

            Ok(())
        } else {
            bail!("No saved state found in database")
        }
    }

    /// Load an existing ProteomeIndex from a RocksDB path
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        // Create RocksDB options
        let mut opts = Options::default();
        opts.create_if_missing(false); // Don't create if missing

        // Open the database
        let db = DB::open(&opts, path)?;

        // Try to load state to get configuration
        let serialized = db.get(b"index_state")?;
        if let Some(data) = serialized {
            let state: ProteomeIndexState = bincode::deserialize(&data)?;

            let hash_function = get_hash_function_from_moltype(&state.moltype)?;
            let encoding_fn = get_encoding_fn_from_moltype(&state.moltype)?;

            let index = Self {
                db,
                signatures: Arc::new(Mutex::new(HashMap::new())),
                combined_minhash: Arc::new(Mutex::new(state.combined_minhash.clone())),
                aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
                encoding_fn,
                moltype: state.moltype,
                ksize: state.ksize,
                minhash_ksize: state.ksize * 3,
                scaled: state.scaled,
                stats: ProteomeIndexKmerStats { idf: HashMap::new(), frequency: HashMap::new() },
                seed: state.seed,
            };

            // Load the actual state
            index.load_state()?;

            Ok(index)
        } else {
            bail!("No saved state found in database")
        }
    }

    /// Get the number of signatures in the index
    pub fn signature_count(&self) -> usize {
        self.signatures.lock().unwrap().len()
    }

    /// Get the combined minhash size
    pub fn combined_minhash_size(&self) -> usize {
        self.combined_minhash.lock().unwrap().size()
    }

    /// Compare this index with another for equivalency
    pub fn is_equivalent_to(&self, other: &ProteomeIndex) -> Result<bool> {
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
        if self.seed != other.seed {
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

        // Compare signatures
        let self_signatures = self.signatures.lock().unwrap();
        let other_signatures = other.signatures.lock().unwrap();

        for (md5, self_sig) in self_signatures.iter() {
            if let Some(other_sig) = other_signatures.get(md5) {
                let self_mins = self_sig.signature().get_minhash().mins();
                let other_mins = other_sig.signature().get_minhash().mins();
                if self_mins != other_mins {
                    return Ok(false);
                }
            } else {
                return Ok(false);
            }
        }

        // Compare combined minhashes
        let self_combined = self.combined_minhash.lock().unwrap();
        let other_combined = other.combined_minhash.lock().unwrap();

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
        println!("  Seed: {}", self.seed);
        println!("  Number of signatures: {}", self.signature_count());
        println!("  Combined minhash size: {}", self.combined_minhash_size());
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
    ///         42,       // seed
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
    pub fn create_protein_signature(&self, sequence: &str, name: &str) -> Result<ProteinSignature> {
        // Validate and resolve ambiguity if needed
        let processed_sequence = self.aa_ambiguity.validate_and_resolve(sequence)?;

        // Create a new protein signature
        let mut protein_sig =
            ProteinSignature::new(name, self.ksize, self.scaled, &self.moltype, self.seed)?;

        // Add the protein sequence to the signature
        protein_sig.add_protein(processed_sequence.as_bytes())?;

        // Process the k-mers to get detailed k-mer information
        self.process_kmers(&processed_sequence, &mut protein_sig)?;

        // Return the processed signature (don't store it yet)
        Ok(protein_sig)
    }

    pub fn process_kmers(
        &self,
        sequence: &str,
        protein_signature: &mut ProteinSignature,
    ) -> Result<()> {
        let ksize = self.ksize as usize;
        let seed = self.seed as u64;
        let hashvals = &protein_signature.signature().get_minhash().to_vec();

        for i in 0..sequence.len().saturating_sub(ksize - 1) {
            let kmer = &sequence[i..i + ksize];

            // Process the k-mer to get encoded version
            let (encoded_kmer, original_kmer) =
                encode_kmer_with_encoding_fn(kmer, self.encoding_fn)?;

            // Get the hash from the minhash implementation
            let hashval = _hash_murmur(encoded_kmer.as_bytes(), seed);

            // If this hashval is in the minhash, then save its k-mer positions
            if hashvals.contains(&hashval) {
                let kmer_info =
                    protein_signature.kmer_infos_mut().entry(hashval).or_insert_with(|| KmerInfo {
                        ksize: ksize,
                        hashval: hashval,
                        encoded_kmer: encoded_kmer.clone(),
                        original_kmer_to_position: HashMap::new(),
                    });

                kmer_info
                    .original_kmer_to_position
                    .entry(original_kmer)
                    .or_insert_with(Vec::new)
                    .push(i);
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
    pub fn store_signatures(&self, protein_signatures: Vec<ProteinSignature>) -> Result<()> {
        // Collect the minhash data from new signatures before storing them
        let new_hashes_and_abunds: Vec<(u64, u64)> = protein_signatures
            .iter()
            .flat_map(|sig| {
                let minhash = sig.signature().get_minhash();
                // If abundance is tracked, use to_vec_abunds, else use mins with abundance 1
                if let Some(abunds) = minhash.abunds() {
                    minhash.mins().into_iter().zip(abunds.into_iter()).collect::<Vec<_>>()
                } else {
                    minhash.mins().into_iter().map(|h| (h, 1)).collect::<Vec<_>>()
                }
            })
            .collect();

        // Store all signatures in the signatures map
        {
            let mut signatures_map = self.signatures.lock().unwrap();
            for protein_signature in protein_signatures {
                let md5sum = protein_signature.signature().md5sum.clone();
                signatures_map.insert(md5sum.to_string(), protein_signature);
            }
        }

        // Update the combined minhash with only the new signatures
        let mut combined_minhash = self.combined_minhash.lock().unwrap();
        combined_minhash.add_many_with_abund(&new_hashes_and_abunds)?;

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
    ) -> Result<()> {
        use needletail::parse_fastx_file;
        use rayon::prelude::*;

        if progress_interval > 0 {
            println!("Reading FASTA file...");
        }

        // Open and parse the FASTA file using needletail
        let mut reader = parse_fastx_file(&fasta_path)?;

        // Convert records into a vector for parallel processing
        let mut records = Vec::new();
        let mut record_count = 0;
        while let Some(record) = reader.next() {
            let record = record?;
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
        let signatures: Result<Vec<ProteinSignature>> = records
            .par_iter()
            .map(|(seq_bytes, id_bytes)| {
                let sequence = std::str::from_utf8(seq_bytes)?;
                let name = std::str::from_utf8(id_bytes)?;

                // Create protein signature for each sequence
                let result = self.create_protein_signature(sequence, name);
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
