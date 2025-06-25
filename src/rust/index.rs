use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, Mutex};

use anyhow::Result;
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

    /// Add a single protein sequence as a signature and process its k-mers
    ///
    /// This method creates a protein signature from the given sequence, processes its k-mers
    /// to extract detailed position information, and stores it in the index. The signature
    /// is also merged into the combined minhash for statistical analysis.
    ///
    /// # Arguments
    ///
    /// * `sequence` - The protein sequence as a string
    /// * `_name` - The name/identifier for the protein (currently unused but reserved for future use)
    ///
    /// # Returns
    ///
    /// Returns `Ok(())` on success, or an error if the operation fails.
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
    ///     index.add_protein_sequence(sequence, "test_protein")?;
    ///     
    ///     // The protein is now indexed and its k-mers are processed
    ///     Ok(())
    /// }
    /// ```
    pub fn add_protein_sequence(&self, sequence: &str, _name: &str) -> Result<()> {
        // Create a new protein signature
        let mut protein_sig =
            ProteinSignature::new(self.ksize, self.scaled, &self.moltype, self.seed)?;

        // Add the protein sequence to the signature
        protein_sig.add_protein(sequence.as_bytes())?;

        // Process the k-mers to get detailed k-mer information
        let processed_signature = self.process_protein_kmers(sequence, &protein_sig.signature())?;

        // Get the MD5 hash of the signature for storage
        let md5sum = protein_sig.signature().get_md5sum();

        // Store the processed signature in the signatures map
        {
            let mut signatures = self.signatures.lock().unwrap();
            signatures.insert(md5sum.to_string(), processed_signature);
        }

        Ok(())
    }

    pub fn process_protein_kmers<S: SignatureAccess>(
        &self,
        sequence: &str,
        signature: &S,
    ) -> Result<ProteinSignature> {
        let ksize = self.ksize as usize;
        let seed = self.seed as u64;
        let hashvals = &signature.get_minhash().to_vec();

        let mut protein_signature =
            ProteinSignature::new(self.ksize, self.scaled, &self.moltype, self.seed)?;

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
                    protein_signature.kmer_infos.entry(hashval).or_insert_with(|| KmerInfo {
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

        Ok(protein_signature)
    }
}
