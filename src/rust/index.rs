use anyhow::{Context, Result};
use needletail::{parse_fastx_file, parse_fastx_stdin, parser::SequenceRecord};
use rayon::prelude::*;
use rocksdb::{Options, DB};
use serde::{Deserialize, Serialize};

use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::str;
use std::sync::{Arc, Mutex};

use crate::aminoacid::AminoAcidAmbiguity;
use crate::encoding::{get_encoding_fn_from_moltype, get_hash_function_from_moltype};
use crate::kmer::{KmerInfo, KmerPosition, KmerStats};
use crate::uniprot::UniProtEntry;
use sourmash::cmd::ComputeParameters;
use sourmash::collection::Collection;
use sourmash::manifest::Manifest;
use sourmash::signature::Signature;
use sourmash::storage::{FSStorage, InnerStorage, SigStore};
use sourmash_plugin_branchwater::search_significance::{
    compute_inverse_document_frequency, get_hash_frequencies, Normalization,
};
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

#[derive(Clone, Serialize, Deserialize)]
struct SerializableSignature {
    location: String,
    name: String,
    md5sum: String,
    minhash: KmerMinHash,
}

impl From<SmallSignature> for SerializableSignature {
    fn from(sig: SmallSignature) -> Self {
        Self {
            location: sig.location,
            name: sig.name,
            md5sum: sig.md5sum,
            minhash: sig.minhash,
        }
    }
}

impl From<SigStore> for SerializableSignature {
    fn from(sig: SigStore) -> Self {
        Self {
            location: sig.filename().clone(),
            name: sig.name().clone(),
            md5sum: sig.md5sum().to_string(),
            minhash: sig.minhash().unwrap().clone(),
        }
    }
}

// Represents the serializable state of ProteomeIndex
#[derive(Serialize, Deserialize)]
struct ProteomeIndexState {
    // Store signatures instead of Collection
    signatures: Vec<SerializableSignature>,
    combined_minhash: KmerMinHash,
    protein_info: HashMap<String, ProteinKmerInfo>,
    moltype: String,
    ksize: u32,
    scaled: u32,
}

// Represents a single protein's k-mer information
#[derive(Clone, Serialize, Deserialize)]
pub struct ProteinKmerInfo {
    // The original protein sequence and metadata
    pub protein: UniProtEntry,
    // Map of hash -> encoded k-mer for this protein
    pub hash_to_encoded: HashMap<u64, String>,
    // Map of encoded k-mer -> original k-mer info
    pub encoded_to_originals: HashMap<String, Vec<KmerInfo>>,
    // The protein's signature for searching
    pub signature: SerializableSignature,
}

impl ProteinKmerInfo {
    pub fn new(protein: UniProtEntry, signature: SmallSignature) -> Self {
        Self {
            protein,
            hash_to_encoded: HashMap::new(),
            encoded_to_originals: HashMap::new(),
            signature: SerializableSignature::from(signature),
        }
    }
}

pub struct ProteomeIndex {
    // RocksDB instance for persistent storage
    db: DB,

    // Collection of all signatures for searching
    collection: Arc<Mutex<Collection>>,

    // Combined minhash of all proteins for statistics
    combined_minhash: Arc<Mutex<KmerMinHash>>,

    // Map of signature md5 -> protein k-mer information
    protein_info: Arc<Mutex<HashMap<String, ProteinKmerInfo>>>,

    // Amino acid ambiguity handler
    aa_ambiguity: Arc<AminoAcidAmbiguity>,

    // Protein encoding function
    encoding_fn: fn(u8) -> u8,

    // Add moltype field for serialization
    moltype: String,

    // Add ksize field for serialization
    ksize: u32,

    // Add scaled field for serialization
    scaled: u32,
}

impl ProteomeIndex {
    pub fn new<P: AsRef<Path>>(path: P, ksize: u32, scaled: u32, moltype: &str) -> Result<Self> {
        // Create RocksDB options
        let mut opts = Options::default();
        opts.create_if_missing(true);

        // Open the database
        let db = DB::open(&opts, path)?;

        let hash_function = match get_hash_function_from_moltype(moltype) {
            Ok(value) => value,
            Err(value) => return value,
        };

        let encoding_fn = match get_encoding_fn_from_moltype(moltype) {
            Ok(value) => value,
            Err(value) => return value,
        };

        // Create the minhash sketch
        let minhash = KmerMinHash::new(
            scaled,
            ksize,
            hash_function,
            42,   // seed
            true, // track_abundance
            0,    // num (use scaled instead)
        );

        // Create an empty collection with storage
        let manifest = Manifest::default();
        let storage = InnerStorage::new(
            FSStorage::builder()
                .fullpath("".into())
                .subdir("".into())
                .build(),
        );
        let collection = Collection::new(manifest, storage);

        Ok(Self {
            db,
            collection: Arc::new(Mutex::new(collection)),
            combined_minhash: Arc::new(Mutex::new(minhash)),
            protein_info: Arc::new(Mutex::new(HashMap::new())),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            encoding_fn,
            moltype: moltype.to_string(),
            ksize,
            scaled,
        })
    }

    pub fn add_protein(&mut self, protein: UniProtEntry) -> Result<()> {
        // Get parameters from combined_minhash
        let minhash_guard = self.combined_minhash.lock().unwrap();
        let ksize = minhash_guard.ksize() as u32;
        let scaled = minhash_guard.scaled();
        drop(minhash_guard);

        // Create compute parameters for this protein's signature
        let params = ComputeParameters::builder()
            .ksizes(vec![ksize])
            .scaled(scaled)
            .protein(true)
            .num_hashes(0) // use scaled instead
            .build();

        // Create a signature for this protein
        let mut sig = Signature::from_params(&params);
        sig.set_name(&protein.id);

        // Convert to SmallSignature
        let small_sig = SmallSignature {
            location: sig.filename().clone(),
            name: sig.name().clone().unwrap(),
            md5sum: sig.md5sum().to_string(),
            minhash: sig.minhash().unwrap().clone(),
        };
        let md5sum = small_sig.md5sum.to_string();

        // Create protein-specific k-mer info
        let mut protein_info = ProteinKmerInfo::new(protein.clone(), small_sig);

        // Process k-mers for this protein
        self.process_protein_kmers(&mut protein_info)?;

        // Update global structures
        {
            let mut collection = self.collection.lock().unwrap();
            let mut sigs = vec![sig];
            let new_collection = Collection::from_sigs(sigs)?;
            *collection = new_collection;

            let mut combined = self.combined_minhash.lock().unwrap();
            combined.add_sequence(protein.sequence.as_bytes(), true)?;

            let mut proteins = self.protein_info.lock().unwrap();
            proteins.insert(md5sum.clone(), protein_info.clone());
        }

        // Store in RocksDB
        self.db.put(
            format!("protein:{}", md5sum).as_bytes(),
            bincode::serialize(&protein_info)?,
        )?;

        Ok(())
    }

    fn process_protein_kmers(&self, protein_info: &mut ProteinKmerInfo) -> Result<()> {
        let sequence = protein_info.protein.sequence.as_bytes();
        let minhash_guard = self.combined_minhash.lock().unwrap();
        let ksize = minhash_guard.ksize() as usize;
        drop(minhash_guard);

        for i in 0..sequence.len().saturating_sub(ksize - 1) {
            let kmer = &sequence[i..i + ksize];

            // Process the k-mer to get encoded version
            let (encoded_kmer, original_kmer) = self.process_kmer(kmer);

            // Create position information
            let kmer_position = KmerPosition {
                protein: protein_info.protein.clone(),
                position: i,
            };

            // Get the hash from the minhash implementation
            let mut minhash_guard = self.combined_minhash.lock().unwrap();
            minhash_guard.add_sequence(kmer, true)?;
            let hash = minhash_guard.mins().last().copied().unwrap_or(0);
            drop(minhash_guard);

            // Update protein-specific mappings
            protein_info
                .hash_to_encoded
                .insert(hash, encoded_kmer.clone());

            let kmer_infos = protein_info
                .encoded_to_originals
                .entry(encoded_kmer)
                .or_insert_with(Vec::new);

            // Find or create KmerInfo for this original k-mer
            if let Some(info) = kmer_infos
                .iter_mut()
                .find(|info| info.original_kmer == original_kmer)
            {
                info.positions.push(kmer_position);
            } else {
                kmer_infos.push(KmerInfo {
                    original_kmer,
                    positions: vec![kmer_position],
                });
            }
        }

        Ok(())
    }

    /// Process a k-mer to get its encoded version
    fn process_kmer(&self, kmer: &[u8]) -> (String, String) {
        // Convert the k-mer to a string
        let kmer_str = str::from_utf8(kmer).unwrap();

        // Create the encoded k-mer by applying the encoding function
        let encoded: Vec<u8> = kmer.iter().map(|&b| (self.encoding_fn)(b)).collect();
        let encoded_str = str::from_utf8(&encoded).unwrap();

        (encoded_str.to_string(), kmer_str.to_string())
    }

    /// Get all information about a hash value
    pub fn get_kmer_info(&self, hash: u64) -> Option<Vec<(String, Vec<KmerInfo>)>> {
        let proteins = self.protein_info.lock().unwrap();
        let mut results = Vec::new();

        for (md5sum, protein_info) in proteins.iter() {
            if let Some(encoded) = protein_info.hash_to_encoded.get(&hash) {
                if let Some(infos) = protein_info.encoded_to_originals.get(encoded) {
                    results.push((md5sum.clone(), infos.clone()));
                }
            }
        }

        if results.is_empty() {
            None
        } else {
            Some(results)
        }
    }

    /// Get statistics for a hash value
    pub fn get_kmer_stats(&self, hash: u64) -> Result<Option<KmerStats>> {
        let key = format!("stats:{}", hash);
        if let Some(value) = self.db.get(key.as_bytes())? {
            Ok(Some(bincode::deserialize(&value)?))
        } else {
            Ok(None)
        }
    }

    /// Compute and store k-mer statistics
    pub fn compute_statistics(&self) -> Result<()> {
        let minhash = self.combined_minhash.lock().unwrap();
        let collection = self.collection.lock().unwrap();

        // Get all signatures from the collection
        let mut signatures = Vec::new();
        for (_, record) in collection.iter() {
            if let Ok(sig) = collection.sig_from_record(record) {
                let serialized = SerializableSignature::from(sig);
                let small_sig = SmallSignature {
                    location: serialized.location,
                    name: serialized.name,
                    md5sum: serialized.md5sum,
                    minhash: serialized.minhash,
                };
                signatures.push(small_sig);
            }
        }

        // Compute hash frequencies using L1 normalization
        let frequencies = get_hash_frequencies(&minhash, Some(Normalization::L1));

        // Compute inverse document frequency with smoothing
        let idf = compute_inverse_document_frequency(&minhash, &signatures, Some(true));

        // Store statistics for each hash
        for hash in minhash.mins() {
            let stats = KmerStats {
                idf: *idf.get(&hash).unwrap_or(&0.0),
                frequency: *frequencies.get(&hash).unwrap_or(&0.0),
            };

            let key = format!("stats:{}", hash);
            let value = bincode::serialize(&stats)?;
            self.db.put(key.as_bytes(), value)?;
        }

        Ok(())
    }

    /// Process a single protein sequence
    pub fn process_sequence(&self, sequence: &[u8], protein: UniProtEntry) -> Result<()> {
        // Validate sequence - check for invalid amino acids
        for (i, c) in sequence.iter().enumerate() {
            if !self.aa_ambiguity.is_valid_aa(*c as char) {
                return Err(anyhow::anyhow!(
                    "Invalid amino acid '{}' at position {}",
                    *c as char,
                    i + 1
                ));
            }
        }

        // Create a temporary mutable copy of self to add the protein
        let mut opts = Options::default();
        opts.create_if_missing(true);
        let mut temp_self = Self {
            db: DB::open(&opts, self.db.path())?,
            collection: self.collection.clone(),
            combined_minhash: self.combined_minhash.clone(),
            protein_info: self.protein_info.clone(),
            aa_ambiguity: self.aa_ambiguity.clone(),
            encoding_fn: self.encoding_fn,
            moltype: self.moltype.clone(),
            ksize: self.ksize,
            scaled: self.scaled,
        };

        // Add the protein
        temp_self.add_protein(protein)?;

        Ok(())
    }

    /// Get all hashes from the combined minhash
    pub fn get_hashes(&self) -> Vec<u64> {
        let minhash = self.combined_minhash.lock().unwrap();
        minhash.mins().to_vec()
    }

    /// Save the index state to RocksDB
    pub fn save(&self) -> Result<()> {
        // Extract signatures from collection
        let collection = self.collection.lock().unwrap();
        let mut signatures = Vec::new();

        // Convert each signature in the collection to a SmallSignature
        for (_, record) in collection.iter() {
            if let Ok(sig) = collection.sig_from_record(record) {
                signatures.push(SerializableSignature::from(sig));
            }
        }

        // Create a serializable state
        let state = ProteomeIndexState {
            signatures: signatures
                .into_iter()
                .map(SerializableSignature::from)
                .collect(),
            combined_minhash: self.combined_minhash.lock().unwrap().clone(),
            protein_info: self.protein_info.lock().unwrap().clone(),
            moltype: self.moltype.clone(),
            ksize: self.ksize,
            scaled: self.scaled,
        };

        // Serialize and store the state using bincode 1.3
        let serialized = bincode::serialize(&state)?;
        self.db.put("index_state", serialized)?;

        Ok(())
    }

    /// Load the index state from RocksDB
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut opts = Options::default();
        opts.create_if_missing(false);
        let db = DB::open(&opts, path)?;

        // Load the serialized state using bincode 1.3
        let state: ProteomeIndexState = match db.get("index_state")? {
            Some(bytes) => bincode::deserialize(&bytes)?,
            None => return Err(anyhow::anyhow!("No index state found in database")),
        };

        // Get the encoding function
        let encoding_fn = match get_encoding_fn_from_moltype(&state.moltype) {
            Ok(value) => value,
            Err(value) => return value,
        };

        // Reconstruct collection from signatures
        let mut collection = Collection::new(
            Manifest::default(),
            InnerStorage::new(
                FSStorage::builder()
                    .fullpath("".into())
                    .subdir("".into())
                    .build(),
            ),
        );

        // Convert SerializableSignature to Signature
        let mut sigs = Vec::new();
        for serialized_sig in state.signatures {
            // Create a new signature with the same parameters
            let mut sig = Signature::from_params(
                &ComputeParameters::builder()
                    .ksizes(vec![state.ksize])
                    .scaled(state.scaled)
                    .protein(true)
                    .num_hashes(0)
                    .build(),
            );

            // Set the signature properties
            sig.set_name(&serialized_sig.name);
            sig.set_filename(&serialized_sig.location);

            sigs.push(sig);
        }
        collection = Collection::from_sigs(sigs)?;

        Ok(Self {
            db,
            collection: Arc::new(Mutex::new(collection)),
            combined_minhash: Arc::new(Mutex::new(state.combined_minhash)),
            protein_info: Arc::new(Mutex::new(state.protein_info)),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            encoding_fn,
            moltype: state.moltype,
            ksize: state.ksize,
            scaled: state.scaled,
        })
    }

    /// Close the database properly
    pub fn close(self) -> Result<()> {
        // Save the state before closing
        self.save()?;
        drop(self.db);
        Ok(())
    }

    /// Process a list of protein FASTA files
    pub fn process_protein_files<P: AsRef<Path>>(&mut self, paths: &[P]) -> Result<()> {
        for path in paths {
            // Create a FASTX reader from the file or stdin
            let mut fastx_reader = if path.as_ref() == Path::new("-") {
                parse_fastx_stdin().context("Failed to parse FASTA/FASTQ data from stdin")?
            } else {
                parse_fastx_file(path).context("Failed to open file for FASTA/FASTQ data")?
            };

            let mut record_count: u64 = 0;

            // Parse records and add sequences to signatures
            while let Some(record_result) = fastx_reader.next() {
                let record = record_result.context("Failed to read a record from input")?;
                let id = String::from_utf8_lossy(record.id()).to_string();
                let seq_bytes = record.seq();
                let sequence = String::from_utf8_lossy(seq_bytes.as_ref()).to_string();

                let protein = UniProtEntry {
                    id,
                    sequence: sequence.clone(),
                    features: Vec::new(),
                    ..Default::default()
                };

                self.process_sequence(seq_bytes.as_ref(), protein)?;
                record_count += 1;
            }
            println!("Processed {} sequence records", record_count);
        }
        Ok(())
    }

    /// Process a UniProt XML file
    pub fn process_uniprot_xml<P: AsRef<Path>>(&mut self, path: P) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut parser = quick_xml::Reader::from_reader(reader);
        let mut buf = Vec::new();

        let mut current_protein = None;
        let mut current_sequence = String::new();

        loop {
            match parser.read_event_into(&mut buf) {
                Ok(quick_xml::events::Event::Start(ref e)) => match e.name().as_ref() {
                    b"entry" => {
                        current_protein = Some(UniProtEntry::default());
                    }
                    b"sequence" => {
                        current_sequence.clear();
                    }
                    _ => {}
                },
                Ok(quick_xml::events::Event::Text(e)) => {
                    if let Some(_protein) = current_protein.as_mut() {
                        if !current_sequence.is_empty() {
                            current_sequence.push_str(&e.unescape()?.into_owned());
                        }
                    }
                }
                Ok(quick_xml::events::Event::End(ref e)) => match e.name().as_ref() {
                    b"entry" => {
                        if let Some(mut protein) = current_protein.take() {
                            protein.sequence = current_sequence.clone();
                            self.process_sequence(current_sequence.as_bytes(), protein)?;
                        }
                        current_sequence.clear();
                    }
                    _ => {}
                },
                Ok(quick_xml::events::Event::Eof) => break,
                Err(e) => return Err(anyhow::anyhow!("Error parsing XML: {}", e)),
                _ => {}
            }
        }
        Ok(())
    }
}
