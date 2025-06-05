use anyhow::{bail, Result};
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use rand::Rng;
use rayon::prelude::*;
use rocksdb::{Options, DB};
use serde::{Deserialize, Serialize};
use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::sync::{Arc, Mutex};

use crate::aminoacid::AminoAcidAmbiguity;
use crate::uniprot::{self, UniProtSequence};

/// Represents the position of a k-mer in a protein sequence
#[derive(Debug, Clone)]
pub struct KmerPosition {
    protein: UniProtSequence,
    position: usize, // 0-based position in sequence
}

impl KmerPosition {
    /// Get all features that overlap with this k-mer
    pub fn overlapping_features(&self, kmer_length: usize) -> Vec<&uniprot::ProteinFeature> {
        let kmer_start = self.position + 1; // Convert to 1-based position
        let kmer_end = kmer_start + kmer_length - 1;

        self.protein
            .features
            .iter()
            .filter(|feature| {
                // Check if the feature overlaps with the k-mer
                !(feature.end < kmer_start || feature.start > kmer_end)
            })
            .collect()
    }
}

/// Represents information about a k-mer occurrence
#[derive(Debug, Clone)]
pub struct KmerInfo {
    original_kmer: String,
    positions: Vec<KmerPosition>,
}

/// Statistics for k-mer frequency analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerStats {
    idf: f64,               // Inverse document frequency
    frequency: f64,         // Raw frequency
    l1_norm_frequency: f64, // L1 normalized frequency
    l2_norm_frequency: f64, // L2 normalized frequency
}

pub struct ProteomeIndex {
    // RocksDB instance for persistent storage
    db: DB,
    // First level: hash -> encoded k-mer
    hash_to_encoded: Arc<Mutex<HashMap<u64, String>>>,
    // Second level: encoded k-mer -> original k-mers and their positions
    encoded_to_originals: Arc<Mutex<HashMap<String, Vec<KmerInfo>>>>,
    // Store the underlying MinHash sketch
    minhash: Arc<Mutex<KmerMinHash>>,
    // Amino acid ambiguity handler
    aa_ambiguity: Arc<AminoAcidAmbiguity>,
    // Protein encoding function
    encoding_fn: fn(u8) -> u8,
}

#[derive(Debug)]
pub enum IndexError {
    IoError(io::Error),
}

impl From<io::Error> for IndexError {
    fn from(error: io::Error) -> Self {
        IndexError::IoError(error)
    }
}

impl ProteomeIndex {
    pub fn new<P: AsRef<Path>>(
        ksize: u32,
        scaled: u32,
        _moltype: &str,
        encoding_type: &str,
        db_path: P,
    ) -> Result<Self> {
        // Create RocksDB options
        let mut opts = Options::default();
        opts.create_if_missing(true);
        opts.set_compression_type(rocksdb::DBCompressionType::Lz4);

        // Open RocksDB
        let db = DB::open(&opts, db_path)?;

        // Create MinHash sketch
        let hash_function = match encoding_type {
            "dayhoff" => HashFunctions::Murmur64Dayhoff,
            "hp" => HashFunctions::Murmur64Hp,
            _ => HashFunctions::Murmur64Protein,
        };

        let minhash = KmerMinHash::new(
            scaled as u64,
            ksize,
            hash_function,
            42,    // seed
            false, // track_abundance
            0,     // num
        );

        // Select encoding function
        let encoding_fn = match encoding_type {
            "dayhoff" => aa_to_dayhoff,
            "hp" => aa_to_hp,
            _ => |c| c,
        };

        Ok(ProteomeIndex {
            db,
            hash_to_encoded: Arc::new(Mutex::new(HashMap::new())),
            encoded_to_originals: Arc::new(Mutex::new(HashMap::new())),
            minhash: Arc::new(Mutex::new(minhash)),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            encoding_fn,
        })
    }

    /// Process a UniProt XML file
    pub fn process_uniprot_xml<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let proteins = UniProtSequence::from_xml(path)?;

        // Process proteins in parallel
        proteins
            .par_iter()
            .try_for_each(|protein| self.process_sequence(&protein.sequence, protein.clone()))?;

        // Compute and store statistics
        self.compute_and_store_statistics()?;

        Ok(())
    }

    /// Process a sequence with its protein information
    pub fn process_sequence(&self, sequence: &str, protein: UniProtSequence) -> Result<()> {
        // First validate the sequence
        self.aa_ambiguity.validate_sequence(sequence)?;

        let k = {
            let minhash = self.minhash.lock().unwrap();
            minhash.ksize() as usize
        };

        // Process k-mers in parallel
        (0..sequence.len().saturating_sub(k - 1))
            .into_par_iter()
            .filter(|&i| i + k <= sequence.len())
            .for_each(|i| {
                let kmer = &sequence[i..i + k];
                self.add_kmer(kmer, i, protein.clone());
            });

        Ok(())
    }

    /// Process a list of protein FASTA files
    pub fn process_protein_files<P: AsRef<Path> + Sync>(&self, files: &[P]) -> Result<()> {
        // Process files in parallel
        files
            .par_iter()
            .try_for_each(|file_path| self.process_fasta_file(file_path.as_ref()))?;

        // Compute and store statistics
        self.compute_and_store_statistics()?;

        Ok(())
    }

    fn process_fasta_file(&self, path: &Path) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        // Buffer for collecting sequences
        let mut sequences = Vec::new();
        let mut current_sequence = String::new();
        let mut current_id = String::new();

        // First pass: collect sequences
        for line in reader.lines() {
            let line = line?;

            if line.starts_with('>') {
                if !current_sequence.is_empty() {
                    sequences.push((current_id.clone(), current_sequence.clone()));
                    current_sequence.clear();
                }
                current_id = line[1..].trim().to_string();
            } else {
                current_sequence.push_str(line.trim());
            }
        }

        // Don't forget the last sequence
        if !current_sequence.is_empty() {
            sequences.push((current_id.clone(), current_sequence.clone()));
        }

        // Second pass: process sequences in parallel
        sequences.into_par_iter().try_for_each(|(id, seq)| {
            self.process_sequence(
                &seq,
                UniProtSequence {
                    id,
                    accession: String::new(), // We don't have accession in FASTA
                    sequence: seq.clone(),
                    features: Vec::new(),
                },
            )
            .map_err(|e| anyhow::anyhow!("Error processing sequence: {}", e))
        })?;

        Ok(())
    }

    /// Process a k-mer to get its encoded version
    fn process_kmer(&self, kmer: &str) -> (String, String) {
        // First resolve ambiguities
        let resolved: String = kmer
            .chars()
            .map(|aa| self.aa_ambiguity.resolve_ambiguity(aa))
            .collect();

        // Then encode using the selected encoding function
        let encoded: String = resolved
            .bytes()
            .map(|aa| (self.encoding_fn)(aa))
            .map(|aa| aa as char)
            .collect();

        (encoded, kmer.to_string())
    }

    /// Add a k-mer with its position information
    fn add_kmer(&self, kmer: &str, position: usize, protein: UniProtSequence) {
        // Process the k-mer to get encoded version
        let (encoded_kmer, original_kmer) = self.process_kmer(kmer);

        // Get the hash from the minhash implementation
        let hash = {
            let mut minhash = self.minhash.lock().unwrap();
            minhash.add_sequence(encoded_kmer.as_bytes(), true).unwrap();
            minhash.to_vec().last().copied().unwrap_or(0)
        };

        // Create position information
        let kmer_position = KmerPosition { protein, position };

        // Update hash -> encoded mapping
        {
            let mut hash_to_encoded = self.hash_to_encoded.lock().unwrap();
            hash_to_encoded.insert(hash, encoded_kmer.clone());
        }

        // Update encoded -> originals mapping
        {
            let mut encoded_to_originals = self.encoded_to_originals.lock().unwrap();
            let kmer_infos = encoded_to_originals
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
    }

    /// Get all information about a hash value
    pub fn get_kmer_info(&self, hash: u64) -> Option<(String, Vec<KmerInfo>)> {
        let hash_to_encoded = self.hash_to_encoded.lock().unwrap();
        let encoded_to_originals = self.encoded_to_originals.lock().unwrap();

        hash_to_encoded.get(&hash).and_then(|encoded| {
            encoded_to_originals
                .get(encoded)
                .map(|infos| (encoded.clone(), infos.clone()))
        })
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

    /// Get all hashes in the index
    pub fn get_mins(&self) -> Vec<u64> {
        let minhash = self.minhash.lock().unwrap();
        minhash.to_vec()
    }

    /// Compute and store k-mer statistics
    fn compute_and_store_statistics(&self) -> Result<()> {
        let minhash = self.minhash.lock().unwrap();
        let hashes = minhash.to_vec();

        // TODO: Implement frequency computations
        // For now, just store empty statistics
        for hash in hashes {
            let stats = KmerStats {
                idf: 0.0,
                frequency: 0.0,
                l1_norm_frequency: 0.0,
                l2_norm_frequency: 0.0,
            };

            let key = format!("stats:{}", hash);
            let value = bincode::serialize(&stats)?;
            self.db.put(key.as_bytes(), value)?;
        }

        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub enum Normalization {
    L1,
    L2,
}
