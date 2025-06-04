use anyhow::{bail, Result};
use rand::seq::SliceRandom;
use rand::thread_rng;
use rayon::prelude::*;
use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// Standard amino acids and their properties
const STANDARD_AA: [char; 20] = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'Y',
];

/// Represents amino acid ambiguity codes and their possible resolutions
#[derive(Debug)]
struct AminoAcidAmbiguity {
    /// Maps ambiguous amino acid codes to their possible resolutions
    replacements: HashMap<char, Vec<char>>,
}

impl AminoAcidAmbiguity {
    fn new() -> Self {
        let mut replacements = HashMap::new();

        // Standard amino acids map to themselves
        for aa in STANDARD_AA.iter() {
            replacements.insert(*aa, vec![*aa]);
        }

        // Add ambiguous amino acid codes
        replacements.insert('X', STANDARD_AA.to_vec()); // Unknown - any amino acid
        replacements.insert('B', vec!['D', 'N']); // Aspartic acid or Asparagine
        replacements.insert('Z', vec!['E', 'Q']); // Glutamic acid or Glutamine
        replacements.insert('J', vec!['I', 'L']); // Isoleucine or Leucine

        AminoAcidAmbiguity { replacements }
    }

    fn is_valid_aa(&self, aa: char) -> bool {
        self.replacements.contains_key(&aa)
    }

    fn resolve_ambiguity(&self, aa: char) -> char {
        if let Some(possible_aas) = self.replacements.get(&aa) {
            // If there's only one possibility, use it
            if possible_aas.len() == 1 {
                possible_aas[0]
            } else {
                // Otherwise randomly choose one of the possibilities
                *possible_aas.choose(&mut thread_rng()).unwrap_or(&aa)
            }
        } else {
            // If not found in replacements, return the original character
            // (this shouldn't happen if is_valid_aa is called first)
            aa
        }
    }
}

pub struct ProteomeIndex {
    // Store hash -> encoded k-mer -> original k-mer mappings
    hash_to_kmers: Arc<Mutex<HashMap<u64, (String, String)>>>,
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
    ParseError(String),
}

impl From<io::Error> for IndexError {
    fn from(error: io::Error) -> Self {
        IndexError::IoError(error)
    }
}

impl ProteomeIndex {
    pub fn new(ksize: u32, scaled: u64, moltype: &str, encoding_type: &str) -> Self {
        // Create a new KmerMinHash with protein parameters
        let hash_function = match encoding_type {
            "dayhoff" => HashFunctions::Murmur64Dayhoff,
            "hp" => HashFunctions::Murmur64Hp,
            _ => HashFunctions::Murmur64Protein,
        };

        let minhash = KmerMinHash::new(
            scaled as u64,
            ksize,
            hash_function.clone(),
            42,    // seed
            false, // track_abundance
            0,     // num
        );

        // Select the encoding function based on type
        let encoding_fn = match encoding_type {
            "dayhoff" => aa_to_dayhoff,
            "hp" => aa_to_hp,
            _ => |c| c, // Raw encoding - no transformation
        };

        ProteomeIndex {
            hash_to_kmers: Arc::new(Mutex::new(HashMap::new())),
            minhash: Arc::new(Mutex::new(minhash)),
            aa_ambiguity: Arc::new(AminoAcidAmbiguity::new()),
            encoding_fn,
        }
    }

    /// Validates a protein sequence and returns an error if invalid characters are found
    fn validate_sequence(&self, sequence: &str) -> Result<()> {
        for (i, c) in sequence.chars().enumerate() {
            if !self.aa_ambiguity.is_valid_aa(c) {
                bail!("Invalid amino acid '{}' found at position {}", c, i + 1);
            }
        }
        Ok(())
    }

    /// Resolves ambiguous amino acids in a k-mer and returns both the resolved and encoded versions
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

    pub fn add_kmer(&self, kmer: &str) {
        // Process the k-mer to get encoded and original versions
        let (encoded_kmer, original_kmer) = self.process_kmer(kmer);

        // Get the hash from the minhash implementation
        let hash = {
            let mut minhash = self.minhash.lock().unwrap();
            minhash.add_sequence(encoded_kmer.as_bytes(), true).unwrap();
            minhash.to_vec().last().copied().unwrap_or(0)
        };

        // Store both the encoded and original k-mers
        let mut hash_to_kmers = self.hash_to_kmers.lock().unwrap();
        hash_to_kmers.insert(hash, (encoded_kmer, original_kmer));
    }

    pub fn get_kmers(&self, hash: u64) -> Option<(String, String)> {
        let hash_to_kmers = self.hash_to_kmers.lock().unwrap();
        hash_to_kmers.get(&hash).cloned()
    }

    pub fn process_sequence_chunk(&self, chunk: &str) -> Result<()> {
        // First validate the entire chunk
        self.validate_sequence(chunk)?;

        let k = {
            let minhash = self.minhash.lock().unwrap();
            minhash.ksize() as usize
        };

        // Process k-mers sequentially for now to avoid locking issues
        for i in 0..chunk.len().saturating_sub(k - 1) {
            if i + k <= chunk.len() {
                let kmer = &chunk[i..i + k];
                self.add_kmer(kmer);
            }
        }

        Ok(())
    }

    fn process_fasta_file(&self, path: &Path) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        // Buffer for sequence data
        let mut current_sequence = String::new();
        let mut line_number = 0;

        // Process each line
        for line in reader.lines() {
            line_number += 1;
            let line = line?;

            if line.starts_with('>') {
                // Process previous sequence if it exists
                if !current_sequence.is_empty() {
                    self.process_sequence_chunk(&current_sequence)
                        .map_err(|e| {
                            anyhow::anyhow!(
                                "Error in sequence ending at line {}: {}",
                                line_number - 1,
                                e
                            )
                        })?;
                    current_sequence.clear();
                }
            } else {
                // Append sequence data
                current_sequence.push_str(line.trim());
            }
        }

        // Process the last sequence
        if !current_sequence.is_empty() {
            self.process_sequence_chunk(&current_sequence)
                .map_err(|e| {
                    anyhow::anyhow!("Error in sequence ending at line {}: {}", line_number, e)
                })?;
        }

        Ok(())
    }

    /// Process a list of protein FASTA files
    pub fn process_protein_files<P: AsRef<Path>>(&self, files: &[P]) -> Result<()> {
        // Process files sequentially for now
        for file_path in files {
            self.process_fasta_file(file_path.as_ref())?;
        }
        Ok(())
    }

    pub fn add_sequence(&mut self, seq: &str) {
        // Process sequence in chunks to avoid memory issues with very long sequences
        const CHUNK_SIZE: usize = 1_000_000;

        seq.as_bytes()
            .chunks(CHUNK_SIZE)
            .filter_map(|chunk| std::str::from_utf8(chunk).ok())
            .for_each(|chunk| {
                if let Err(e) = self.process_sequence_chunk(chunk) {
                    eprintln!("Error processing sequence: {}", e);
                }
            });
    }

    pub fn add_hash(&mut self, hash: u64) {
        let mut minhash = self.minhash.lock().unwrap();
        minhash.add_hash(hash);
    }

    pub fn get_mins(&self) -> Vec<u64> {
        let minhash = self.minhash.lock().unwrap();
        minhash.to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    const TEST_FASTA: &str = "tests/testdata/index/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";

    #[test]
    fn test_proteome_index_creation() {
        let index = ProteomeIndex::new(7, 100, "protein", "raw");
        assert!(index.get_mins().is_empty());
    }

    #[test]
    fn test_proteome_index_with_dayhoff() {
        let index = ProteomeIndex::new(7, 100, "protein", "dayhoff");
        assert!(index.get_mins().is_empty());
    }

    #[test]
    fn test_proteome_index_with_hp() {
        let index = ProteomeIndex::new(7, 100, "protein", "hp");
        assert!(index.get_mins().is_empty());
    }

    #[test]
    fn test_process_protein_files() -> Result<()> {
        let index = ProteomeIndex::new(7, 100, "protein", "raw");
        let files = vec![PathBuf::from(TEST_FASTA)];
        index.process_protein_files(&files)?;

        // After processing, we should have some hashes
        let mins = index.get_mins();
        assert!(!mins.is_empty());

        Ok(())
    }

    #[test]
    fn test_invalid_sequence() {
        let index = ProteomeIndex::new(7, 100, "protein", "raw");
        let result = index.process_sequence_chunk("ACDEFGHIKLMNPQRSTVWY"); // Valid sequence
        assert!(result.is_ok());

        let result = index.process_sequence_chunk("ACDEFGHIKLMNPQRSTVWY1"); // Invalid character '1'
        assert!(result.is_err());
        let err = result.unwrap_err();
        assert!(err.to_string().contains("Invalid amino acid"));
    }

    #[test]
    fn test_ambiguous_amino_acids() {
        let index = ProteomeIndex::new(7, 100, "protein", "raw");

        // Test valid ambiguous amino acids
        let result = index.process_sequence_chunk("ACDEFXBZJ"); // X, B, Z, J are valid ambiguity codes
        assert!(result.is_ok());

        // When we get kmers containing ambiguous amino acids, they should be resolved
        let kmers = index.get_mins();
        assert!(!kmers.is_empty());
    }

    #[test]
    fn test_different_encodings() -> Result<()> {
        // Create indices with different encodings
        let raw_index = ProteomeIndex::new(7, 100, "protein", "raw");
        let dayhoff_index = ProteomeIndex::new(10, 100, "protein", "dayhoff");
        let hp_index = ProteomeIndex::new(15, 100, "protein", "hp");

        // Process the same sequence with each encoding
        let sequence = "ACDEFGHIKLMNPQRSTVWY";
        raw_index.process_sequence_chunk(sequence)?;
        dayhoff_index.process_sequence_chunk(sequence)?;
        hp_index.process_sequence_chunk(sequence)?;

        // Get the hashes from each index
        let raw_hashes = raw_index.get_mins();
        let dayhoff_hashes = dayhoff_index.get_mins();
        let hp_hashes = hp_index.get_mins();

        // Different encodings should produce different hash sets
        assert_ne!(raw_hashes, dayhoff_hashes);
        assert_ne!(raw_hashes, hp_hashes);
        assert_ne!(dayhoff_hashes, hp_hashes);

        Ok(())
    }

    #[test]
    fn test_kmer_retrieval() -> Result<()> {
        let index = ProteomeIndex::new(3, 100, "protein", "raw");

        // Add a simple sequence
        index.process_sequence_chunk("ACDEF")?;

        // Get the hashes
        let mins = index.get_mins();
        assert!(!mins.is_empty());

        // For each hash, we should be able to get back the original kmer
        for hash in mins {
            let kmer_info = index.get_kmers(hash);
            assert!(kmer_info.is_some());
            let (encoded, original) = kmer_info.unwrap();
            assert_eq!(encoded.len(), 3); // ksize = 3
            assert_eq!(original.len(), 3);
        }

        Ok(())
    }

    #[test]
    fn test_scaled_values() {
        // Test with different scaled values
        let index_100 = ProteomeIndex::new(7, 100, "protein", "raw");
        let index_1000 = ProteomeIndex::new(7, 1000, "protein", "raw");

        let sequence = "ACDEFGHIKLMNPQRSTVWY";
        index_100.process_sequence_chunk(sequence).unwrap();
        index_1000.process_sequence_chunk(sequence).unwrap();

        // Higher scaled value should result in fewer hashes
        let hashes_100 = index_100.get_mins();
        let hashes_1000 = index_1000.get_mins();
        assert!(hashes_100.len() >= hashes_1000.len());
    }
}
