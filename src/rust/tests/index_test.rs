use super::*;
use anyhow::Result;
use std::path::PathBuf;
use tempfile::tempdir;

use crate::index::ProteomeIndex;
use crate::uniprot::UniProtSequence;

const TEST_FASTA: &str =
    "tests/testdata/index/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";

#[test]
fn test_proteome_index_creation() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(7, 100, "protein", "raw", dir.path())?;
    let sequence = "ACDEF";
    let result = index.process_sequence(
        sequence,
        UniProtSequence {
            id: String::new(),
            accession: String::new(),
            sequence: sequence.to_string(),
            features: Vec::new(),
        },
    );
    assert!(result.is_ok());
    assert!(index.get_mins().is_empty());
    Ok(())
}

#[test]
fn test_proteome_index_with_dayhoff() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(7, 100, "protein", "dayhoff", dir.path())?;
    assert!(index.get_mins().is_empty());
    Ok(())
}

#[test]
fn test_proteome_index_with_hp() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(7, 100, "protein", "hp", dir.path())?;
    assert!(index.get_mins().is_empty());
    Ok(())
}

#[test]
fn test_process_protein_files() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(7, 100, "protein", "raw", dir.path())?;
    let files = vec![PathBuf::from(TEST_FASTA)];
    index.process_protein_files(&files)?;

    // After processing, we should have some hashes
    let mins = index.get_mins();
    assert!(!mins.is_empty());

    Ok(())
}

#[test]
fn test_invalid_sequence() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(7, 100, "protein", "raw", dir.path())?;
    let result = index.process_sequence(
        "ACDEFGHIKLMNPQRSTVWY",
        UniProtSequence {
            id: String::new(),
            accession: String::new(),
            sequence: "ACDEFGHIKLMNPQRSTVWY".to_string(),
            features: Vec::new(),
        },
    ); // Valid sequence
    assert!(result.is_ok());

    let result = index.process_sequence(
        "ACDEFGHIKLMNPQRSTVWY1",
        UniProtSequence {
            id: String::new(),
            accession: String::new(),
            sequence: "ACDEFGHIKLMNPQRSTVWY1".to_string(),
            features: Vec::new(),
        },
    ); // Invalid character '1'
    assert!(result.is_err());
    let err = result.unwrap_err();
    assert!(err.to_string().contains("Invalid amino acid"));
    Ok(())
}

#[test]
fn test_ambiguous_amino_acids() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(7, 100, "protein", "raw", dir.path())?;

    // Test valid ambiguous amino acids
    let result = index.process_sequence(
        "ACDEFXBZJ",
        UniProtSequence {
            id: String::new(),
            accession: String::new(),
            sequence: "ACDEFXBZJ".to_string(),
            features: Vec::new(),
        },
    ); // X, B, Z, J are valid ambiguity codes
    assert!(result.is_ok());

    // When we get kmers containing ambiguous amino acids, they should be resolved
    let kmers = index.get_mins();
    assert!(!kmers.is_empty());
    Ok(())
}

#[test]
fn test_different_encodings() -> Result<()> {
    let dir = tempdir().unwrap();
    // Create indices with different encodings
    let raw_index = ProteomeIndex::new(7, 100, "protein", "raw", dir.path().join("protein_db"))?;
    let dayhoff_index =
        ProteomeIndex::new(7, 100, "protein", "dayhoff", dir.path().join("dayhoff_db"))?;
    let hp_index = ProteomeIndex::new(7, 100, "protein", "hp", dir.path().join("hp_db"))?;

    // Process the same sequence with each encoding
    let sequence = "ACDEFGHIKLMNPQRSTVWY";
    let protein = UniProtSequence {
        id: String::new(),
        accession: String::new(),
        sequence: sequence.to_string(),
        features: Vec::new(),
    };
    raw_index.process_sequence(sequence, protein.clone())?;
    dayhoff_index.process_sequence(sequence, protein.clone())?;
    hp_index.process_sequence(sequence, protein)?;

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
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(3, 100, "protein", "raw", dir.path())?;

    // Add a simple sequence
    index.process_sequence(
        "ACDEF",
        UniProtSequence {
            id: String::new(),
            accession: String::new(),
            sequence: "ACDEF".to_string(),
            features: Vec::new(),
        },
    )?;

    // Get the hashes
    let mins = index.get_mins();
    assert!(!mins.is_empty());

    // For each hash, we should be able to get back the original kmer
    for hash in mins {
        let kmer_info = index.get_kmer_info(hash);
        assert!(kmer_info.is_some());
        let (encoded, original) = kmer_info.unwrap();
        assert_eq!(encoded.len(), 3); // ksize = 3
        assert!(!original.is_empty());
    }

    Ok(())
}

#[test]
fn test_scaled_values() -> Result<()> {
    let dir = tempdir().unwrap();
    // Test with different scaled values
    let index_100 = ProteomeIndex::new(7, 100, "protein", "raw", dir.path().join("index_100"))?;
    let index_1000 = ProteomeIndex::new(7, 1000, "protein", "raw", dir.path().join("index_1000"))?;

    let sequence = "ACDEFGHIKLMNPQRSTVWY";
    let protein = UniProtSequence {
        id: String::new(),
        accession: String::new(),
        sequence: sequence.to_string(),
        features: Vec::new(),
    };
    index_100.process_sequence(sequence, protein.clone())?;
    index_1000.process_sequence(sequence, protein)?;

    // Higher scaled value should result in fewer hashes
    let hashes_100 = index_100.get_mins();
    let hashes_1000 = index_1000.get_mins();
    assert!(hashes_100.len() >= hashes_1000.len());
    Ok(())
}
