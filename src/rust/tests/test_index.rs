use anyhow::Result;
use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};
use std::path::PathBuf;
use tempfile::tempdir;

use crate::index::ProteomeIndex;
use crate::tests::test_fixtures::*;
use crate::uniprot::UniProtEntry;

#[test]
fn test_proteome_index_creation() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(
        dir.path().to_path_buf(), // Convert &Path to PathBuf to satisfy AsRef<Path>
        15,                       // ksize
        1,                        // scaled
        "hp",                     // hydrophobic-polar encoding
    )?;
    let sequence = TEST_PROTEIN;
    let result = index.process_sequence(
        sequence,
        UniProtEntry {
            id: String::new(),
            accession: String::new(),
            sequence: sequence.to_string(),
            features: Vec::new(),
        },
    );
    println!("hashes: {:?}", index.get_hashes());
    println!("number of hashes: {:?}", index.get_hashes().len());
    assert!(result.is_ok());
    assert!(index.get_hashes().is_empty());
    Ok(())
}

#[test]
fn test_proteome_index_with_dayhoff() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(
        dir.path().to_path_buf(),
        7,
        100,
        HashFunctions::Murmur64Dayhoff,
        aa_to_dayhoff,
    )?;
    assert!(index.get_hashes().is_empty());
    Ok(())
}

#[test]
fn test_proteome_index_with_hp() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(
        dir.path().to_path_buf(),
        7,
        100,
        HashFunctions::Murmur64Hp,
        aa_to_hp,
    )?;
    assert!(index.get_hashes().is_empty());
    Ok(())
}

#[test]
fn test_process_protein_files() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(
        dir.path().to_path_buf(),
        7,
        100,
        HashFunctions::Murmur64Protein,
        |x| x,
    )?;
    let files = vec![PathBuf::from(TEST_FASTA)];
    index.process_protein_files(&files)?;

    // After processing, we should have some hashes
    let mins = index.get_hashes();
    assert!(!mins.is_empty());

    Ok(())
}

#[test]
fn test_invalid_sequence() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(
        dir.path().to_path_buf(),
        7,
        100,
        HashFunctions::Murmur64Protein,
        |x| x,
    )?;
    let result = index.process_sequence(
        "ACDEFGHIKLMNPQRSTVWY",
        UniProtEntry {
            id: String::new(),
            accession: String::new(),
            sequence: "ACDEFGHIKLMNPQRSTVWY".to_string(),
            features: Vec::new(),
        },
    ); // Valid sequence
    assert!(result.is_ok());

    let result = index.process_sequence(
        "ACDEFGHIKLMNPQRSTVWY1",
        UniProtEntry {
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
    let index = ProteomeIndex::new(
        dir.path().to_path_buf(),
        7,
        100,
        HashFunctions::Murmur64Protein,
        |x| x,
    )?;

    // Test valid ambiguous amino acids
    let result = index.process_sequence(
        "ACDEFXBZJ",
        UniProtEntry {
            id: String::new(),
            accession: String::new(),
            sequence: "ACDEFXBZJ".to_string(),
            features: Vec::new(),
        },
    ); // X, B, Z, J are valid ambiguity codes
    assert!(result.is_ok());

    // When we get kmers containing ambiguous amino acids, they should be resolved
    let kmers = index.get_hashes();
    assert!(!kmers.is_empty());
    Ok(())
}

#[test]
fn test_different_encodings() -> Result<()> {
    let dir = tempdir().unwrap();
    // Create indices with different encodings
    let raw_index = ProteomeIndex::new(
        dir.path().join("protein_db"),
        7,
        100,
        HashFunctions::Murmur64Protein,
        |x| x,
    )?;
    let dayhoff_index = ProteomeIndex::new(
        dir.path().join("dayhoff_db"),
        7,
        100,
        HashFunctions::Murmur64Dayhoff,
        aa_to_dayhoff,
    )?;
    let hp_index = ProteomeIndex::new(
        dir.path().join("hp_db"),
        7,
        100,
        HashFunctions::Murmur64Hp,
        aa_to_hp,
    )?;

    // Process the same sequence with each encoding
    let sequence = "ACDEFGHIKLMNPQRSTVWY";
    let protein = UniProtEntry {
        id: String::new(),
        accession: String::new(),
        sequence: sequence.to_string(),
        features: Vec::new(),
    };
    raw_index.process_sequence(sequence, protein.clone())?;
    dayhoff_index.process_sequence(sequence, protein.clone())?;
    hp_index.process_sequence(sequence, protein)?;

    // Get the hashes from each index
    let raw_hashes = raw_index.get_hashes();
    let dayhoff_hashes = dayhoff_index.get_hashes();
    let hp_hashes = hp_index.get_hashes();

    // Different encodings should produce different hash sets
    assert_ne!(raw_hashes, dayhoff_hashes);
    assert_ne!(raw_hashes, hp_hashes);
    assert_ne!(dayhoff_hashes, hp_hashes);

    Ok(())
}

#[test]
fn test_kmer_retrieval() -> Result<()> {
    let dir = tempdir().unwrap();
    let index = ProteomeIndex::new(
        dir.path().to_path_buf(),
        3,
        100,
        HashFunctions::Murmur64Protein,
        |x| x,
    )?;

    // Add a simple sequence
    index.process_sequence(
        "ACDEF",
        UniProtEntry {
            id: String::new(),
            accession: String::new(),
            sequence: "ACDEF".to_string(),
            features: Vec::new(),
        },
    )?;

    // Get the hashes
    let mins = index.get_hashes();
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
    let index_100 = ProteomeIndex::new(
        dir.path().join("index_100"),
        7,
        100,
        HashFunctions::Murmur64Protein,
        |x| x,
    )?;
    let index_1000 = ProteomeIndex::new(
        dir.path().join("index_1000"),
        7,
        1000,
        HashFunctions::Murmur64Protein,
        |x| x,
    )?;

    let sequence = "ACDEFGHIKLMNPQRSTVWY";
    let protein = UniProtEntry {
        id: String::new(),
        accession: String::new(),
        sequence: sequence.to_string(),
        features: Vec::new(),
    };
    index_100.process_sequence(sequence, protein.clone())?;
    index_1000.process_sequence(sequence, protein)?;

    // Higher scaled value should result in fewer hashes
    let hashes_100 = index_100.get_hashes();
    let hashes_1000 = index_1000.get_hashes();
    assert!(hashes_100.len() >= hashes_1000.len());
    Ok(())
}
