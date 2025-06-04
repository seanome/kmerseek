use super::*;
use std::path::PathBuf;
use tempfile::tempdir;

const TEST_FASTA: &str =
    "tests/testdata/index/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";

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
    let dayhoff_index = ProteomeIndex::new(7, 100, "protein", "dayhoff");
    let hp_index = ProteomeIndex::new(7, 100, "protein", "hp");

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
