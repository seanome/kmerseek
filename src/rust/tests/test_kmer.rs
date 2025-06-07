use anyhow::Result;
use sourmash::encodings::{aa_to_dayhoff, aa_to_hp, HashFunctions};
use std::path::PathBuf;
use tempfile::tempdir;

use crate::kmer::KmerPosition;
use crate::tests::test_fixtures::*;
use crate::uniprot::ProteinFeature;

#[test]
fn test_kmer_position() -> Result<()> {
    use crate::uniprot::{ProteinFeature, UniProtEntry};

    let protein = UniProtEntry {
        id: "test_protein".to_string(),
        accession: "TEST123".to_string(),
        sequence: TEST_PROTEIN.to_string(),
        features: vec![
            ProteinFeature {
                feature_type: "Feature 1".to_string(),
                start: 1, // 1-based position
                end: 5,   // covers PLANT
                description: "Test feature 1".to_string(),
            },
            ProteinFeature {
                feature_type: "Feature 2".to_string(),
                start: 4, // 1-based position
                end: 8,   // covers TAND
                description: "Test feature 2".to_string(),
            },
        ],
    };

    let kmer_pos = KmerPosition {
        protein,
        position: 2, // 0-based position, points to 'ANT'
    };

    let overlapping = kmer_pos.overlapping_features(3); // kmer length of 3
    assert_eq!(overlapping.len(), 2); // Both features overlap with 'ANT'
    assert_eq!(overlapping[0].name, "Feature 1");
    assert_eq!(overlapping[1].name, "Feature 2");

    Ok(())
}
