use crate::uniprot::UniProtEntry;
use anyhow::Result;

use tempfile::tempdir;

use crate::index::ProteomeIndex;
use crate::tests::test_fixtures::{TEST_FASTA};

// #[test]
// fn test_proteome_index_creation() -> Result<()> {
//     let dir = tempdir()?;
//     let files = vec![PathBuf::from(TEST_FASTA)];

//     // Test with Dayhoff encoding
//     let mut index = ProteomeIndex::new(
//         dir.path().join("dayhoff.db"),
//         3, // ksize
//         1, // scaled
//         "dayhoff",
//     )?;
//     index.process_protein_files(&files)?;

//     // Test with HP encoding
//     let mut index = ProteomeIndex::new(
//         dir.path().join("hp.db"),
//         3, // ksize
//         1, // scaled
//         "hp",
//     )?;
//     index.process_protein_files(&files)?;

//     // Test with raw protein encoding
//     let mut index = ProteomeIndex::new(
//         dir.path().join("protein.db"),
//         3, // ksize
//         1, // scaled
//         "protein",
//     )?;
//     index.process_protein_files(&files)?;

//     // Test sequence validation
//     let result = index.process_sequence(TEST_PROTEIN_INVALID.as_bytes(), Default::default());
//     assert!(result.is_err());

//     // Test k-mer info retrieval
//     let kmer_info = index.get_kmer_info(42); // Use a known hash value
//     if let Some(info) = kmer_info {
//         for (md5sum, infos) in info {
//             println!("MD5: {}", md5sum);
//             for kmer in infos {
//                 println!("  K-mer: {}", kmer.original_kmer);
//                 for pos in kmer.positions {
//                     println!("    Position: {}", pos.position);
//                 }
//             }
//         }
//     }

//     // Test statistics
//     let mut index_100 = ProteomeIndex::new(
//         dir.path().join("stats_100.db"),
//         3,   // ksize
//         100, // scaled
//         "protein",
//     )?;
//     index_100.process_protein_files(&files)?;
//     index_100.compute_statistics()?;

//     let mut index_1000 = ProteomeIndex::new(
//         dir.path().join("stats_1000.db"),
//         3,    // ksize
//         1000, // scaled
//         "protein",
//     )?;
//     index_1000.process_protein_files(&files)?;
//     index_1000.compute_statistics()?;

//     // Compare statistics between different scaled values
//     for hash in index_100.get_hashes() {
//         if let (Some(stats_100), Some(stats_1000)) = (
//             index_100.get_kmer_stats(hash)?,
//             index_1000.get_kmer_stats(hash)?,
//         ) {
//             println!(
//                 "Hash {} - Scaled 100: {:.2} {:.2}, Scaled 1000: {:.2} {:.2}",
//                 hash, stats_100.idf, stats_100.frequency, stats_1000.idf, stats_1000.frequency
//             );
//         }
//     }

//     Ok(())
// }

#[test]
fn test_proteome_index_creation_raw() -> Result<()> {
    let dir = tempdir()?;
    let files = vec![TEST_FASTA];

    // Test with protein encoding
    let mut index = ProteomeIndex::new(
        dir.path().join("protein.db"),
        7, // ksize
        1, // scaled
        "protein",
    )?;
    index.process_protein_files(&files)?;
    index.compute_statistics()?;

    // Test k-mer info retrieval
    let kmer_info = index.get_kmer_info(42); // Use a known hash value
    if let Some(info) = kmer_info {
        for (md5sum, infos) in info {
            println!("MD5: {}", md5sum);
            for kmer in infos {
                println!("  K-mer: {}", kmer.original_kmer);
                for pos in kmer.positions {
                    println!("    Position: {}", pos.position);
                }
            }
        }
    }

    // Compare statistics between different scaled values
    for hash in index.get_hashes() {
        if let Some(stats) = index.get_kmer_stats(hash)? {
            println!(
                "Hash {} - IDF: {:.2}, Frequency: {:.2}",
                hash, stats.idf, stats.frequency,
            );
        }
    }

    Ok(())
}

// #[test]
// fn test_proteome_index_with_dayhoff() -> Result<()> {
//     let dir = tempdir()?;
//     let mut index = ProteomeIndex::new(
//         dir.path().join("dayhoff.db"),
//         7,   // ksize
//         100, // scaled
//         "dayhoff",
//     )?;
//     assert!(index.get_hashes().is_empty());
//     Ok(())
// }

// #[test]
// fn test_proteome_index_with_hp() -> Result<()> {
//     let dir = tempdir()?;
//     let mut index = ProteomeIndex::new(
//         dir.path().join("hp.db"),
//         7,   // ksize
//         100, // scaled
//         "hp",
//     )?;
//     assert!(index.get_hashes().is_empty());
//     Ok(())
// }

// #[test]
// fn test_proteome_index_with_protein() -> Result<()> {
//     let dir = tempdir()?;
//     let mut index = ProteomeIndex::new(
//         dir.path().join("protein.db"),
//         7,   // ksize
//         100, // scaled
//         "protein",
//     )?;
//     assert!(index.get_hashes().is_empty());
//     Ok(())
// }

#[test]
fn test_single_protein_addition() -> Result<()> {
    let dir = tempdir()?;

    // Create a simple test protein
    let test_protein = UniProtEntry {
        id: "TEST1".to_string(),
        sequence: "MVKVGVNG".to_string(), // Simple 8 amino acid sequence
        features: Vec::new(),
        ..Default::default()
    };

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("single_protein.db"),
        3, // small ksize for test
        1, // scaled
        "protein",
    )?;

    // Add the protein - clone it since process_sequence takes ownership
    index.process_sequence(test_protein.sequence.as_bytes(), test_protein.clone())?;

    // Verify the protein was added by checking if we can get k-mer info
    // Get the first hash from the index
    let hashes = index.get_hashes();
    assert!(
        !hashes.is_empty(),
        "Should have at least one hash after adding protein"
    );

    // Try to get k-mer info for the first hash
    if let Some(kmer_info) = index.get_kmer_info(hashes[0]) {
        assert!(!kmer_info.is_empty(), "Should have k-mer info for the hash");
    } else {
        panic!("Failed to get k-mer info for hash");
    }

    // Explicitly close the index
    index.close()?;

    Ok(())
}
