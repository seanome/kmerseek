use crate::uniprot::UniProtEntry;
use anyhow::Result;
use sourmash::cmd::ComputeParameters;
use sourmash::signature::Signature;
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

use tempfile::tempdir;

use crate::index::ProteomeIndex;
use crate::tests::test_fixtures::{TEST_FASTA, TEST_PROTEIN};

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
        sequence: TEST_PROTEIN.to_string(),
        features: Vec::new(),
        ..Default::default()
    };

    // Create index with minimal parameters
    let mut index = ProteomeIndex::new(
        dir.path().join("single_protein.db"),
        7, // small ksize for test
        1, // scaled
        "protein",
    )?;

    // Add the protein directly
    index.add_uniprot_entry(test_protein.clone())?;

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

#[test]
fn test_kmer_encoding() -> Result<()> {
    let dir = tempdir()?;

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("kmer_test.db"),
        11, // small ksize for test
        1,  // scaled
        "protein",
    )?;

    // Test k-mer encoding
    let test_kmer = b"LIVINGALIVE"; // 3-mer
    let (encoded_kmer, original_kmer) = index.encode_kmer(test_kmer);

    // The original k-mer should be exactly what we put in
    assert_eq!(original_kmer, "LIVINGALIVE");

    // The encoded k-mer should be the same as the original k-mer
    assert_eq!(encoded_kmer, "LIVINGALIVE");

    Ok(())
}

#[test]
fn test_kmer_encoding_dayhoff() -> Result<()> {
    let dir = tempdir()?;

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("kmer_test.db"),
        11, // small ksize for test
        1,  // scaled
        "dayhoff",
    )?;

    // Test k-mer encoding
    let test_kmer = b"LIVINGALIVE"; // 3-mer
    let (encoded_kmer, original_kmer) = index.encode_kmer(test_kmer);

    // The original k-mer should be exactly what we put in
    assert_eq!(original_kmer, "LIVINGALIVE");

    // The encoded k-mer should be dayhoff-encoded k-mer
    assert_eq!(encoded_kmer, "eeeecbbeeec");

    Ok(())
}

#[test]
fn test_kmer_encoding_hp() -> Result<()> {
    let dir = tempdir()?;

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("kmer_test.db"),
        11, // small ksize for test
        1,  // scaled
        "hp",
    )?;

    // Test k-mer encoding
    let test_kmer = b"LIVINGALIVE"; // 3-mer
    let (encoded_kmer, original_kmer) = index.encode_kmer(test_kmer);

    // The original k-mer should be exactly what we put in
    assert_eq!(original_kmer, "LIVINGALIVE");

    // The encoded k-mer should be hp-encoded k-mer
    assert_eq!(encoded_kmer, "hhhhphhhhhp");

    Ok(())
}

#[test]
fn test_process_protein_kmers() -> Result<()> {
    let dir = tempdir()?;

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("protein_test.db"),
        3, // small ksize for test
        1, // scaled=1 to capture all kmers
        "protein",
    )?;

    let sequence = TEST_PROTEIN;

    // Create test signature
    let mut sig = Signature::from_params(
        &ComputeParameters::builder()
            .ksizes(vec![3])
            .scaled(1)
            .protein(true)
            .num_hashes(0)
            .build(),
    );
    sig.set_name("test1");

    let small_sig = SmallSignature {
        location: sig.filename().clone(),
        name: sig.name().clone().unwrap(),
        md5sum: sig.md5sum().to_string(),
        minhash: sig.minhash().unwrap().clone(),
    };

    // Process kmers
    let kmer_signature = index.process_protein_kmers(sequence, &small_sig)?;

    println!("{}", kmer_signature.signature.name);
    println!("{:?}", kmer_signature.kmer_infos.keys());
    for (hash, kmer_info) in kmer_signature.kmer_infos {
        println!("Hash: {}", hash);
        println!("Kmer info: {:?}", kmer_info);
    }

    // Should have 4 kmers (length 6 - ksize 3 + 1)
    assert_eq!(kmer_signature.kmer_infos.len(), 4);

    Ok(())
}
