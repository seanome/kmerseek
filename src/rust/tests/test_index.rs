use anyhow::Result;
use sourmash::signature::SigsTrait;

use tempfile::tempdir;

use crate::index::ProteomeIndex;
use crate::signature::ProteinSignature;
use crate::tests::test_fixtures::{TEST_FASTA_CONTENT, TEST_FASTA_GZ, TEST_PROTEIN};
use crate::SEED;
use std::collections::HashMap;

/// Keeping the tests for ProteomeIndex in a separate file because they're more like integration tests
/// than unit tests with all the moltype testing. Also, it's a lot of tests!

#[test]
fn test_process_kmers_moltype_protein() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "protein";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("dayhoff_test.db"),
        protein_ksize, // protein ksize
        1,             // scaled=1 to capture all kmers
        moltype,
        SEED,
    )?;

    let sequence = TEST_PROTEIN;

    // Create a protein signature
    let mut protein_sig = ProteinSignature::new(
        "test_protein",
        protein_ksize,
        1, // scaled
        moltype,
        SEED,
    )?;

    // Add the sequence
    protein_sig.add_protein(sequence.as_bytes())?;
    println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

    // Process kmers
    index.process_kmers(sequence, &mut protein_sig)?;

    println!("{}", protein_sig.signature().name);
    println!("{:?}", protein_sig.kmer_infos().keys());
    let kmer_count = protein_sig.kmer_infos().len();

    // Should have 17 kmers (length 21 - ksize 5 + 1)
    assert_eq!(kmer_count, 17);

    // Print all kmer infos for debugging
    println!("\nKmer Info Details:");
    println!("Hash\t\t\\tProtein\tOrig\tPos");
    println!("----------------------------------------");
    for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
        for (original_kmer, positions) in &kmer_info.original_kmer_to_position {
            println!("{}\t{}\t{}\t{:?}", hash, kmer_info.encoded_kmer, original_kmer, positions);
        }
    }
    println!("----------------------------------------\n");

    // Create expected hashmap of kmer info
    let raw_data = [
        // Hash               Original  Position
        (2140811952770908281, ("GENQM", [14])),
        (4381446250900425522, ("ENQME", [15])),
        (5798339600059429290, ("DANIM", [7])),
        (7681438632487987439, ("ANIMA", [8])),
        (12896310179337320481, ("LANTA", [1])),
        (2542642819229379552, ("NTAND", [3])),
        (11965201914550078735, ("TANDA", [4])),
        (5893010049374798421, ("PLANT", [0])),
        (110005740849399217, ("NDANI", [6])),
        (3791883307084689782, ("LGENQ", [13])),
        (14610011480386804007, ("ALGEN", [12])),
        (6941015416212662126, ("ANTAN", [2])),
        (12636705882654324958, ("NQMES", [16])),
        (11154024130290913208, ("IMALG", [10])),
        (1225702037828834387, ("MALGE", [11])),
        (12274863873578753245, ("NIMAL", [9])),
        (13616372540306653069, ("ANDAN", [5])),
    ];

    // Convert raw data into the required format with proper string types
    let expected_kmers: HashMap<_, _> = raw_data
        .into_iter()
        .map(|(hash, (kmer, positions))| {
            let mut original_map = HashMap::new();
            original_map.insert(kmer.to_string(), positions.to_vec());
            (hash, (kmer.to_string(), original_map))
        })
        .collect();

    // Verify each kmer info matches expected values
    for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
        let (expected_kmer, expected_positions) =
            expected_kmers.get(hash).expect(&format!("Missing expected hash {}", hash));

        // Verify the k-mer
        assert_eq!(
            &kmer_info.encoded_kmer, expected_kmer,
            "K-mer mismatch for hash {}: expected {}, got {}",
            hash, expected_kmer, kmer_info.encoded_kmer
        );

        // Verify positions for each original k-mer
        for (original_kmer, positions) in &kmer_info.original_kmer_to_position {
            let expected_pos = expected_positions
                .get(original_kmer)
                .expect(&format!("Missing positions for k-mer {}", original_kmer));
            assert_eq!(
                positions, expected_pos,
                "Position mismatch for k-mer {}: expected {:?}, got {:?}",
                original_kmer, expected_pos, positions
            );
        }
    }

    Ok(())
}

#[test]
fn test_process_kmers_moltype_dayhoff() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("dayhoff_test.db"),
        protein_ksize, // protein ksize
        1,             // scaled=1 to capture all kmers
        "dayhoff",
        SEED,
    )?;

    let sequence = TEST_PROTEIN;

    // Create a protein signature
    let mut protein_sig = ProteinSignature::new(
        "test_protein",
        protein_ksize,
        1, // scaled
        "dayhoff",
        SEED,
    )?;

    // Add the sequence
    protein_sig.add_protein(sequence.as_bytes())?;
    println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

    // Process kmers
    index.process_kmers(sequence, &mut protein_sig)?;

    println!("{}", protein_sig.signature().name);
    let hashvals = protein_sig.kmer_infos().keys().collect::<Vec<_>>();
    println!("{:?}", hashvals);
    let kmer_count = protein_sig.kmer_infos().len();

    // Should have 17 kmers (length 21 - ksize 5 + 1)
    assert_eq!(kmer_count, 17);

    // Print all kmer infos for debugging
    println!("\nKmer Info Details:");
    println!("Hash\t\t\tDayhoff\tOrig\tPos");
    println!("----------------------------------------");
    for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
        for (original_kmer, positions) in &kmer_info.original_kmer_to_position {
            println!("{}\t{}\t{}\t{:?}", hash, kmer_info.encoded_kmer, original_kmer, positions);
        }
    }
    println!("----------------------------------------\n");

    // Define the raw data without string conversions
    let raw_data = [
        (17444159595263538048, ("ceebe", "NIMAL", [9])),
        (2945598193614695589, ("cccec", "ENQME", [15])),
        (4548757849819812604, ("bbccb", "TANDA", [4])),
        (6463872878592804545, ("ebccc", "LGENQ", [13])),
        (4030406117949362159, ("cbcee", "DANIM", [7])),
        (7014407397606522347, ("ebcbb", "LANTA", [1])),
        (5045972850709227854, ("bebcb", "PLANT", [0])),
        (11417072151730334367, ("bcbbc", "ANTAN", [2])),
        (13574922562423607435, ("bceeb", "ANIMA", [8])),
        (15050500149255106627, ("bccce", "GENQM", [14])),
        (5430883729707969951, ("eebeb", "IMALG", [10])),
        (13894194422852851851, ("bebcc", "ALGEN", [12])),
        (9604281550621775790, ("bccbc", "ANDAN", [5])),
        (6161374941338912337, ("ccecb", "NQMES", [16])),
        (655307631517862365, ("ccbce", "NDANI", [6])),
        (360995089333906261, ("ebebc", "MALGE", [11])),
        (15056713696431004031, ("cbbcc", "NTAND", [3])),
    ];

    // Convert raw strings to owned types and create the HashMap
    let expected_kmers: HashMap<_, _> = raw_data
        .into_iter()
        .map(|(hash, (encoded, original, positions))| {
            let mut original_map = HashMap::new();
            original_map.insert(original.to_string(), positions.to_vec());
            (hash, (encoded.to_string(), original_map))
        })
        .collect();

    // Verify all expected hashes are present
    let expected_hashes: Vec<_> = expected_kmers.keys().collect();
    let actual_hashes: Vec<_> = protein_sig.kmer_infos().keys().collect();
    assert_eq!(
        expected_hashes.len(),
        actual_hashes.len(),
        "Number of hashes mismatch: expected {}, got {}",
        expected_hashes.len(),
        actual_hashes.len()
    );

    for hash in expected_hashes {
        assert!(protein_sig.kmer_infos().contains_key(hash), "Missing expected hash {}", hash);
    }

    // Verify each kmer info matches expected values
    for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
        let (expected_encoded, expected_originals) =
            expected_kmers.get(hash).expect(&format!("Missing expected hash {}", hash));

        // Verify the encoded k-mer
        assert_eq!(
            &kmer_info.encoded_kmer, expected_encoded,
            "Encoded k-mer mismatch for hash {}: expected {}, got {}",
            hash, expected_encoded, kmer_info.encoded_kmer
        );

        // Verify each original k-mer and its positions
        for (original_kmer, expected_positions) in expected_originals {
            let positions = protein_sig
                .kmer_infos()
                .get(hash)
                .unwrap()
                .original_kmer_to_position
                .get(original_kmer)
                .expect(&format!("Missing positions for k-mer {}", original_kmer));
            assert_eq!(
                positions, expected_positions,
                "Position mismatch for k-mer {}: expected {:?}, got {:?}",
                original_kmer, expected_positions, positions
            );
        }
    }

    Ok(())
}

#[test]
fn test_process_kmers_moltype_hp() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "hp";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("dayhoff_test.db"),
        protein_ksize, // protein ksize
        1,             // scaled=1 to capture all kmers
        moltype,
        SEED,
    )?;

    let sequence = TEST_PROTEIN;

    // Create a protein signature
    let mut protein_sig = ProteinSignature::new(
        "test_protein",
        protein_ksize,
        1, // scaled
        moltype,
        SEED,
    )?;

    // Add the sequence
    protein_sig.add_protein(sequence.as_bytes())?;
    println!("small_sig.minhash.to_vec(): {:?}", protein_sig.signature().minhash.to_vec());

    // Process kmers
    index.process_kmers(sequence, &mut protein_sig)?;

    println!("{}", protein_sig.signature().name);
    let hashvals = protein_sig.kmer_infos().keys().collect::<Vec<_>>();
    println!("{:?}", hashvals);
    let kmer_count = protein_sig.kmer_infos().len();

    // // Should have 14 kmers (length 21 - ksize 5 + 1), but a few duplicates
    assert_eq!(kmer_count, 14);

    // Print all kmer infos for debugging
    println!("\nKmer Info Details:");
    println!("Hash\t\t\tHP\tOrig\tPos");
    println!("----------------------------------------");
    for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
        for (original_kmer, positions) in &kmer_info.original_kmer_to_position {
            println!("{}\t{}\t{}\t{:?}", hash, kmer_info.encoded_kmer, original_kmer, positions);
        }
    }
    println!("----------------------------------------\n");

    // Define test data in a more readable format
    let kmer_data: HashMap<u64, (String, HashMap<String, Vec<usize>>)> = vec![
        // Single k-mer cases
        (17248460043117039725, ("hhhhp", vec!["MALGE"], vec![11])),
        (5673218808929106268, ("phhhh", vec!["NIMAL"], vec![9])),
        (16969835101383990681, ("hhpph", vec!["LANTA"], vec![1])),
        (7345312524621807974, ("pphph", vec!["NDANI"], vec![6])),
        (16370543730027378051, ("phpph", vec!["TANDA"], vec![4])),
        (3278382041688965244, ("hphhh", vec!["ANIMA"], vec![8])),
        (8541583772724823208, ("hhhhh", vec!["IMALG"], vec![10])),
        (16158526221854164806, ("hppph", vec!["GENQM"], vec![14])),
        (11553019557737058697, ("hhppp", vec!["LGENQ"], vec![13])),
        (9081059129327932468, ("ppphp", vec!["ENQME"], vec![15])),
        (2863220259252354754, ("phphh", vec!["DANIM"], vec![7])),
        // Multiple original protein k-mer sequences mapping to same HP encoding
        (4230974618842309829, ("hhhpp", vec!["PLANT", "ALGEN"], vec![0, 12])),
        (13058023948041027181, ("pphpp", vec!["NQMES", "NTAND"], vec![16, 3])),
        (4144736064335623701, ("hpphp", vec!["ANDAN", "ANTAN"], vec![5, 2])),
    ]
    .into_iter()
    .map(|(hash, (encoded, originals, positions))| {
        let mut original_map = HashMap::new();
        for (i, orig) in originals.into_iter().enumerate() {
            original_map.insert(orig.to_string(), vec![positions[i]]);
        }
        (hash, (encoded.to_string(), original_map))
    })
    .collect();

    // Verify all expected hashes are present
    let expected_hashes: Vec<_> = kmer_data.keys().collect();
    let actual_hashes: Vec<_> = protein_sig.kmer_infos().keys().collect();
    assert_eq!(
        expected_hashes.len(),
        actual_hashes.len(),
        "Number of hashes mismatch: expected {}, got {}",
        expected_hashes.len(),
        actual_hashes.len()
    );

    for hash in expected_hashes {
        assert!(protein_sig.kmer_infos().contains_key(hash), "Missing expected hash {}", hash);
    }

    // Verify each kmer info matches expected values
    for (hash, kmer_info) in protein_sig.kmer_infos().iter() {
        let (expected_encoded, expected_originals) =
            kmer_data.get(hash).expect(&format!("Missing expected hash {}", hash));

        // Verify the encoded k-mer
        assert_eq!(
            &kmer_info.encoded_kmer, expected_encoded,
            "Encoded k-mer mismatch for hash {}: expected {}, got {}",
            hash, expected_encoded, kmer_info.encoded_kmer
        );

        // Verify each original k-mer and its positions
        for (original_kmer, expected_positions) in expected_originals {
            let positions = protein_sig
                .kmer_infos()
                .get(hash)
                .unwrap()
                .original_kmer_to_position
                .get(original_kmer)
                .expect(&format!("Missing original k-mer {} for hash {}", original_kmer, hash));
            assert_eq!(
                positions, expected_positions,
                "Position mismatch for k-mer {}: expected {:?}, got {:?}",
                original_kmer, expected_positions, positions
            );
        }

        // Verify no unexpected original k-mers
        assert_eq!(
            kmer_info.original_kmer_to_position.len(),
            expected_originals.len(),
            "Expected {} original k-mer mappings for hash {}, got {}",
            expected_originals.len(),
            hash,
            kmer_info.original_kmer_to_position.len()
        );
    }

    Ok(())
}

#[test]
fn test_create_protein_signature_moltype_protein() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "protein";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("protein_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    let sequence = TEST_PROTEIN;
    let name = "test_protein";

    // Add the protein sequence to the index and get the signature
    let signature = index.create_protein_signature(sequence, name)?;

    // Verify the signature has the expected number of k-mers
    assert_eq!(signature.kmer_infos().len(), 17, "Expected 17 k-mers for the test protein");

    // Verify some specific k-mers are present
    let expected_hash = 5893010049374798421; // Hash for "PLANT"
    assert!(
        signature.kmer_infos().contains_key(&expected_hash),
        "Expected k-mer hash {} to be present",
        expected_hash
    );

    // Store the signature in the index
    index.store_signatures(vec![signature])?;

    // Verify the signature was added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 1, "Expected 1 signature to be stored");
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        assert!(combined_minhash.size() == 17, "Combined minhash should contain 17 hashes");
    }

    Ok(())
}

#[test]
fn test_create_protein_signature_moltype_dayhoff() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "dayhoff";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("dayhoff_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    let sequence = TEST_PROTEIN;
    let name = "test_protein";

    // Add the protein sequence to the index and get the signature
    let signature = index.create_protein_signature(sequence, name)?;

    // Verify the signature has the expected number of k-mers
    assert_eq!(signature.kmer_infos().len(), 17, "Expected 17 k-mers for the test protein");

    // Verify some specific k-mers are present
    let expected_hash = 5045972850709227854; // Hash for "PLANT" in Dayhoff encoding ("bebcb")
    assert!(
        signature.kmer_infos().contains_key(&expected_hash),
        "Expected k-mer hash {} to be present",
        expected_hash
    );

    // Store the signature in the index
    index.store_signatures(vec![signature])?;

    // Verify the signature was added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 1, "Expected 1 signature to be stored");
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        assert!(combined_minhash.size() == 17, "Combined minhash should contain hashes");
    }

    Ok(())
}

#[test]
fn test_create_protein_signature_moltype_hp() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "hp";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("hp_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    let sequence = TEST_PROTEIN;
    let name = "test_protein";

    // Add the protein sequence to the index and get the signature
    let signature = index.create_protein_signature(sequence, name)?;

    // Verify the signature has the expected number of k-mers
    assert_eq!(signature.kmer_infos().len(), 14, "Expected 14 k-mers for the test protein");

    // Verify some specific k-mers are present
    let expected_hash = 4230974618842309829; // Hash for "PLANT" in HP encoding ("hhhpp")
    assert!(
        signature.kmer_infos().contains_key(&expected_hash),
        "Expected k-mer hash {} to be present",
        expected_hash
    );

    // Store the signature in the index
    index.store_signatures(vec![signature])?;

    // Verify the signature was added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 1, "Expected 1 signature to be stored");
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        assert!(combined_minhash.size() == 14, "Combined minhash should contain 14 hashes");
    }

    Ok(())
}

#[test]
fn test_process_fasta_moltype_protein() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "protein";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("fasta_protein_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    // Create a temporary FASTA file for testing with distinct sequences
    let fasta_content = TEST_FASTA_CONTENT;
    let fasta_path = dir.path().join("test.fasta");
    std::fs::write(&fasta_path, fasta_content)?;

    // Process the FASTA file
    index.process_fasta(&fasta_path)?;

    // Verify the signatures were added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 2, "Expected 2 signatures to be stored");

        // Verify each signature has the expected number of k-mers
        for (md5sum, stored_signature) in signatures.iter() {
            if md5sum == "f7661cd829e75c0d" {
                assert!(
                    stored_signature.kmer_infos().len() == 7,
                    "LIVINGALIVE should have 7 protein 5-mers"
                );
            } else if md5sum == "7641839ad508ab8" {
                assert!(
                    stored_signature.kmer_infos().len() == 17,
                    "PLANTANDANIMALGENQMES should have 17 protein 5-mers"
                );
            } else {
                println!("md5sum: {}", md5sum);
                println!("Name: {}", stored_signature.signature().name);
                println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                assert!(false, "Unknown md5sum: {}", md5sum);
            }
        }
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        println!("combined_minhash.size(): {}", combined_minhash.size());
        assert!(combined_minhash.size() == 24, "Combined minhash should contain 24 hashes");
    }

    Ok(())
}

#[test]
fn test_process_fasta_moltype_dayhoff() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "dayhoff";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("fasta_dayhoff_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    // Create a temporary FASTA file for testing with distinct sequences
    let fasta_content = TEST_FASTA_CONTENT;
    let fasta_path = dir.path().join("test.fasta");
    std::fs::write(&fasta_path, fasta_content)?;

    // Process the FASTA file
    index.process_fasta(&fasta_path)?;

    // Verify the signatures were added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 2, "Expected 2 signatures to be stored");

        // Verify each signature has the expected number of k-mers
        for (md5sum, stored_signature) in signatures.iter() {
            if md5sum == "a963d06839b6d6a9" {
                assert!(
                    stored_signature.kmer_infos().len() == 7,
                    "LIVINGALIVE should have 7 dayhoff 5-mers"
                );
            } else if md5sum == "84d7545d531dcf51" {
                assert!(
                    stored_signature.kmer_infos().len() == 17,
                    "PLANTANDANIMALGENQMES should have 17 dayhoff 5-mers"
                );
            } else {
                println!("md5sum: {}", md5sum);
                println!("Name: {}", stored_signature.signature().name);
                println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                assert!(false, "Unknown md5sum: {}", md5sum);
            }
        }
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        println!("combined_minhash.size(): {}", combined_minhash.size());
        assert!(combined_minhash.size() == 24, "Combined minhash should contain 24 hashes");
    }

    Ok(())
}

#[test]
fn test_process_fasta_moltype_hp() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "hp";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("fasta_hp_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    // Create a temporary FASTA file for testing with distinct sequences
    let fasta_content = TEST_FASTA_CONTENT;
    let fasta_path = dir.path().join("test.fasta");
    std::fs::write(&fasta_path, fasta_content)?;

    // Process the FASTA file
    index.process_fasta(&fasta_path)?;

    // Verify the signatures were added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 2, "Expected 2 signatures to be stored");

        // Verify each signature has the expected number of k-mers
        for (md5sum, stored_signature) in signatures.iter() {
            if md5sum == "24ca8d939672666b" {
                assert!(
                    stored_signature.kmer_infos().len() == 6,
                    "LIVINGALIVE should have 6 hp 5-mers"
                );
            } else if md5sum == "668d7173d661287b" {
                assert!(
                    stored_signature.kmer_infos().len() == 14,
                    "PLANTANDANIMALGENQMES should have 14 hp 5-mers"
                );
            } else {
                println!("md5sum: {}", md5sum);
                println!("Name: {}", stored_signature.signature().name);
                println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
                assert!(false, "Unknown md5sum: {}", md5sum);
            }
        }
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        println!("combined_minhash.size(): {}", combined_minhash.size());
        assert!(combined_minhash.size() == 16, "Combined minhash should contain 16 hashes");
    }

    Ok(())
}

#[test]
fn test_process_fasta_gz_moltype_protein() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "protein";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("fasta_gz_protein_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    // Process the FASTA file
    index.process_fasta(TEST_FASTA_GZ)?;

    // Verify the signatures were added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 25, "Expected 25 signatures to be stored");

        // Check a few signatures
        for (md5sum, stored_signature) in signatures.iter() {
            println!("\n---\nmd5sum: {}", md5sum);
            println!("Name: {}", stored_signature.signature().name);
            println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
            if md5sum == "4d565dee9c8de9db" {
                assert!(
                    stored_signature.kmer_infos().len() == 474,
                    "sp|O43236|SEPT4_HUMAN should have 474 protein 5-mers"
                );
            }
            if md5sum == "4da1f84ad8be618e" {
                assert!(
                    stored_signature.kmer_infos().len() == 235,
                    "sp|P10415|BCL2_HUMAN should have 235 protein 5-mers"
                );
            }
        }
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        println!("combined_minhash.size(): {}", combined_minhash.size());
        assert!(
            combined_minhash.size() == 9049,
            "Combined minhash should contain 9049 protein 5-mer hashes"
        );
    }

    Ok(())
}

#[test]
fn test_process_fasta_gz_moltype_dayhoff() -> Result<()> {
    let dir = tempdir()?;

    let protein_ksize = 5;
    let moltype = "dayhoff";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("fasta_gz_dayhoff_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    // Process the FASTA file
    index.process_fasta(TEST_FASTA_GZ)?;

    // Verify the signatures were added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 25, "Expected 25 signatures to be stored");

        // Check a few signatures
        for (md5sum, stored_signature) in signatures.iter() {
            println!("\n---\nmd5sum: {}", md5sum);
            println!("Name: {}", stored_signature.signature().name);
            println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
            if md5sum == "fc27dcd533217385" {
                assert!(
                    stored_signature.kmer_infos().len() == 433,
                    "sp|O43236|SEPT4_HUMAN should have 433 dayhoff 5-mers"
                );
            }
            if md5sum == "3206706fa14185e7" {
                assert!(
                    stored_signature.kmer_infos().len() == 204,
                    "sp|P10415|BCL2_HUMAN should have 204 dayhoff 5-mers"
                );
            }
        }
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        println!("combined_minhash.size(): {}", combined_minhash.size());
        assert!(
            combined_minhash.size() == 2730,
            "Combined minhash should contain 2730 dayhoff 5-mer hashes"
        );
    }

    Ok(())
}

#[test]
fn test_process_fasta_gz_moltype_hp() -> Result<()> {
    let dir = tempdir()?;

    // Need a higher k-mer size because otherwise the binary space of 5-mers is saturated (other tests use 5-mers)
    // If k=5, then all signatures have ~32 (=2^5) 5-mers:
    // - Not unique -> each md5sum is identical -> "25 signatures" fails
    // - "Combined minhash" only contains 32 hashes, which isn't an interesting test
    // If k=12, then all signatures have ~4096 (=2^12) 12-mers:
    // - Unique -> "25 signatures" passes
    // - "Combined minhash" has an upper bound of 4096 hashes, which is a lot more interesting
    let protein_ksize = 12;
    let moltype = "hp";

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("fasta_gz_hp_test.db"),
        protein_ksize,
        1, // scaled=1 to capture all kmers for testing
        moltype,
        SEED,
    )?;

    // Process the FASTA file
    index.process_fasta(TEST_FASTA_GZ)?;

    // Verify the signatures were added to the signatures map
    {
        let signatures = index.get_signatures().lock().unwrap();
        assert_eq!(signatures.len(), 25, "Expected 25 signatures to be stored");

        // Check a few signatures
        for (md5sum, stored_signature) in signatures.iter() {
            println!("\n---\nmd5sum: {}", md5sum);
            println!("Name: {}", stored_signature.signature().name);
            println!("Len of Kmer infos: {}", stored_signature.kmer_infos().len());
            if md5sum == "38ffedf9d3ec7cec" {
                assert!(
                    stored_signature.kmer_infos().len() == 452,
                    "sp|O43236|SEPT4_HUMAN should have 452 hp 12-mers"
                );
            }
            if md5sum == "204716e4d80eb350" {
                assert!(
                    stored_signature.kmer_infos().len() == 220,
                    "sp|P10415|BCL2_HUMAN should have 220 hp 12-mers"
                );
            }
        }
    }

    // Verify the combined minhash was updated
    {
        let combined_minhash = index.get_combined_minhash().lock().unwrap();
        println!("combined_minhash.size(): {}", combined_minhash.size());
        assert!(
            combined_minhash.size() == 3549,
            "Combined minhash should contain 3549 hp 12-mer hashes"
        );
    }

    Ok(())
}
