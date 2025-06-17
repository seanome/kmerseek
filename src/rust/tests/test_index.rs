use anyhow::Result;
use sourmash::encodings::HashFunctions;
use sourmash::signature::SigsTrait;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash_plugin_branchwater::utils::multicollection::SmallSignature;

use tempfile::tempdir;

use crate::index::ProteomeIndex;
use crate::tests::test_fixtures::{TEST_KMER, TEST_PROTEIN};
use crate::SEED;
use std::collections::HashMap;

#[test]
fn test_kmer_encoding() -> Result<()> {
    let dir = tempdir()?;

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("kmer_test.db"),
        11, // small ksize for test
        1,  // scaled
        "protein",
        SEED,
    )?;

    // Test k-mer encoding
    let (encoded_kmer, original_kmer) = index.encode_kmer(TEST_KMER);

    // The original k-mer should be exactly what we put in
    assert_eq!(original_kmer, TEST_KMER);

    // The encoded k-mer should be the same as the original k-mer
    assert_eq!(encoded_kmer, TEST_KMER);

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
        SEED,
    )?;

    // Test k-mer encoding
    let (encoded_kmer, original_kmer) = index.encode_kmer(TEST_KMER)?;

    // The original k-mer should be exactly what we put in
    assert_eq!(original_kmer, TEST_KMER);

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
        SEED,
    )?;

    // Test k-mer encoding
    let (encoded_kmer, original_kmer) = index.encode_kmer(TEST_KMER)?;

    // The original k-mer should be exactly what we put in
    assert_eq!(original_kmer, TEST_KMER);

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
        1, // ksize=1 for protein (will be multiplied by 3)
        1, // scaled=1 to capture all kmers
        "protein",
        SEED,
    )?;

    let sequence = TEST_PROTEIN;

    // Create a minhash sketch for protein
    let minhash = KmerMinHash::new(
        1, // scaled
        5, // ksize
        HashFunctions::Murmur64Protein,
        SEED, // seed
        true, // track_abundance
        0,    // num (use scaled instead)
    );

    let mut small_sig = SmallSignature {
        location: "test1".to_string(),
        name: "test1".to_string(),
        md5sum: "test1".to_string(),
        minhash: minhash,
    };
    small_sig.minhash.add_sequence(sequence.as_bytes(), true)?;

    // Process kmers
    let kmer_signature = index.process_protein_kmers(sequence, &small_sig)?;

    println!("{}", kmer_signature.signature.name);
    println!("{:?}", kmer_signature.kmer_infos.keys());
    let kmer_count = kmer_signature.kmer_infos.len();

    // Should have 17 kmers (length 21 - ksize 5 + 1)
    assert_eq!(kmer_count, 17);

    // Print all kmer infos for debugging
    println!("\nKmer Info Details:");
    println!("Hash\t\t\\tProtein\tOrig\tPos");
    println!("----------------------------------------");
    for (hash, kmer_info) in &kmer_signature.kmer_infos {
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
    for (hash, kmer_info) in &kmer_signature.kmer_infos {
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
fn test_process_protein_kmers_dayhoff() -> Result<()> {
    let dir = tempdir()?;

    let ksize = 5;

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("dayhoff_test.db"),
        ksize, // ksize=1 for protein (will be multiplied by 3)
        1,     // scaled=1 to capture all kmers
        "dayhoff",
        SEED,
    )?;

    let sequence = TEST_PROTEIN;

    // Create a minhash sketch for protein
    let minhash = KmerMinHash::new(
        1,     // scaled
        ksize, // ksize
        HashFunctions::Murmur64Dayhoff,
        SEED, // seed
        true, // track_abundance
        0,    // num (use scaled instead)
    );

    let mut small_sig = SmallSignature {
        location: "test1".to_string(),
        name: "test1".to_string(),
        md5sum: "test1".to_string(),
        minhash: minhash,
    };
    small_sig.minhash.add_sequence(sequence.as_bytes(), true)?;

    // Process kmers
    let kmer_signature = index.process_protein_kmers(sequence, &small_sig)?;

    println!("{}", kmer_signature.signature.name);
    let hashvals = kmer_signature.kmer_infos.keys().collect::<Vec<_>>();
    println!("{:?}", hashvals);
    let kmer_count = kmer_signature.kmer_infos.len();

    // Should have 17 kmers (length 21 - ksize 5 + 1)
    assert_eq!(kmer_count, 17);

    // Print all kmer infos for debugging
    println!("\nKmer Info Details:");
    println!("Hash\t\t\tDayhoff\tOrig\tPos");
    println!("----------------------------------------");
    for (hash, kmer_info) in &kmer_signature.kmer_infos {
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
    let actual_hashes: Vec<_> = kmer_signature.kmer_infos.keys().collect();
    assert_eq!(
        expected_hashes.len(),
        actual_hashes.len(),
        "Number of hashes mismatch: expected {}, got {}",
        expected_hashes.len(),
        actual_hashes.len()
    );

    for hash in expected_hashes {
        assert!(kmer_signature.kmer_infos.contains_key(hash), "Missing expected hash {}", hash);
    }

    // Verify each kmer info matches expected values
    for (hash, kmer_info) in &kmer_signature.kmer_infos {
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
            let positions = kmer_info
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
fn test_process_protein_kmers_hp() -> Result<()> {
    let dir = tempdir()?;

    let ksize = 5;

    // Create index with minimal parameters
    let index = ProteomeIndex::new(
        dir.path().join("hp_test.db"),
        ksize, // ksize=1 for protein (will be multiplied by 3)
        1,     // scaled=1 to capture all kmers
        "hp",
        SEED,
    )?;

    let sequence = TEST_PROTEIN;

    // Create a minhash sketch for protein
    let minhash = KmerMinHash::new(
        1,     // scaled
        ksize, // ksize
        HashFunctions::Murmur64Dayhoff,
        SEED, // seed
        true, // track_abundance
        0,    // num (use scaled instead)
    );

    let mut small_sig = SmallSignature {
        location: "test1".to_string(),
        name: "test1".to_string(),
        md5sum: "test1".to_string(),
        minhash: minhash,
    };
    small_sig.minhash.add_sequence(sequence.as_bytes(), true)?;

    // Process kmers
    let kmer_signature = index.process_protein_kmers(sequence, &small_sig)?;

    println!("{}", kmer_signature.signature.name);
    let hashvals = kmer_signature.kmer_infos.keys().collect::<Vec<_>>();
    println!("{:?}", hashvals);
    let kmer_count = kmer_signature.kmer_infos.len();

    // // Should have 14 kmers (length 21 - ksize 5 + 1), but a few duplicates
    assert_eq!(kmer_count, 14);

    // Print all kmer infos for debugging
    println!("\nKmer Info Details:");
    println!("Hash\t\t\tHP\tOrig\tPos");
    println!("----------------------------------------");
    for (hash, kmer_info) in &kmer_signature.kmer_infos {
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
    let actual_hashes: Vec<_> = kmer_signature.kmer_infos.keys().collect();
    assert_eq!(
        expected_hashes.len(),
        actual_hashes.len(),
        "Number of hashes mismatch: expected {}, got {}",
        expected_hashes.len(),
        actual_hashes.len()
    );

    for hash in expected_hashes {
        assert!(kmer_signature.kmer_infos.contains_key(hash), "Missing expected hash {}", hash);
    }

    // Verify each kmer info matches expected values
    for (hash, kmer_info) in &kmer_signature.kmer_infos {
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
            let positions = kmer_info
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
