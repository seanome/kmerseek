use kmerseek::index::ProteomeIndex;
use kmerseek::signature::SEED;
use tempfile::tempdir;

fn main() -> anyhow::Result<()> {
    println!("Testing ProteomeIndex equivalence functionality...");

    let temp_dir = tempdir()?;
    let db_path1 = temp_dir.path().join("index1.db");
    let db_path2 = temp_dir.path().join("index2.db");

    // Create two indices with the same parameters
    println!("Creating two indices with identical parameters...");
    let index1 = ProteomeIndex::new(&db_path1, 5, 1, "protein", SEED, false)?;
    let index2 = ProteomeIndex::new(&db_path2, 5, 1, "protein", SEED, false)?;

    // Add the same protein sequences to both indices
    println!("Adding identical protein sequences to both indices...");
    let sequences = vec![
        ("ACDEFGHIKLMNPQRSTVWY", "protein1"),
        ("PLANTANDANIMALGENQMES", "protein2"),
        ("METHIONINELEUCINE", "protein3"),
    ];

    for (seq, name) in &sequences {
        let sig1 = index1.create_protein_signature(seq, name)?;
        let sig2 = index2.create_protein_signature(seq, name)?;

        index1.store_signatures(vec![sig1])?;
        index2.store_signatures(vec![sig2])?;
    }

    // Print statistics for both indices
    println!("\nIndex 1 statistics:");
    index1.print_stats();

    println!("\nIndex 2 statistics:");
    index2.print_stats();

    // Test equivalence
    println!("\nTesting equivalence...");
    let are_equivalent = index1.is_equivalent_to(&index2)?;
    println!("Indices are equivalent: {}", are_equivalent);

    // Create a third index with different parameters
    println!("\nCreating a third index with different parameters...");
    let index3 =
        ProteomeIndex::new(&temp_dir.path().join("index3.db"), 10, 1, "protein", SEED, false)?;
    let sig3 = index3.create_protein_signature("ACDEFGHIKLMNPQRSTVWY", "test_protein")?;
    index3.store_signatures(vec![sig3])?;

    // Test that different indices are not equivalent
    let are_equivalent_3 = index1.is_equivalent_to(&index3)?;
    println!("Index 1 equivalent to Index 3 (different ksize): {}", are_equivalent_3);
    assert!(!are_equivalent_3, "Indices with different parameters should not be equivalent");

    // Test with different sequences
    println!("\nCreating a fourth index with different sequences...");
    let index4 =
        ProteomeIndex::new(&temp_dir.path().join("index4.db"), 5, 1, "protein", SEED, false)?;
    let sig4 = index4.create_protein_signature("DIFFERENTSEQUENCE", "different")?;
    index4.store_signatures(vec![sig4])?;

    let are_equivalent_4 = index1.is_equivalent_to(&index4)?;
    println!("Index 1 and Index 4 are equivalent: {}", are_equivalent_4);

    println!("\nEquivalence testing completed successfully!");
    Ok(())
}
