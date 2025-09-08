use kmerseek::index::ProteomeIndex;
use tempfile::tempdir;

fn main() -> anyhow::Result<()> {
    let dir = tempdir()?;

    println!("=== ProteomeIndex Builder Pattern Demo ===\n");

    // Example 1: Using the builder pattern with explicit path
    println!("1. Creating index with builder pattern (explicit path):");
    let index1 = ProteomeIndex::builder()
        .path(dir.path().join("explicit.db"))
        .ksize(5)
        .scaled(1)
        .moltype("protein")
        .build()?;
    println!(
        "   ✓ Created index with ksize={}, scaled={}, moltype={}",
        index1.ksize(),
        index1.scaled(),
        index1.moltype(),
    );

    // Example 2: Using the builder pattern with auto filename generation
    println!("\n2. Creating index with builder pattern (auto filename):");
    let index2 = ProteomeIndex::builder()
        .path(dir.path().join("data.fasta"))
        .ksize(8)
        .scaled(10)
        .moltype("hp")
        .build_with_auto_filename()?;
    println!(
        "   ✓ Created index with ksize={}, scaled={}, moltype={}",
        index2.ksize(),
        index2.scaled(),
        index2.moltype(),
    );

    // Example 3: Using the builder pattern with raw sequence storage
    println!("\n3. Creating index with builder pattern (raw sequences):");
    let index3 = ProteomeIndex::builder()
        .path(dir.path().join("raw_sequences.db"))
        .ksize(6)
        .scaled(5)
        .moltype("protein")
        .store_raw_sequences(true)
        .build()?;
    println!(
        "   ✓ Created index with ksize={}, scaled={}, moltype={}, store_raw_sequences={}",
        index3.ksize(),
        index3.scaled(),
        index3.moltype(),
        index3.store_raw_sequences()
    );

    // Example 4: Using convenience methods
    println!("\n4. Using convenience methods:");
    let index4 = ProteomeIndex::new(dir.path().join("convenience.db"), 7, 2, "protein", false)?;
    println!(
        "   ✓ Created index with ksize={}, scaled={}, moltype={}",
        index4.ksize(),
        index4.scaled(),
        index4.moltype(),
    );

    // Example 5: Using convenience method with auto filename
    println!("\n5. Using convenience method with auto filename:");
    let index5 = ProteomeIndex::new_with_auto_filename(
        dir.path().join("auto_convenience.fasta"),
        9,
        3,
        "hp",
        true,
    )?;
    println!(
        "   ✓ Created index with ksize={}, scaled={}, moltype={}, store_raw_sequences={}",
        index5.ksize(),
        index5.scaled(),
        index5.moltype(),
        index5.store_raw_sequences()
    );

    println!("\n=== All examples completed successfully! ===");

    Ok(())
}
