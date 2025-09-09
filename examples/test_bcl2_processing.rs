use kmerseek::index::ProteomeIndex;
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    println!("Processing BCL2 FASTA file with automatic filename generation...");

    // Define the FASTA file path
    let fasta_path = PathBuf::from("tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz");

    // Ensure the FASTA file exists
    if !fasta_path.exists() {
        eprintln!("BCL2 FASTA file not found at {:?}", fasta_path);
        eprintln!("Please ensure the test data is available.");
        return Ok(());
    }

    // Test different parameter combinations - run multiple times for better profiling
    let test_cases = vec![
        (7, 1, "protein", "BCL2 with protein encoding, k=7, scaled=1"),
        (10, 1, "dayhoff", "BCL2 with dayhoff encoding, k=10, scaled=1"),
        (14, 1, "hp", "BCL2 with hp encoding, k=14, scaled=1"),
    ];

    // Run multiple iterations to get more profiling data
    let iterations = 5;

    for iteration in 0..iterations {
        println!("\n=== ITERATION {} ===", iteration + 1);

        for (ksize, scaled, moltype, description) in &test_cases {
            println!("\n--- {} ---", description);

            // Create index with automatic filename generation
            let auto_index = ProteomeIndex::new_with_auto_filename(
                &fasta_path,
                *ksize,
                *scaled,
                moltype,
                false,
            )?;

            // Show the generated filename
            let generated_filename = auto_index.generate_filename(
                "bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz",
            );
            println!("Generated filename: {}", generated_filename);

            // Process the FASTA file
            println!("Processing FASTA file...");
            auto_index.process_fasta(&fasta_path, 10, 1000)?;

            // Print statistics
            println!("Index statistics:");
            auto_index.print_stats();

            // Show some example signatures
            let signatures = auto_index.get_signatures();
            let sig_map = signatures;
            println!("First few signatures:");
            for (i, entry) in sig_map.iter().take(3).enumerate() {
                let md5 = entry.key();
                let sig = entry.value();
                println!("  {}. {} (md5: {})", i + 1, sig.signature().name, md5);
            }
        }
    }

    println!("\nBCL2 processing examples completed successfully!");
    println!("\nNote: The automatic filename generation creates predictable filenames based on:");
    println!("  - Original FASTA filename");
    println!("  - Molecular type (hp, protein, dayhoff)");
    println!("  - K-mer size");
    println!("  - Scaled value");
    println!("  - .kmerseek.rocksdb extension");

    Ok(())
}
