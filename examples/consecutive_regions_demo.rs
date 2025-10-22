use kmerseek::search::ProteinSearcher;
use kmerseek::{errors::IndexResult, ProteomeIndex};
use tempfile::TempDir;

fn main() -> IndexResult<()> {
    println!("=== Consecutive Regions Demo: BCL2 vs CED9 ===");
    println!("This demo shows how different k-mer sizes affect the number of consecutive regions found.\n");

    // Create temporary directory for test data
    let temp_dir = TempDir::new()?;
    let temp_path = temp_dir.path();

    // Create test FASTA files with BCL2 and CED9 sequences
    let query_fasta = temp_path.join("ced9.fasta");
    let target_fasta = temp_path.join("bcl2.fasta");

    // Create CED9 sequence (query)
    std::fs::write(&query_fasta, ">sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
MTRCTADNSLTNPAYRRRTMATGEMKEFLGIKGTEPTDFGINSDAQDLPSPSRQASTRRM
SIGESIDGKINDWEEPRLDIEGFVVDYFTHRIRQNGMEWFGAPGLPCGVQPEHEMMRVMG
TIFEKKHAENFETFCEQLLAVPRISFSLYQDVVRTVGNAQTDQCPMSYGRLIGLISFGGF
VAAKMMESVELQGQVRNLFVYTSLFIKTRIRNNWKEHNRSWDDFMTLGKQMKEDYERAEA
EKVGRRKQNRRWSMIGAGVTAGAIGIVGVVVCGRMMFSLK")?;

    // Create BCL2 sequence (target)
    std::fs::write(
        &target_fasta,
        ">sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=BCL2 PE=1 SV=2
MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPA
ASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLH
LTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEY
LNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK",
    )?;

    println!("Created test FASTA files:");
    println!("  Query (CED9): {}", query_fasta.display());
    println!("  Target (BCL2): {}", target_fasta.display());

    // Test with k=16 (should find fewer, longer regions)
    println!("\n{}", "=".repeat(60));
    println!("TESTING WITH K=16 (Larger k-mers, fewer consecutive regions)");
    println!("{}", "=".repeat(60));

    test_kmer_size(&query_fasta, &target_fasta, 16, "hp", temp_path)?;

    // Test with k=10 (should find more, shorter regions)
    println!("\n{}", "=".repeat(60));
    println!("TESTING WITH K=10 (Smaller k-mers, more consecutive regions)");
    println!("{}", "=".repeat(60));

    test_kmer_size(&query_fasta, &target_fasta, 10, "hp", temp_path)?;

    println!("\n{}", "=".repeat(60));
    println!("SUMMARY");
    println!("{}", "=".repeat(60));
    println!("• k=16: Finds fewer consecutive regions, but each region is longer");
    println!("• k=10: Finds more consecutive regions, but each region is shorter");
    println!("• Both approaches identify the same biological similarity");
    println!("• The choice of k affects the granularity of the analysis");

    println!("\nDemo completed successfully!");
    Ok(())
}

fn test_kmer_size(
    query_fasta: &std::path::Path,
    target_fasta: &std::path::Path,
    ksize: u32,
    moltype: &str,
    temp_path: &std::path::Path,
) -> IndexResult<()> {
    println!("\n--- Configuration ---");
    println!("K-mer size: {}", ksize);
    println!("Encoding: {}", moltype);
    println!("Scaled: 1");

    // Create target index
    let target_index_path = temp_path.join(format!("target_index_k{}", ksize));
    let target_index = ProteomeIndex::new(
        &target_index_path,
        ksize,
        1, // scaled
        moltype,
        true, // store_raw_sequences
    )?;

    target_index.process_fasta(target_fasta, 1000, 1000)?;
    println!("Target index created with {} signatures", target_index.signature_count());

    // Create searcher
    let searcher = ProteinSearcher::new(target_index);

    // Create query index
    let query_index = ProteomeIndex::new_with_auto_filename(
        query_fasta,
        ksize,
        1, // scaled
        moltype,
        true, // store_raw_sequences
    )?;

    query_index.process_fasta(query_fasta, 1000, 1000)?;

    // Get query signatures
    let query_signatures: Vec<_> =
        query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

    println!("Found {} query signatures", query_signatures.len());

    // Perform basic search
    let results = searcher.search(&query_signatures)?;
    println!("\n--- Basic Search Results ---");
    println!("Found {} matches", results.len());

    if !results.is_empty() {
        let result = &results[0];
        println!("Top match:");
        println!("  Query: {}", result.query_name);
        println!("  Target: {}", result.match_name);
        println!("  Containment: {:.6}", result.containment);
        println!("  Jaccard: {:.6}", result.jaccard);
        println!("  Intersecting k-mers: {}", result.intersect_hashes);
    }

    // Test all consecutive regions
    println!("\n--- All Consecutive Regions ---");
    let detailed_results = searcher.search_with_all_consecutive_regions(&query_signatures)?;
    println!("Found {} consecutive regions", detailed_results.len());

    if !detailed_results.is_empty() {
        println!("\nTop 5 consecutive regions (sorted by length):");
        for (i, result) in detailed_results.iter().take(5).enumerate() {
            println!("  Region {}: Length {} characters", i + 1, result.length);
            println!("    Query:   {} ({}-{})", result.query, result.query_start, result.query_end);
            println!(
                "    Target:  {} ({}-{})",
                result.r#match, result.match_start, result.match_end
            );
            println!("    Encoded: {}", result.encoded);
            println!();
        }

        if detailed_results.len() > 5 {
            println!("  ... and {} more regions", detailed_results.len() - 5);
        }

        // Show statistics
        let total_length: u32 = detailed_results.iter().map(|r| r.length).sum();
        let avg_length = total_length as f64 / detailed_results.len() as f64;
        let max_length = detailed_results.iter().map(|r| r.length).max().unwrap_or(0);
        let min_length = detailed_results.iter().map(|r| r.length).min().unwrap_or(0);

        println!("--- Statistics ---");
        println!("Total regions: {}", detailed_results.len());
        println!("Total length: {} characters", total_length);
        println!("Average length: {:.1} characters", avg_length);
        println!("Max length: {} characters", max_length);
        println!("Min length: {} characters", min_length);
    }

    Ok(())
}
