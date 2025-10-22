use kmerseek::search::ProteinSearcher;
use kmerseek::{errors::IndexResult, ProteomeIndex};
use tempfile::TempDir;

fn main() -> IndexResult<()> {
    println!("=== Protein Signature Search Demo ===");

    // Create temporary directory for test data
    let temp_dir = TempDir::new()?;
    let temp_path = temp_dir.path();

    // Create test FASTA files
    let query_fasta = temp_path.join("query.fasta");
    let target_fasta = temp_path.join("target.fasta");

    // Create query sequence (CED9-like sequence)
    std::fs::write(&query_fasta, ">sp|P41958|CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1
MTRCTADNSLTNPAYRRRTMATGEMKEFLGIKGTEPTDFGINSDAQDLPSPSRQASTRRM
SIGESIDGKINDWEEPRLDIEGFVVDYFTHRIRQNGMEWFGAPGLPCGVQPEHEMMRVMG
TIFEKKHAENFETFCEQLLAVPRISFSLYQDVVRTVGNAQTDQCPMSYGRLIGLISFGGF
VAAKMMESVELQGQVRNLFVYTSLFIKTRIRNNWKEHNRSWDDFMTLGKQMKEDYERAEA
EKVGRRKQNRRWSMIGAGVTAGAIGIVGVVVCGRMMFSLK")?;

    // Create target sequences (BCL2 family proteins)
    std::fs::write(
        &target_fasta,
        ">sp|P10415|BCL2_HUMAN Apoptosis regulator Bcl-2 OS=Homo sapiens OX=9606 GN=BCL2 PE=1 SV=2
MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPA
ASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLH
LTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEY
LNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK",
    )?;

    println!("Created test FASTA files:");
    println!("  Query: {}", query_fasta.display());
    println!("  Target: {}", target_fasta.display());

    // Create target index
    let target_index_path = temp_path.join("target_index");
    println!("\nCreating target index...");
    let target_index = ProteomeIndex::new(
        &target_index_path,
        10,   // ksize
        1,    // scaled
        "hp", // moltype
        true, // store_raw_sequences
    )?;

    target_index.process_fasta(&target_fasta, 1000, 1000)?;
    println!("Target index created with {} signatures", target_index.signature_count());

    // Create searcher
    let searcher = ProteinSearcher::new(target_index);

    // Create query index
    println!("\nProcessing query sequences...");
    let query_index = ProteomeIndex::new_with_auto_filename(
        &query_fasta,
        10,   // ksize
        1,    // scaled
        "hp", // moltype
        true, // store_raw_sequences
    )?;

    query_index.process_fasta(&query_fasta, 1000, 1000)?;

    // Get query signatures
    let query_signatures: Vec<_> =
        query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

    println!("Found {} query signatures", query_signatures.len());

    // Perform basic search first
    println!("\nPerforming basic search...");
    let basic_results = searcher.search(&query_signatures)?;
    println!("Found {} basic matches", basic_results.len());

    // Display basic results
    println!("\n=== Basic Search Results ===");
    for (i, result) in basic_results.iter().enumerate() {
        println!("\nMatch {}:", i + 1);
        println!("  Query: {}", result.query_name);
        println!("  Target: {}", result.match_name);
        println!("  Containment: {:.6}", result.containment);
        println!("  Jaccard: {:.6}", result.jaccard);
        println!("  Intersecting k-mers: {}", result.intersect_hashes);
        println!("  TF-IDF: {:.6}", result.tfidf);
        println!("  Overlap probability: {:.6}", result.overlap_probability);
    }

    // Perform detailed search with all consecutive regions
    println!("\nPerforming detailed search (all consecutive regions)...");
    let detailed_results = searcher.search_with_all_consecutive_regions(&query_signatures)?;
    println!("Found {} detailed matches across all consecutive regions", detailed_results.len());

    // Display detailed results
    println!("\n=== Detailed Search Results (All Consecutive Regions) ===");
    for (i, result) in detailed_results.iter().enumerate() {
        println!("\nMatch {}:", i + 1);
        println!("  Query:   {} ({}-{})", result.query, result.query_start, result.query_end);
        println!("  Encoded: {}", result.encoded);
        println!("  Target:  {} ({}-{})", result.r#match, result.match_start, result.match_end);
        println!("  Length: {}", result.length);
        println!("  Match: {}", result.match_name);
        println!();
    }

    // Display search statistics
    println!("\n=== Search Statistics ===");
    let stats = searcher.stats();
    println!("Total signatures in database: {}", stats.total_signatures);
    println!("Unique k-mers in database: {}", stats.kmer_frequencies.len());

    println!("\nDemo completed successfully!");
    Ok(())
}
