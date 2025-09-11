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
    std::fs::write(&query_fasta, ">CED9_CAEEL Apoptosis regulator ced-9 OS=Caenorhabditis elegans OX=6239 GN=ced-9 PE=1 SV=1\nMSIGESIDGKINDWEEPGIVGVVVCGRMMFSLKRLDIEGFVVDYFTHRILFVYTSLFIKTRIRNNKVGRRKQNRRWSMIGAGVTARKQNRRWSMIGAGVTA")?;

    // Create target sequences (BCL2 family proteins)
    std::fs::write(&target_fasta, ">BNIP2_HUMAN BCL2/adenovirus E1B 19 kDa protein-interacting protein 2 OS=Homo sapiens OX=9606 GN=BNIP2 PE=1 SV=1\nSIEADILAITGPEDQPLLAVTRPFISSKFSQK\n>ASPP2_HUMAN Apoptosis-stimulating of p53 protein 2 OS=Homo sapiens OX=9606 GN=TP53BP2 PE=1 SV=2\nTIIHREDEDEIEWWWA\n>BAK_HUMAN Bcl-2 homologous antagonist/killer OS=Homo sapiens OX=9606 GN=BAK1 PE=1 SV=1\nHQQEQEAEGVAAPADP\n>BBC3_HUMAN Bcl-2-binding component 3, isoforms 1/2 OS=Homo sapiens OX=9606 GN=BBC3 PE=1 SV=1\nAPAAPTLLPAAYLCAPT\n>FBX10_HUMAN F-box only protein 10 OS=Homo sapiens OX=9606 GN=FBXO10 PE=1 SV=3\nPNWPNQPDVEPESWREAAGIYILYHGNPVVSGN")?;

    println!("Created test FASTA files:");
    println!("  Query: {}", query_fasta.display());
    println!("  Target: {}", target_fasta.display());

    // Create target index
    let target_index_path = temp_path.join("target_index");
    println!("\nCreating target index...");
    let target_index = ProteomeIndex::new(
        &target_index_path,
        16,    // ksize
        5,     // scaled
        "hp",  // moltype
        false, // store_raw_sequences
    )?;

    target_index.process_fasta(&target_fasta, 1000, 1000)?;
    println!("Target index created with {} signatures", target_index.signature_count());

    // Create searcher
    let searcher = ProteinSearcher::new(target_index);

    // Create query index
    println!("\nProcessing query sequences...");
    let query_index = ProteomeIndex::new_with_auto_filename(
        &query_fasta,
        16,    // ksize
        5,     // scaled
        "hp",  // moltype
        false, // store_raw_sequences
    )?;

    query_index.process_fasta(&query_fasta, 1000, 1000)?;

    // Get query signatures
    let query_signatures: Vec<_> =
        query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

    println!("Found {} query signatures", query_signatures.len());

    // Perform search
    println!("\nPerforming search...");
    let results = searcher.search(&query_signatures)?;

    println!("Found {} matches", results.len());

    // Display results
    println!("\n=== Search Results ===");
    for (i, result) in results.iter().enumerate() {
        println!("\nMatch {}:", i + 1);
        println!("  Query: {}", result.query_name);
        println!("  Target: {}", result.match_name);
        println!("  Containment: {:.6}", result.containment);
        println!("  Jaccard: {:.6}", result.jaccard);
        println!("  Intersecting k-mers: {}", result.intersect_hashes);
        println!("  Max containment: {:.6}", result.max_containment);
        println!("  Query containment ANI: {:.6}", result.query_containment_ani);
        println!("  Match containment ANI: {:.6}", result.match_containment_ani);
    }

    // Calculate TF-IDF for the query
    if !query_signatures.is_empty() {
        println!("\n=== TF-IDF Analysis ===");
        let tfidf = searcher.calculate_tfidf(&query_signatures[0]);
        println!("TF-IDF score for query: {:.6}", tfidf);

        // Calculate overlap probabilities with top matches
        println!("\n=== Overlap Probability Analysis ===");
        let target_signatures: Vec<_> =
            searcher.index().get_signatures().iter().map(|entry| entry.value().clone()).collect();

        for result in results.iter().take(3) {
            if let Some(target_sig) =
                target_signatures.iter().find(|sig| sig.signature().name == result.match_name)
            {
                let prob = searcher.calculate_overlap_probability(&query_signatures[0], target_sig);
                println!("Overlap probability with {}: {:.6}", result.match_name, prob);
            }
        }
    }

    // Display search statistics
    println!("\n=== Search Statistics ===");
    let stats = searcher.stats();
    println!("Total signatures in database: {}", stats.total_signatures);
    println!("Unique k-mers in database: {}", stats.kmer_frequencies.len());

    // Test k-mer extraction mode
    println!("\n=== Testing K-mer Extraction Mode ===");
    let detailed_results = searcher.search_with_kmer_extraction(&query_signatures)?;
    println!("Found {} detailed matches", detailed_results.len());

    for result in &detailed_results {
        println!("Detailed match:");
        println!("  Query: {} ({}-{})", result.query, result.query_start, result.query_end);
        println!("  Target: {} ({}-{})", result.r#match, result.match_start, result.match_end);
        println!("  Encoded: {}", result.encoded);
        println!("  Length: {}", result.length);
    }

    println!("\nDemo completed successfully!");
    Ok(())
}
