use tempfile::TempDir;

use kmerseek::search::ProteinSearcher;
use kmerseek::{errors::IndexResult, ProteomeIndex};

/// Test search functionality similar to the Python tests
#[test]
fn test_search_basic() -> IndexResult<()> {
    // Create temporary directory for test data
    let temp_dir = TempDir::new()?;
    let temp_path = temp_dir.path();

    // Create a simple test FASTA file for query (protein sequence)
    let query_fasta = temp_path.join("query.fasta");
    std::fs::write(
        &query_fasta,
        ">test_query\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL",
    )?;

    // Create a simple test FASTA file for target (protein sequence)
    let target_fasta = temp_path.join("target.fasta");
    std::fs::write(
        &target_fasta,
        ">test_target\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL",
    )?;

    // Create target index
    let target_index_path = temp_path.join("target_index");
    let target_index = ProteomeIndex::new(
        &target_index_path,
        10,    // ksize
        1,     // scaled
        "hp",  // moltype
        false, // store_raw_sequences
    )?;

    target_index.process_fasta(&target_fasta, 1000, 1000)?;

    // Create searcher
    let searcher = ProteinSearcher::new(target_index);

    // Create query index
    let query_index = ProteomeIndex::new_with_auto_filename(
        &query_fasta,
        10,    // ksize
        1,     // scaled
        "hp",  // moltype
        false, // store_raw_sequences
    )?;

    query_index.process_fasta(&query_fasta, 1000, 1000)?;

    // Get query signatures
    let query_signatures: Vec<_> =
        query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

    assert!(!query_signatures.is_empty(), "Should have at least one query signature");

    // Perform search using the new comprehensive search method
    let results = searcher.search(&query_signatures)?;

    // Should find at least one match (exact match)
    assert!(!results.is_empty(), "Should find at least one match");

    // Check that the first result has reasonable values
    let first_result = &results[0];
    assert_eq!(first_result.query_name, "test_query");
    assert_eq!(first_result.match_name, "test_target");
    assert!(first_result.containment > 0.0);
    assert!(first_result.jaccard > 0.0);
    assert!(first_result.intersect_hashes > 0);
    // Check that TF-IDF and overlap probability are included
    assert!(first_result.tfidf >= 0.0);
    assert!((0.0..=1.0).contains(&first_result.overlap_probability));

    Ok(())
}

/// Test TF-IDF calculation
#[test]
fn test_tfidf_calculation() -> IndexResult<()> {
    let temp_dir = TempDir::new()?;
    let temp_path = temp_dir.path();

    // Create target FASTA with multiple sequences (protein sequences)
    let target_fasta = temp_path.join("target.fasta");
    std::fs::write(&target_fasta, ">seq1\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL\n>seq2\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL")?;

    // Create target index
    let target_index_path = temp_path.join("target_index");
    let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp", false)?;

    target_index.process_fasta(&target_fasta, 1000, 1000)?;

    // Create searcher
    let searcher = ProteinSearcher::new(target_index);

    // Create query FASTA
    let query_fasta = temp_path.join("query.fasta");
    std::fs::write(
        &query_fasta,
        ">query\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL",
    )?;

    let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;

    query_index.process_fasta(&query_fasta, 1000, 1000)?;

    let query_signatures: Vec<_> =
        query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

    // Calculate TF-IDF
    let tfidf = searcher.calculate_tfidf(&query_signatures[0]);
    assert!(tfidf >= 0.0, "TF-IDF should be non-negative");

    Ok(())
}

/// Test overlap probability calculation
#[test]
fn test_overlap_probability() -> IndexResult<()> {
    let temp_dir = TempDir::new()?;
    let temp_path = temp_dir.path();

    // Create target FASTA
    let target_fasta = temp_path.join("target.fasta");
    std::fs::write(
        &target_fasta,
        ">target\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL",
    )?;

    let target_index_path = temp_path.join("target_index");
    let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp", false)?;

    target_index.process_fasta(&target_fasta, 1000, 1000)?;

    let searcher = ProteinSearcher::new(target_index);

    // Create query FASTA
    let query_fasta = temp_path.join("query.fasta");
    std::fs::write(
        &query_fasta,
        ">query\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL",
    )?;

    let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;

    query_index.process_fasta(&query_fasta, 1000, 1000)?;

    let query_signatures: Vec<_> =
        query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

    let target_signatures: Vec<_> =
        searcher.index().get_signatures().iter().map(|entry| entry.value().clone()).collect();

    // Calculate overlap probability
    let prob = searcher.calculate_overlap_probability(&query_signatures[0], &target_signatures[0]);
    assert!((0.0..=1.0).contains(&prob), "Probability should be between 0 and 1");

    Ok(())
}

/// Test search result structure matches expected format
#[test]
fn test_search_result_structure() -> IndexResult<()> {
    let temp_dir = TempDir::new()?;
    let temp_path = temp_dir.path();

    // Create test data
    let query_fasta = temp_path.join("query.fasta");
    std::fs::write(
        &query_fasta,
        ">test_query\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL",
    )?;

    let target_fasta = temp_path.join("target.fasta");
    std::fs::write(
        &target_fasta,
        ">test_target\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL",
    )?;

    // Create indices
    let target_index_path = temp_path.join("target_index");
    let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp", false)?;
    target_index.process_fasta(&target_fasta, 1000, 1000)?;

    let searcher = ProteinSearcher::new(target_index);

    let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;
    query_index.process_fasta(&query_fasta, 1000, 1000)?;

    let query_signatures: Vec<_> =
        query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

    let results = searcher.search(&query_signatures)?;

    if !results.is_empty() {
        let result = &results[0];

        // Test that all required fields are present and have reasonable values
        assert!(!result.query_name.is_empty());
        assert!(!result.query_md5.is_empty());
        assert!(!result.match_name.is_empty());
        assert!(!result.match_md5.is_empty());
        assert!(!result.moltype.is_empty());

        assert!(result.containment >= 0.0 && result.containment <= 1.0);
        assert!(result.jaccard >= 0.0 && result.jaccard <= 1.0);
        assert!(result.max_containment >= 0.0 && result.max_containment <= 1.0);
        assert!(result.intersect_hashes > 0);
        assert!(result.ksize > 0);
        assert!(result.scaled > 0);

        // Test abundance statistics
        assert!(result.average_abund >= 0.0);
        assert!(result.median_abund >= 0.0);
        assert!(result.std_abund >= 0.0);

        // Test ANI values
        assert!(result.query_containment_ani >= 0.0 && result.query_containment_ani <= 1.0);
        assert!(result.match_containment_ani >= 0.0 && result.match_containment_ani <= 1.0);
        assert!(result.average_containment_ani >= 0.0 && result.average_containment_ani <= 1.0);
        assert!(result.max_containment_ani >= 0.0 && result.max_containment_ani <= 1.0);

        // Test weighted metrics
        assert!(result.n_weighted_found > 0);
        assert!(result.total_weighted_hashes > 0);
        assert!(
            result.containment_target_in_query >= 0.0 && result.containment_target_in_query <= 1.0
        );
        assert!(result.f_weighted_target_in_query >= 0.0);
    }

    Ok(())
}

/// Test that search results are sorted by containment score
#[test]
fn test_search_results_sorted() -> IndexResult<()> {
    let temp_dir = TempDir::new()?;
    let temp_path = temp_dir.path();

    // Create target with multiple sequences of different similarity
    let target_fasta = temp_path.join("target.fasta");
    std::fs::write(&target_fasta, ">exact_match\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL\n>partial_match\nATCGATCGATCGATCA\n>no_match\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL")?;

    let target_index_path = temp_path.join("target_index");
    let target_index = ProteomeIndex::new(&target_index_path, 10, 1, "hp", false)?;
    target_index.process_fasta(&target_fasta, 1000, 1000)?;

    let searcher = ProteinSearcher::new(target_index);

    // Create query
    let query_fasta = temp_path.join("query.fasta");
    std::fs::write(
        &query_fasta,
        ">query\nMKLLILTCLVAVALARPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNEL",
    )?;

    let query_index = ProteomeIndex::new_with_auto_filename(&query_fasta, 10, 1, "hp", false)?;
    query_index.process_fasta(&query_fasta, 1000, 1000)?;

    let query_signatures: Vec<_> =
        query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

    let results = searcher.search(&query_signatures)?;

    // Results should be sorted by containment (descending)
    for i in 1..results.len() {
        assert!(
            results[i - 1].containment >= results[i].containment,
            "Results should be sorted by containment score"
        );
    }

    Ok(())
}
