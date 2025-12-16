use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::tempdir;

use crate::search::SearchResultWithRegionCsv;
use crate::tests::test_fixtures::{TEST_CED9_FASTA, TEST_FASTA_GZ};

#[test]
fn test_cli_help() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("kmerseek-rust")?;
    cmd.arg("--help");
    cmd.assert().success().stdout(predicate::str::contains(
        "Efficient protein domain annotation search with reduced amino acid k-mers",
    ));

    Ok(())
}

#[test]
fn test_cli_index_help() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("kmerseek-rust")?;
    cmd.args(["index", "--help"]);
    cmd.assert().success().stdout(predicate::str::contains("Index a FASTA file"));

    Ok(())
}

#[test]
fn test_cli_index_basic() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let output_path = temp_dir.path().join("test_output.db");

    let mut cmd = Command::cargo_bin("kmerseek-rust")?;
    cmd.args([
        "index",
        "--input",
        "tests/testdata/fasta/ced9.fasta",
        "--output",
        output_path.to_str().unwrap(),
        "--ksize",
        "5",
        "--encoding",
        "protein",
    ]);

    cmd.assert().success().stdout(predicate::str::contains("Indexing completed successfully!"));

    // Check that the output database was created
    assert!(output_path.exists());

    Ok(())
}

#[test]
fn test_cli_index_gzipped() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let output_path = temp_dir.path().join("test_output_gz.db");

    let mut cmd = Command::cargo_bin("kmerseek-rust")?;
    cmd.args([
        "index",
        "--input",
        TEST_FASTA_GZ,
        "--output",
        output_path.to_str().unwrap(),
        "--ksize",
        "10",
        "--encoding",
        "hp",
    ]);

    cmd.assert().success().stdout(predicate::str::contains("Indexing completed successfully!"));

    // Check that the output database was created
    assert!(output_path.exists());

    Ok(())
}

#[test]
fn test_cli_index_different_encodings() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;

    for encoding in ["protein", "dayhoff", "hp"] {
        let output_path = temp_dir.path().join(format!("test_output_{}.db", encoding));

        let mut cmd = Command::cargo_bin("kmerseek-rust")?;
        cmd.args([
            "index",
            "--input",
            TEST_CED9_FASTA,
            "--output",
            output_path.to_str().unwrap(),
            "--ksize",
            "8",
            "--encoding",
            encoding,
        ]);

        cmd.assert().success().stdout(predicate::str::contains("Indexing completed successfully!"));

        // Check that the output database was created
        assert!(output_path.exists());
    }

    Ok(())
}

#[test]
fn test_cli_index_nonexistent_file() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let output_path = temp_dir.path().join("test_error.db");

    let mut cmd = Command::cargo_bin("kmerseek-rust")?;
    cmd.args([
        "index",
        "--input",
        "nonexistent.fasta",
        "--output",
        output_path.to_str().unwrap(),
        "--ksize",
        "5",
    ]);

    cmd.assert().failure().stderr(predicate::str::contains("No such file or directory"));

    Ok(())
}

#[test]
fn test_cli_index_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("kmerseek-rust")?;
    cmd.args(["index", "--ksize", "5"]);

    cmd.assert().failure().stderr(predicate::str::contains("required"));

    Ok(())
}

#[test]
fn test_cli_search_bcl2_ced9() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;

    // Step 1: Create target index from bcl2_first25 FASTA
    let target_index_path = temp_dir.path().join("target_index.db");
    let mut index_cmd = Command::cargo_bin("kmerseek-rust")?;
    index_cmd.args([
        "index",
        "--input",
        TEST_FASTA_GZ,
        "--output",
        target_index_path.to_str().unwrap(),
        "--ksize",
        "12",
        "--scaled",
        "1",
        "--encoding",
        "hp",
    ]);

    index_cmd
        .assert()
        .success()
        .stdout(predicate::str::contains("Indexing completed successfully!"));
    assert!(target_index_path.exists(), "Target index should be created");

    // Step 2: Run search with CED9 as query
    let output_csv = temp_dir.path().join("search_results.csv");
    let mut search_cmd = Command::cargo_bin("kmerseek-rust")?;
    search_cmd.args([
        "search",
        "--query",
        TEST_CED9_FASTA,
        "--target",
        target_index_path.to_str().unwrap(),
        "--output",
        output_csv.to_str().unwrap(),
        "--ksize",
        "12",
        "--scaled",
        "1",
        "--encoding",
        "hp",
    ]);

    search_cmd
        .assert()
        .success()
        .stdout(predicate::str::contains("Found"))
        .stdout(predicate::str::contains("matches"));

    // Step 3: Verify output CSV exists and contains expected results
    assert!(output_csv.exists(), "Search results CSV should be created");

    // Verify CSV file is not empty
    let csv_content = std::fs::read_to_string(&output_csv)?;
    assert!(!csv_content.is_empty(), "CSV file should not be empty");
    assert!(csv_content.lines().count() > 1, "CSV should have at least a header and one data row");

    // Read and verify CSV contents
    // WHY: We deserialize into SearchResultWithRegionCsv which is the same struct used for CSV output.
    // This ensures type safety and matches exactly what the CLI writes to the CSV file.
    // Each matched region gets its own row, so we may have multiple rows per SearchResult.
    let mut reader = csv::Reader::from_path(&output_csv)?;
    let mut results_count = 0;
    let mut found_bcl2 = false;

    for result in reader.deserialize::<SearchResultWithRegionCsv>() {
        let record = match result {
            Ok(r) => r,
            Err(e) => {
                eprintln!("CSV deserialization error: {}", e);
                eprintln!(
                    "CSV file content (first 500 chars): {}",
                    &csv_content[..csv_content.len().min(500)]
                );
                return Err(Box::new(e));
            }
        };
        results_count += 1;

        // Check if we found BCL2_HUMAN in the results
        if record.target_name.contains("BCL2_HUMAN") {
            found_bcl2 = true;

            // Verify query name contains CED9
            assert!(
                record.query_name.contains("CED9_CAEEL"),
                "Query name should contain CED9_CAEEL, got {}",
                record.query_name
            );

            // Verify ksize, scaled, and moltype match
            assert_eq!(record.ksize, 12, "Ksize should be 12");
            assert_eq!(record.scaled, 1, "Scaled should be 1");
            assert_eq!(record.moltype, "hp", "Moltype should be hp");

            // Verify we have intersecting k-mers
            assert!(
                record.n_intersecting_hashes > 0,
                "Should have intersecting k-mers between CED9 and BCL2, got {}",
                record.n_intersecting_hashes
            );

            // Verify similarity metrics are valid
            assert!(
                record.containment > 0.0 && record.containment <= 1.0,
                "Containment should be in [0, 1], got {}",
                record.containment
            );
            assert!(
                record.jaccard > 0.0 && record.jaccard <= 1.0,
                "Jaccard should be in [0, 1], got {}",
                record.jaccard
            );

            // Verify TF-IDF is meaningful (should not be 0 with multiple signatures)
            assert!(record.tfidf >= 0.0, "TF-IDF should be non-negative, got {}", record.tfidf);
        }
    }

    // Verify we found at least 25 matches (should be 25 SearchResults, but each may have multiple
    // matched regions, so we'll have more rows than 25)
    assert!(
        results_count >= 25,
        "Should find at least 25 rows (one per SearchResult, possibly more if matched regions exist), got {}",
        results_count
    );

    // Verify we found BCL2_HUMAN in the results
    assert!(found_bcl2, "Should find a match with BCL2_HUMAN in the search results");

    Ok(())
}
