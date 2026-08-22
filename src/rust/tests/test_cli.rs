use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::tempdir;

use approx::assert_relative_eq;

use crate::hp_alphabets::HpAlphabet;
use crate::reduced_alphabets::ReducedAlphabet;
use crate::search::SearchResultCsv;
use crate::tests::test_fixtures::{TEST_CED9_FASTA, TEST_FASTA_GZ};

#[test]
fn test_cli_help() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.arg("--help");
    cmd.assert().success().stdout(predicate::str::contains(
        "Efficient protein domain annotation search with reduced amino acid k-mers",
    ));

    Ok(())
}

#[test]
fn test_cli_index_help() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.args(["index", "--help"]);
    cmd.assert().success().stdout(predicate::str::contains("Index a FASTA file"));

    Ok(())
}

#[test]
fn test_cli_index_basic() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let output_path = temp_dir.path().join("test_output.db");

    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.args([
        "index",
        "--input",
        "tests/testdata/fasta/ced9.fasta",
        "--output",
        output_path.to_str().unwrap(),
        "--ksize",
        "5",
        "--alphabet",
        "protein20",
    ]);

    cmd.assert().success().stderr(predicate::str::contains("Indexing completed successfully!"));

    // Check that the output database was created
    assert!(output_path.exists());

    Ok(())
}

#[test]
fn test_cli_index_gzipped() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let output_path = temp_dir.path().join("test_output_gz.db");

    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.args([
        "index",
        "--input",
        TEST_FASTA_GZ,
        "--output",
        output_path.to_str().unwrap(),
        "--ksize",
        "10",
        "--alphabet",
        "hp_lehninger2",
    ]);

    cmd.assert().success().stderr(predicate::str::contains("Indexing completed successfully!"));

    // Check that the output database was created
    assert!(output_path.exists());

    Ok(())
}

#[test]
fn test_cli_index_every_alphabet() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;

    // Every alphabet, not a sample: a wrong table or hash function shows up as an indexing
    // failure. Built from the alphabet lists rather than hardcoded, so a newly added
    // alphabet is covered without editing this test.
    let alphabets: Vec<String> = ["protein20".to_string(), "dayhoff6".to_string()]
        .into_iter()
        .chain(HpAlphabet::all_named().iter().map(HpAlphabet::to_moltype))
        .chain(ReducedAlphabet::all().iter().map(ReducedAlphabet::to_moltype))
        .collect();
    assert_eq!(alphabets.len(), 18, "every alphabet must be exercised here");

    for alphabet in &alphabets {
        let output_path = temp_dir.path().join(format!("test_output_{}.db", alphabet));

        let mut cmd = Command::cargo_bin("kmerseek")?;
        cmd.args([
            "index",
            "--input",
            TEST_CED9_FASTA,
            "--output",
            output_path.to_str().unwrap(),
            "--ksize",
            "8",
            "--alphabet",
            alphabet,
        ]);

        cmd.assert().success().stderr(predicate::str::contains("Indexing completed successfully!"));

        // Check that the output database was created
        assert!(output_path.exists());
    }

    Ok(())
}

/// Low-complexity removal must round-trip: `index --remove-low-complexity`
/// persists the setting, and `search` picks it up from the index without the
/// user restating it. A mismatch would silently skew containment, so this
/// asserts on the reported value rather than just on exit status.
#[test]
fn test_cli_remove_low_complexity_round_trips_to_search() -> Result<(), Box<dyn std::error::Error>>
{
    let temp_dir = tempdir()?;

    for (flag, expected) in [(true, "true"), (false, "false")] {
        let index_path = temp_dir.path().join(format!("lc_{}.db", flag));

        let mut index_cmd = Command::cargo_bin("kmerseek")?;
        index_cmd.args([
            "index",
            "--input",
            TEST_FASTA_GZ,
            "--output",
            index_path.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp_lehninger2",
        ]);
        if flag {
            index_cmd.arg("--remove-low-complexity");
        }
        index_cmd
            .assert()
            .success()
            .stderr(predicate::str::contains("Indexing completed successfully!"));

        // Search must recover the setting from the index, with no flag passed here.
        let mut search_cmd = Command::cargo_bin("kmerseek")?;
        search_cmd.args([
            "search",
            "--query",
            TEST_CED9_FASTA,
            "--target",
            index_path.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp_lehninger2",
        ]);
        let state = if expected == "true" { "REMOVED" } else { "KEPT" };
        search_cmd
            .assert()
            .success()
            .stderr(predicate::str::contains(format!(
                "Index: low-complexity k-mers were {state} when it was built"
            )))
            .stderr(predicate::str::contains(format!(
                "This search: low-complexity k-mers are {state} from query sketches"
            )))
            .stderr(predicate::str::contains("(matching the index)"))
            // Agreeing with the index must not warn.
            .stderr(predicate::str::contains("WARNING").not());
    }

    Ok(())
}

/// Overriding the index's setting at search time is allowed but warned about,
/// because containment is only comparable when both sides drop the same k-mers.
#[test]
fn test_cli_search_remove_low_complexity_override_warns() -> Result<(), Box<dyn std::error::Error>>
{
    let temp_dir = tempdir()?;
    let index_path = temp_dir.path().join("kept.db");

    let mut index_cmd = Command::cargo_bin("kmerseek")?;
    index_cmd.args([
        "index",
        "--input",
        TEST_FASTA_GZ,
        "--output",
        index_path.to_str().unwrap(),
        "--ksize",
        "12",
        "--alphabet",
        "hp_lehninger2",
    ]);
    index_cmd.assert().success();

    let mut search_cmd = Command::cargo_bin("kmerseek")?;
    search_cmd.args([
        "search",
        "--query",
        TEST_CED9_FASTA,
        "--target",
        index_path.to_str().unwrap(),
        "--ksize",
        "12",
        "--alphabet",
        "hp_lehninger2",
        "--remove-low-complexity",
    ]);

    search_cmd
        .assert()
        .success()
        .stderr(predicate::str::contains("Index: low-complexity k-mers were KEPT"))
        .stderr(predicate::str::contains(
            "This search: low-complexity k-mers are REMOVED from query sketches",
        ))
        .stderr(predicate::str::contains("WARNING: this disagrees with the index"));

    Ok(())
}

/// The results file records the setting per row, so a CSV is self-describing
/// without needing the command line that produced it.
#[test]
fn test_cli_search_csv_records_remove_low_complexity() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;

    for (flag, expected) in [(true, "true"), (false, "false")] {
        let index_path = temp_dir.path().join(format!("csv_{flag}.db"));
        let csv_path = temp_dir.path().join(format!("csv_{flag}.csv"));

        let mut index_cmd = Command::cargo_bin("kmerseek")?;
        index_cmd.args([
            "index",
            "--input",
            TEST_FASTA_GZ,
            "--output",
            index_path.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp_lehninger2",
        ]);
        if flag {
            index_cmd.arg("--remove-low-complexity");
        }
        index_cmd.assert().success();

        let mut search_cmd = Command::cargo_bin("kmerseek")?;
        search_cmd.args([
            "search",
            "--query",
            TEST_CED9_FASTA,
            "--target",
            index_path.to_str().unwrap(),
            "--output",
            csv_path.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp_lehninger2",
        ]);
        search_cmd.assert().success();

        let csv = std::fs::read_to_string(&csv_path)?;
        let mut lines = csv.lines();
        let header: Vec<&str> = lines.next().expect("header row").split(',').collect();
        let col = header
            .iter()
            .position(|h| *h == "remove_low_complexity")
            .expect("remove_low_complexity column");
        // It sits with the other run metadata rather than at the end.
        assert_eq!(header[col - 1], "moltype");

        let first = lines.next().expect("at least one result row");
        assert_eq!(first.split(',').nth(col), Some(expected));
    }

    Ok(())
}

/// Auto-generated names must differ with and without removal, or indexing the
/// same FASTA both ways would silently clobber the first database.
#[test]
fn test_cli_remove_low_complexity_auto_filename_does_not_collide(
) -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let fasta = temp_dir.path().join("ced9.fasta");
    std::fs::copy(TEST_CED9_FASTA, &fasta)?;

    for flag in [false, true] {
        let mut cmd = Command::cargo_bin("kmerseek")?;
        cmd.args([
            "index",
            "--input",
            fasta.to_str().unwrap(),
            "--ksize",
            "8",
            "--alphabet",
            "hp_lehninger2",
        ]);
        if flag {
            cmd.arg("--remove-low-complexity");
        }
        cmd.assert().success();
    }

    // `--encoding hp` is stored under the alphabet's current name, so that is what the
    // generated filename carries.
    let kept_all = temp_dir.path().join("ced9.fasta.hp_lehninger2.k8.scaled1.kmerseek.rocksdb");
    let removed = temp_dir
        .path()
        .join("ced9.fasta.hp_lehninger2.k8.scaled1.nolowcomplexity.kmerseek.rocksdb");
    assert!(kept_all.exists(), "index keeping every k-mer should keep its historical name");
    assert!(removed.exists(), "index with removal should get its own name");

    Ok(())
}

#[test]
fn test_cli_index_nonexistent_file() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let output_path = temp_dir.path().join("test_error.db");

    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.args([
        "index",
        "--input",
        "nonexistent.fasta",
        "--output",
        output_path.to_str().unwrap(),
        "--ksize",
        "5",
    ]);

    // WHY: The index command now validates file existence and returns a clear, custom
    // error message instead of the lower-level OS error. We assert on the stable,
    // human-friendly prefix ("FASTA file not found") rather than the exact OS string
    // ("No such file or directory"), which can vary across platforms.
    cmd.assert().failure().stderr(predicate::str::contains("FASTA file not found"));

    Ok(())
}

#[test]
fn test_cli_index_missing_required_args() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.args(["index", "--ksize", "5"]);

    cmd.assert().failure().stderr(predicate::str::contains("required"));

    Ok(())
}

#[test]
fn test_cli_search_bcl2_ced9() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;

    // Step 1: Create target index from bcl2_first25 FASTA
    let target_index_path = temp_dir.path().join("target_index.db");
    let mut index_cmd = Command::cargo_bin("kmerseek")?;
    index_cmd.args([
        "index",
        "--input",
        TEST_FASTA_GZ,
        "--output",
        target_index_path.to_str().unwrap(),
        "--ksize",
        "12",
        "--alphabet",
        "hp_lehninger2",
    ]);

    index_cmd
        .assert()
        .success()
        .stderr(predicate::str::contains("Indexing completed successfully!"));
    assert!(target_index_path.exists(), "Target index should be created");

    // Step 2: Run search with CED9 as query
    let output_csv = temp_dir.path().join("search_results.csv");
    let mut search_cmd = Command::cargo_bin("kmerseek")?;
    // WHY: This test validates known ground-truth values for the BCL2/CED9 match (see
    // n_intersecting_hashes/containment/etc. assertions below). The search command's default
    // p-value cutoff (0.05) doesn't work here: this fixture database is only 25 sequences, all
    // BCL2-family, so most k-mers are common across the DB and expected_shared_kmers is
    // inflated - the real poisson_pvalue for BCL2_HUMAN vs CED9 is 0.636 (n_intersecting_hashes
    // 24 vs expected_shared_kmers 25.4, i.e. at-or-below chance in this curated set). p < 0.7
    // is loose enough to keep that match on this specific database without disabling the
    // p-value filter outright.
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
        "--alphabet",
        "hp_lehninger2",
        "--min-shared-kmers",
        "0",
        "--max-pvalue",
        "0.7",
    ]);

    search_cmd.assert().success().stderr(predicate::str::contains("Total matches"));

    // Step 3: Verify output CSV exists and contains expected results
    assert!(output_csv.exists(), "Search results CSV should be created");

    // Verify CSV file is not empty
    let csv_content = std::fs::read_to_string(&output_csv)?;
    assert!(!csv_content.is_empty(), "CSV file should not be empty");
    // WHY: 364 (fully unfiltered) drops to 243 once --max-pvalue 0.7 excludes matches that
    // aren't enriched above chance in this small, BCL2-heavy fixture database (see comment above).
    assert!(
        csv_content.lines().count() == 243,
        "CSV should have 243 rows, found {} rows",
        csv_content.lines().count()
    );

    // Read and verify CSV contents
    // WHY: We deserialize into SearchResultCsv which is the same struct used for CSV output.
    // This ensures type safety and matches exactly what the CLI writes to the CSV file.
    // Each matched region gets its own row, so we may have multiple rows per SearchResult.
    let mut reader = csv::Reader::from_path(&output_csv)?;
    let mut results_count = 0;
    let mut found_bcl2 = false;

    for result in reader.deserialize::<SearchResultCsv>() {
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
            assert_eq!(
                record.moltype, "hp_lehninger2",
                "Moltype should be the normalized name for the hp alphabet"
            );

            // Verify we have intersecting k-mers
            // WHY: These values come from the compare test in search.rs (BCL2_CED9_K12 constant),
            // which uses sourmash sig overlap to get ground truth values for k=12, scaled=1, hp encoding.
            assert_eq!(
                record.n_intersecting_hashes, 24,
                "Should have 24 intersecting k-mers between CED9 and BCL2, got {}",
                record.n_intersecting_hashes
            );

            // Verify similarity metrics match expected values from compare test
            assert_relative_eq!(record.containment, 0.09091, epsilon = 1e-5);
            assert_relative_eq!(record.jaccard, 0.05217, epsilon = 1e-5);
            assert_relative_eq!(record.max_containment, 0.10909, epsilon = 1e-5);
            assert_relative_eq!(record.containment_target_in_query, 0.10909, epsilon = 1e-5);

            // Verify TF-IDF is meaningful (should not be 0 with multiple signatures)
            assert_relative_eq!(record.query_tfidf, 565.119680433367, epsilon = 1e-5);
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
