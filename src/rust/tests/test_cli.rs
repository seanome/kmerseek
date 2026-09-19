use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::tempdir;

use approx::assert_relative_eq;
use rstest::rstest;

use crate::alphabets::Alphabet;
use crate::search::SearchResultCsv;
use crate::tests::test_fixtures::{TEST_BLC2_FASTA, TEST_CED9_FASTA, TEST_FASTA_GZ};

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
    let alphabets: Vec<&str> = Alphabet::all().iter().map(Alphabet::to_moltype).collect();
    assert_eq!(alphabets.len(), 19, "every alphabet must be exercised here");

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

/// `--ksize 0` used to abort with an integer-underflow panic partway through
/// indexing. It should fail cleanly, with a message naming the problem.
#[test]
fn test_cli_index_rejects_zero_ksize() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let output_path = temp_dir.path().join("zero_ksize.db");

    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.args([
        "index",
        "--input",
        TEST_CED9_FASTA,
        "--output",
        output_path.to_str().unwrap(),
        "--ksize",
        "0",
        "--alphabet",
        "hp_lehninger2",
    ]);

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("K-mer size must be greater than 0"))
        .stderr(predicate::str::contains("panicked").not());
    assert!(!output_path.exists(), "a rejected k-mer size should leave no database behind");

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

        // Parsed with a real CSV reader, not split(','): some FASTA descriptions
        // hold a comma, so the writer quotes those fields and a naive split
        // shifts every column after them.
        let mut reader = csv::Reader::from_path(&csv_path)?;
        let header = reader.headers()?.clone();
        let col = header
            .iter()
            .position(|h| h == "remove_low_complexity")
            .expect("remove_low_complexity column");
        // It sits with the other run metadata rather than at the end.
        assert_eq!(&header[col - 1], "moltype");

        let mut rows = 0;
        for record in reader.records() {
            assert_eq!(&record?[col], expected);
            rows += 1;
        }
        // 330 regions: chaining seeds per diagonal merges the pieces a repeated k-mer used to
        // split (362 before).
        assert_eq!(rows, 330, "ced9 against the 25-sequence bcl2 index at k=12");
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

    // `--alphabet hp_lehninger2` is stored under that name, so that is what the
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
    // WHY: --max-pvalue 0.7 excludes matches that aren't enriched above chance in this small,
    // BCL2-heavy fixture database (see comment above). 218 regions plus the header: chaining
    // seeds per diagonal merges the pieces a repeated k-mer used to split, which is why this
    // was 243 before.
    assert!(
        csv_content.lines().count() == 219,
        "CSV should have 219 rows, found {} rows",
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

/// Index the 25-sequence test FASTA at `ksize`/`scaled` in `hp_lehninger2`, search CED9
/// against it with every result filter open, and return the rows hitting BCL2_HUMAN sorted
/// by query start. `scaled` is read back from the database, so `search` takes no flag.
fn index_and_search_bcl2(
    ksize: u32,
    scaled: u32,
) -> Result<Vec<SearchResultCsv>, Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let target_index_path = temp_dir.path().join("target.db");
    Command::cargo_bin("kmerseek")?
        .args([
            "index",
            "--input",
            TEST_FASTA_GZ,
            "--output",
            target_index_path.to_str().unwrap(),
            "--ksize",
            &ksize.to_string(),
            "--scaled",
            &scaled.to_string(),
            "--alphabet",
            "hp_lehninger2",
        ])
        .assert()
        .success()
        .stderr(predicate::str::contains(format!("Scaled: {scaled}")));

    let output_csv = temp_dir.path().join("hits.csv");
    Command::cargo_bin("kmerseek")?
        .args([
            "search",
            "--query",
            TEST_CED9_FASTA,
            "--target",
            target_index_path.to_str().unwrap(),
            "--output",
            output_csv.to_str().unwrap(),
            "--ksize",
            &ksize.to_string(),
            "--alphabet",
            "hp_lehninger2",
            "--min-shared-kmers",
            "0",
            "--max-query-pvalue",
            "1.0",
        ])
        .assert()
        .success()
        .stderr(predicate::str::contains(format!("Scaled: {scaled} (detected: {scaled})")));

    let mut rows = Vec::new();
    for record in csv::Reader::from_path(&output_csv)?.deserialize::<SearchResultCsv>() {
        let record = record?;
        assert_eq!(record.scaled, scaled);
        if record.target_name.contains("BCL2_HUMAN") {
            rows.push(record);
        }
    }
    rows.sort_by_key(|r| (r.region_start, r.target_start));
    Ok(rows)
}

/// `--scaled 5` keeps about a fifth of the k-mers. The same CED9 vs BCL2 search then reports
/// exactly the four regions the unit tests fix for this pair at scaled=5 (see
/// `test_sampled_regions_scaled_5_exact`), each spanning the whole exact match, with
/// `region_n_shared_kmers` counting only the kept k-mers.
#[test]
fn test_cli_index_scaled_5_search_regions() -> Result<(), Box<dyn std::error::Error>> {
    let rows = index_and_search_bcl2(12, 5)?;
    let regions: Vec<_> = rows
        .iter()
        .map(|r| {
            (
                r.region_start,
                r.region_end,
                r.target_start,
                r.target_end,
                r.region_n_shared_kmers,
                r.region_subseq.as_str(),
            )
        })
        .collect();
    assert_eq!(
        regions,
        [
            (145, 159, 130, 144, 1, "FSLYQDVVRTVGNA"),
            (162, 181, 138, 157, 1, "QCPMSYGRLIGLISFGGFV"),
            (253, 266, 80, 93, 1, "MIGAGVTAGAIGI"),
            (267, 280, 200, 213, 2, "GVVVCGRMMFSLK"),
        ]
    );
    for r in &rows {
        assert_eq!(r.n_intersecting_hashes, 5);
    }
    Ok(())
}

/// The BH1 match between CED9 and BCL2, CED9 162..181 against BCL2 138..157:
///
/// ```text
/// Ced9 pr: …RTVGNAQTD QCPMSYGRLIGLISFGGFV AAKMMESVE…
/// Ced9 hp: …pphhphppp pphhphhphhhhhphhhhh hhphhpphp…
///                     |||||||||||||||||||
/// BCL2 hp: …hhphhpphh pphhphhphhhhhphhhhh phpphppph…
/// BCL2 pr: …FATVVEELF RDGVNWGRIVAFFEFGGVM CVESVNREM…
/// ```
///
/// Whenever any of its k-mers survives the cutoff the CLI reports the whole 19-residue
/// stretch, both residue strings and the encoding intact, and `region_n_shared_kmers` counts
/// the survivors: 8, 4 and 1 of its 8 k-mers at k=12, and 5, 5, 3 and 1 of its 5 at k=15.
/// At k=12 and scaled=10 all 8 hash above the cutoff and the match is missed, which is what
/// the README's survival formula predicts for a match this short at that sampling rate.
#[rstest]
#[case::k12_scaled_1(12, 1, Some(8))]
#[case::k12_scaled_2(12, 2, Some(4))]
#[case::k12_scaled_5(12, 5, Some(1))]
#[case::k12_scaled_10(12, 10, None)]
#[case::k15_scaled_1(15, 1, Some(5))]
#[case::k15_scaled_2(15, 2, Some(5))]
#[case::k15_scaled_5(15, 5, Some(3))]
#[case::k15_scaled_10(15, 10, Some(1))]
fn test_cli_landmark_region_across_scaled(
    #[case] ksize: u32,
    #[case] scaled: u32,
    #[case] n_shared: Option<u32>,
) -> Result<(), Box<dyn std::error::Error>> {
    let rows = index_and_search_bcl2(ksize, scaled)?;
    let landmark: Vec<_> =
        rows.iter().filter(|r| r.region_start == 162 && r.target_start == 138).collect();
    let Some(n_shared) = n_shared else {
        assert!(landmark.is_empty(), "k={ksize} scaled={scaled}: {landmark:?}");
        return Ok(());
    };
    assert_eq!(landmark.len(), 1, "k={ksize} scaled={scaled}");
    let r = landmark[0];
    assert_eq!(r.region_end, 181);
    assert_eq!(r.target_end, 157);
    assert_eq!(r.region_subseq, "QCPMSYGRLIGLISFGGFV");
    assert_eq!(r.target_subseq, "RDGVNWGRIVAFFEFGGVM");
    assert_eq!(r.moltype_seq, "pphhphhphhhhhphhhhh");
    assert_eq!(r.region_length, 19);
    assert_eq!(r.region_n_shared_kmers, n_shared);
    Ok(())
}

/// Two indexes are compared sketch for sketch, so `--query-is-index` refuses a query index
/// whose scaled differs from the target's instead of reaching the assertion in
/// `find_matched_regions`.
#[test]
fn test_cli_query_is_index_rejects_scaled_mismatch() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let mut dbs = Vec::new();
    for scaled in ["1", "5"] {
        let db = temp_dir.path().join(format!("scaled_{scaled}.db"));
        Command::cargo_bin("kmerseek")?
            .args([
                "index",
                "--input",
                TEST_FASTA_GZ,
                "--output",
                db.to_str().unwrap(),
                "--ksize",
                "12",
                "--scaled",
                scaled,
                "--alphabet",
                "hp_lehninger2",
            ])
            .assert()
            .success();
        dbs.push(db);
    }
    Command::cargo_bin("kmerseek")?
        .args([
            "search",
            "--query",
            dbs[1].to_str().unwrap(),
            "--query-is-index",
            "--target",
            dbs[0].to_str().unwrap(),
            "--output",
            temp_dir.path().join("hits.csv").to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp_lehninger2",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "Query index (ksize=12, scaled=5, alphabet=hp_lehninger2) was not built with the \
             target's parameters (ksize=12, scaled=1, alphabet=hp_lehninger2)",
        ));
    Ok(())
}

/// A bad `--scaled` is rejected before any database is created.
#[test]
fn test_cli_index_rejects_bad_scaled() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    for (value, message) in [("0", "must be greater than 0"), ("11", "Scaled value too large")] {
        let db = temp_dir.path().join(format!("scaled_{value}.db"));
        Command::cargo_bin("kmerseek")?
            .args([
                "index",
                "--input",
                TEST_FASTA_GZ,
                "--output",
                db.to_str().unwrap(),
                "--scaled",
                value,
            ])
            .assert()
            .failure()
            .stderr(predicate::str::contains(message));
        assert!(!db.exists(), "--scaled {value} must not leave a database behind");
    }
    Ok(())
}

#[test]
fn test_cli_pair_bcl2_ced9_writes_json() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let output = temp_dir.path().join("pair.json");
    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.args([
        "pair",
        "--query",
        TEST_BLC2_FASTA,
        "--target",
        TEST_CED9_FASTA,
        "--ksize",
        "12",
        "--alphabet",
        "hp",
        "--output",
        output.to_str().unwrap(),
    ]);
    cmd.assert().success().stderr(predicate::str::contains(
        "27 shared 12-mers in 13 matched regions (hp_lehninger2)",
    ));

    let report: serde_json::Value = serde_json::from_str(&std::fs::read_to_string(&output)?)?;
    assert_eq!(report["ksize"], 12);
    assert_eq!(report["moltype"], "hp_lehninger2");
    assert_eq!(
        report["query"]["name"].as_str().unwrap().split(' ').next(),
        Some("sp|P10415|BCL2_HUMAN")
    );
    assert_eq!(report["shared_kmers"].as_array().unwrap().len(), 27);
    // The longest region comes first: BH1, 19 residues.
    assert_eq!(report["regions"][0]["query_start"], 138);
    assert_eq!(report["regions"][0]["target_start"], 162);
    assert_eq!(report["regions"][0]["length"], 19);
    Ok(())
}

#[test]
fn test_cli_pair_stdout_and_named_record() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("kmerseek")?;
    cmd.args([
        "pair",
        "--query",
        TEST_CED9_FASTA,
        "--query-name",
        "sp|P41958|CED9_CAEEL",
        "--target",
        TEST_BLC2_FASTA,
        "--ksize",
        "3",
    ]);
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("\"moltype\": \"protein20\""))
        .stdout(predicate::str::contains("\"kmer\": \"APG\""));
    Ok(())
}

/// `--extend-mismatch-penalty` grows regions past their exact seeds and reports the
/// mismatches inside; without it every region is exact and the new column reads 0.
#[test]
fn test_cli_search_extend_mismatch_penalty() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let target_index_path = temp_dir.path().join("target_index.db");
    Command::cargo_bin("kmerseek")?
        .args([
            "index",
            "--input",
            TEST_FASTA_GZ,
            "--output",
            target_index_path.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp",
        ])
        .assert()
        .success();

    let run = |extra: &[&str],
               out: &std::path::Path|
     -> Result<Vec<SearchResultCsv>, Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin("kmerseek")?;
        cmd.args([
            "search",
            "--query",
            TEST_CED9_FASTA,
            "--target",
            target_index_path.to_str().unwrap(),
            "--output",
            out.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp",
            "--min-shared-kmers",
            "0",
            "--max-pvalue",
            "0.7",
        ]);
        cmd.args(extra);
        cmd.assert().success();
        let mut reader = csv::Reader::from_path(out)?;
        Ok(reader.deserialize::<SearchResultCsv>().collect::<Result<Vec<_>, _>>()?)
    };

    let exact = run(&[], &temp_dir.path().join("exact.csv"))?;
    let extended = run(
        &["--extend-mismatch-penalty", "2", "--extend-xdrop", "8"],
        &temp_dir.path().join("extended.csv"),
    )?;

    // Off: the same 218 records (219 lines with the header) as test_cli_search_bcl2_ced9,
    // all exact. 242 before #62 joined regions that repeat a k-mer on one diagonal.
    assert_eq!(exact.len(), 218);
    assert!(exact.iter().all(|r| r.region_n_mismatches == 0));
    assert!(exact.iter().all(|r| r.region_n_shared_kmers == r.region_length - 12 + 1));

    // On: regions only grow or merge, so there are no more rows than before, at least one
    // region now spans a mismatch, and n_shared never exceeds what the span could hold.
    assert!(!extended.is_empty());
    assert!(extended.len() <= exact.len(), "{} vs {}", extended.len(), exact.len());
    assert!(extended.iter().any(|r| r.region_n_mismatches > 0));
    for r in &extended {
        assert!(r.region_length >= 12);
        assert!(r.region_n_shared_kmers <= r.region_length - 12 + 1);
        assert_eq!(r.region_subseq.len() as u32, r.region_length);
        assert_eq!(r.target_subseq.len() as u32, r.region_length);
    }
    let bcl2_exact: Vec<_> =
        exact.iter().filter(|r| r.target_name.contains("BCL2_HUMAN")).collect();
    let bcl2_ext: Vec<_> =
        extended.iter().filter(|r| r.target_name.contains("BCL2_HUMAN")).collect();
    assert!(bcl2_ext.len() <= bcl2_exact.len());
    assert!(
        bcl2_ext.iter().map(|r| r.region_length).max()
            >= bcl2_exact.iter().map(|r| r.region_length).max(),
        "the longest BCL2/CED9 region can only get longer"
    );
    Ok(())
}

/// `kmerseek index` fits the lambda scale and K for its penalty and X-drop and stores them;
/// `kmerseek search` reads them back, lets `--ka-k` override them, refuses a penalty that
/// was never fitted when `--ka-queries 0` forbids fitting one now, and fits one otherwise.
#[test]
fn test_cli_ka_fit_at_index_time_is_reused() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let index_path = temp_dir.path().join("target_index.db");
    Command::cargo_bin("kmerseek")?
        .args([
            "index",
            "--input",
            TEST_FASTA_GZ,
            "--output",
            index_path.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp",
            "--ka-queries",
            "25",
        ])
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "Fitting Karlin-Altschul lambda and K on 25 database sequences (penalty 2, X-drop 8), \
             censored against the same sequences shuffled-dipeptide",
        ))
        .stderr(predicate::str::contains(
            "Closed form at the database's own match probability 0.500: K 0.1631",
        ))
        .stderr(predicate::str::contains(
            "fitted now: 25 database queries, 9838 regions; slope 0.806 per nat of lambda_pair S \
             (1 = closed form holds; closed form 0.481 at the database's match probability \
             0.500), K 0.0115, fit on x 7.5..11.5, rms 0.086; shuffled-dipeptide reference \
             slope 0.760 over the same bins",
        ))
        .stderr(predicate::str::contains(
            "Stored in the index for --extend-mismatch-penalty 2 --extend-xdrop 8",
        ));

    let search =
        |extra: &[&str]| -> Result<assert_cmd::assert::Assert, Box<dyn std::error::Error>> {
            let mut cmd = Command::cargo_bin("kmerseek")?;
            cmd.args([
                "search",
                "--query",
                TEST_CED9_FASTA,
                "--target",
                index_path.to_str().unwrap(),
                "--output",
                temp_dir.path().join("out.csv").to_str().unwrap(),
            ]);
            cmd.args(extra);
            Ok(cmd.assert())
        };

    search(&["--extend-mismatch-penalty", "2"])?.success().stderr(predicate::str::contains(
        "Karlin-Altschul: K 0.0115, lambda scale 0.806 (stored in the index: 25 database queries, 9838 regions",
    ));
    search(&["--extend-mismatch-penalty", "2", "--ka-k", "0.03"])?.success().stderr(
        predicate::str::contains(
            "Karlin-Altschul: K 0.0300, lambda scale 1.000 (--ka-k, closed-form lambda)",
        ),
    );
    search(&["--extend-mismatch-penalty", "3", "--ka-queries", "0"])?
        .failure()
        .stderr(predicate::str::contains("no Karlin-Altschul fit for penalty 3, X-drop 8"));
    search(&["--extend-mismatch-penalty", "3", "--ka-queries", "25"])?.success().stderr(
        predicate::str::contains(
            "Karlin-Altschul: K 0.1264, lambda scale 0.962 (fitted now: 25 database queries, 9854 regions; slope 0.962",
        ),
    );
    Ok(())
}

/// `kmerseek calibrate` adds a Karlin-Altschul fit for one more (penalty, X-drop) pair to
/// a finished index, so a search at that penalty reads it instead of refitting; a pair
/// nobody fitted is still refused when `--ka-queries 0` forbids fitting one now.
#[test]
fn test_cli_calibrate_stores_a_second_penalty() -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let index_path = temp_dir.path().join("target_index.db");
    Command::cargo_bin("kmerseek")?
        .args([
            "index",
            "--input",
            TEST_FASTA_GZ,
            "--output",
            index_path.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp",
            "--ka-queries",
            "0",
        ])
        .assert()
        .success();

    Command::cargo_bin("kmerseek")?
        .args([
            "calibrate",
            "--target",
            index_path.to_str().unwrap(),
            "--extend-mismatch-penalty",
            "1.5",
            "--extend-xdrop",
            "6",
            "--ka-queries",
            "25",
        ])
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "Stored in the index for --extend-mismatch-penalty 1.5 --extend-xdrop 6",
        ));

    let search =
        |extra: &[&str]| -> Result<assert_cmd::assert::Assert, Box<dyn std::error::Error>> {
            let mut cmd = Command::cargo_bin("kmerseek")?;
            cmd.args([
                "search",
                "--query",
                TEST_CED9_FASTA,
                "--target",
                index_path.to_str().unwrap(),
                "--output",
                temp_dir.path().join("out.csv").to_str().unwrap(),
                "--ka-queries",
                "0",
            ]);
            cmd.args(extra);
            Ok(cmd.assert())
        };
    search(&["--extend-mismatch-penalty", "1.5", "--extend-xdrop", "6"])?
        .success()
        .stderr(predicate::str::contains("(stored in the index: 25 database queries"));
    search(&["--extend-mismatch-penalty", "2"])?
        .failure()
        .stderr(predicate::str::contains("no Karlin-Altschul fit for penalty 2, X-drop 8"));

    // --ka-queries 0 is a request to fit nothing, which calibrate cannot honour.
    Command::cargo_bin("kmerseek")?
        .args(["calibrate", "--target", index_path.to_str().unwrap(), "--ka-queries", "0"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("--ka-queries must be at least 1"));
    Ok(())
}
