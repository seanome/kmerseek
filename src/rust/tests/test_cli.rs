use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::tempdir;

use approx::assert_relative_eq;

use crate::alphabets::Alphabet;
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
        assert_eq!(rows, 362, "ced9 against the 25-sequence bcl2 index at k=12");
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

    // Off: the same 242 records (243 lines with the header) as test_cli_search_bcl2_ced9,
    // all exact.
    assert_eq!(exact.len(), 242);
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

/// One calibration query gives too few score bins to fit, so nothing is stored; three
/// give a fit that leans on the seed end of the curve, which the index step says out
/// loud and writes out as a survival curve.
#[test]
fn test_cli_ka_fit_on_few_queries_warns_and_writes_the_survival_curve(
) -> Result<(), Box<dyn std::error::Error>> {
    let temp_dir = tempdir()?;
    let index = |n_queries: &str, null: &str, survival: &std::path::Path| {
        let out = temp_dir.path().join(format!("index_{n_queries}_{null}.db"));
        let mut cmd = Command::cargo_bin("kmerseek")?;
        cmd.args([
            "index",
            "--input",
            TEST_FASTA_GZ,
            "--output",
            out.to_str().unwrap(),
            "--ksize",
            "12",
            "--alphabet",
            "hp",
            "--ka-queries",
            n_queries,
            "--ka-null",
            null,
            "--ka-survival-out",
            survival.to_str().unwrap(),
        ]);
        Ok::<_, Box<dyn std::error::Error>>(cmd.assert())
    };

    let unfit = temp_dir.path().join("unfit.csv");
    index("1", "database", &unfit)?.success().stderr(predicate::str::contains(
        "1 queries gave only 420 regions, too few score bins to fit; nothing stored.",
    ));
    assert!(!unfit.exists(), "no fit, no curve to write");

    // Shuffled queries need no reference to censor against, so none is announced.
    index("3", "shuffled", &temp_dir.path().join("shuffled.csv"))?
        .success()
        .stderr(predicate::str::contains(
            "Fitting Karlin-Altschul lambda and K on 3 shuffled sequences (penalty 2, X-drop 8)...",
        ))
        .stderr(predicate::str::contains(
            "fitted now: 3 shuffled queries, 789 regions; slope 0.896 per nat of lambda_pair S",
        ))
        .stderr(predicate::str::contains("K 0.0191, fit on x 6.5..8.5, rms 0.070\n"));

    let survival = temp_dir.path().join("survival.csv");
    index("3", "database", &survival)?
        .success()
        .stderr(predicate::str::contains(
            "Fitting Karlin-Altschul lambda and K on 3 database sequences (penalty 2, X-drop 8), \
             censored against the same sequences shuffled-dipeptide...",
        ))
        .stderr(predicate::str::contains(
            "fitted now: 3 database queries, 903 regions; slope 0.730 per nat of lambda_pair S \
             (1 = closed form holds; closed form 0.481 at the database's match probability \
             0.500), K 0.0082, fit on x 6.0..9.0, rms 0.184; shuffled-dipeptide reference \
             slope 0.712 over the same bins",
        ))
        .stderr(predicate::str::contains(
            "WARNING: the fit has only 6 bins (x 6.0..9.0) below the relatives at x none.",
        ))
        .stderr(predicate::str::contains(format!(
            "Survival curve written to {}",
            survival.display()
        )));

    let mut reader = csv::Reader::from_path(&survival)?;
    assert_eq!(
        reader.headers()?.iter().collect::<Vec<_>>(),
        [
            "x",
            "n_regions_at_least",
            "fitted_n_regions_at_least",
            "reference_n_regions_at_least",
            "in_fit",
            "slope",
            "k",
            "lambda_analytic",
            "match_probability",
            "null",
            "reference",
            "mismatch_penalty",
            "xdrop",
            "n_queries",
            "query_residues",
            "database_kmers",
            "bin_width",
        ]
    );
    let rows: Vec<csv::StringRecord> = reader.records().collect::<Result<_, _>>()?;
    assert_eq!(rows.len(), 62, "one row per half-nat bin of x from 4.5 to 35");
    assert_eq!((&rows[0][0], &rows[61][0]), ("4.500", "35.000"));
    // The lowest bin holds every region (903 seeds and up); the fit window is x 6.0..9.0.
    assert_eq!(&rows[0][1], "903");
    assert_eq!(&rows[0][4], "false");
    let in_fit: Vec<&str> = rows.iter().filter(|r| &r[4] == "true").map(|r| &r[0]).collect();
    assert_eq!(in_fit, ["6.000", "6.500", "7.000", "7.500", "8.000", "8.500"]);
    // Lowest bin: 903 regions, the fitted line far above at 1955.245 (the seed floor bends
    // the curve there), 864 in the shuffled-dipeptide reference.
    assert_eq!((&rows[0][1], &rows[0][2], &rows[0][3]), ("903", "1955.245", "864"));
    // The top bin is one region, below the fitted line's resolution, and the reference
    // never got there.
    assert_eq!((&rows[61][1], &rows[61][2], &rows[61][3]), ("1", "0.000", ""));
    // Fit constants repeat on every row.
    let constants = |r: &csv::StringRecord| r.iter().skip(5).map(str::to_owned).collect::<Vec<_>>();
    assert!(rows.iter().all(|r| constants(r) == constants(&rows[0])));
    assert_eq!(
        constants(&rows[0]),
        [
            "0.7304804242524211",
            "0.008235243224364983",
            "0.48120853696601995",
            "0.5000011360087127",
            "database",
            "shuffled-dipeptide",
            "2",
            "8",
            "3",
            "762",
            "8340",
            "0.5"
        ]
    );
    Ok(())
}
