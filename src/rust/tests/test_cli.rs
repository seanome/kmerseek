use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::tempdir;

#[test]
fn test_cli_help() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("kmerseek-rust")?;
    cmd.arg("--help");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("A tool for indexing protein sequences from FASTA files"));

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
        "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz",
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
            "tests/testdata/fasta/ced9.fasta",
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
