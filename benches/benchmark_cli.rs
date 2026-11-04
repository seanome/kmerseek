use criterion::{criterion_group, criterion_main, Criterion};
use std::path::Path;
use std::process::Command;
use std::sync::Once;
use tempfile::tempdir;

/// The `kmerseek` binary cargo built for this bench run, not whichever one is on PATH.
const KMERSEEK: &str = env!("CARGO_BIN_EXE_kmerseek");

/// Criterion refuses a sample size below 10.
const SAMPLE_SIZE: usize = 10;

/// `/usr/bin/time` flag that prints peak memory: `-l` on macOS (BSD time), `-v` on Linux
/// (GNU time).
#[cfg(target_os = "macos")]
const TIME_FLAG: &str = "-l";
#[cfg(not(target_os = "macos"))]
const TIME_FLAG: &str = "-v";

/// Arguments for `kmerseek index` on one FASTA file.
fn index_args(input: &str, output: &Path, ksize: u32, alphabet: &str) -> Vec<String> {
    let args = ["index", "--input", input, "--output", output.to_str().unwrap()];
    let mut args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
    args.extend(["--ksize".into(), ksize.to_string(), "--alphabet".into(), alphabet.into()]);
    args
}

/// Run `kmerseek index` and panic if it fails.
fn run_index(input: &str, output: &Path, ksize: u32, alphabet: &str) {
    let status = Command::new(KMERSEEK)
        .args(index_args(input, output, ksize, alphabet))
        .status()
        .expect("Failed to execute kmerseek");
    assert!(status.success(), "kmerseek failed");
}

/// Peak resident memory in KiB from `/usr/bin/time` output. macOS prints
/// "<bytes>  maximum resident set size"; GNU time prints
/// "Maximum resident set size (kbytes): <KiB>".
fn parse_peak_rss_kib(time_output: &str) -> Option<u64> {
    let line = time_output
        .lines()
        .find(|line| line.to_lowercase().contains("maximum resident set size"))?;
    if cfg!(target_os = "macos") {
        line.split_whitespace().next()?.parse::<u64>().ok().map(|bytes| bytes / 1024)
    } else {
        line.rsplit(':').next()?.trim().parse().ok()
    }
}

/// Total size in bytes of every file under `path`. The index is a directory.
fn dir_size_bytes(path: &Path) -> u64 {
    std::fs::read_dir(path)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .map(|p| if p.is_dir() { dir_size_bytes(&p) } else { p.metadata().unwrap().len() })
        .sum()
}

fn benchmark_rust_cli_performance(c: &mut Criterion) {
    let mut group = c.benchmark_group("Rust CLI Indexing Performance");
    group.sample_size(SAMPLE_SIZE);

    // Test data files
    let test_files = [
        ("ced9.fasta", "tests/testdata/fasta/ced9.fasta"),
        (
            "bcl2_first25.fasta.gz",
            "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz",
        ),
    ];

    for (name, file_path) in test_files.iter() {
        // Test different k-mer sizes
        for ksize in [5, 10, 15] {
            // Test different alphabets
            for alphabet in ["protein20", "hp_lehninger2", "dayhoff6"] {
                let benchmark_name = format!("{}_{}_k{}", name, alphabet, ksize);

                group.bench_function(&benchmark_name, |b| {
                    b.iter(|| {
                        let temp_dir = tempdir().unwrap();
                        let output_path = temp_dir.path().join("test_output.db");
                        run_index(file_path, &output_path, ksize, alphabet);
                    });
                });
            }
        }
    }

    group.finish();
}

fn benchmark_memory_usage(c: &mut Criterion) {
    let mut group = c.benchmark_group("Memory Usage");
    group.sample_size(SAMPLE_SIZE);

    let test_file = "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";
    // Print the peak memory once, not on every iteration. It is measured inside `iter` so
    // that `--list` does not run kmerseek.
    let report = Once::new();

    group.bench_function("rust_memory_usage", |b| {
        b.iter(|| {
            let temp_dir = tempdir().unwrap();
            let output_path = temp_dir.path().join("test_output.db");

            let output = Command::new("/usr/bin/time")
                .arg(TIME_FLAG)
                .arg(KMERSEEK)
                .args(index_args(test_file, &output_path, 10, "hp_lehninger2"))
                .output()
                .expect("Failed to execute /usr/bin/time");
            assert!(output.status.success(), "kmerseek failed");

            let stderr = String::from_utf8_lossy(&output.stderr);
            let peak_kib = parse_peak_rss_kib(&stderr).expect("no peak memory in time output");
            report.call_once(|| eprintln!("Rust peak memory: {} KiB", peak_kib));
        });
    });

    group.finish();
}

fn benchmark_output_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("Output Size");
    group.sample_size(SAMPLE_SIZE);

    let test_file = "tests/testdata/fasta/ced9.fasta";
    let report = Once::new();

    group.bench_function("rust_output_size", |b| {
        b.iter(|| {
            let temp_dir = tempdir().unwrap();
            let output_path = temp_dir.path().join("test_output.db");
            run_index(test_file, &output_path, 10, "protein20");

            let size_mb = dir_size_bytes(&output_path) as f64 / (1024.0 * 1024.0);
            report.call_once(|| eprintln!("Rust output size: {:.2} MB", size_mb));
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    benchmark_rust_cli_performance,
    benchmark_memory_usage,
    benchmark_output_size,
);
criterion_main!(benches);
