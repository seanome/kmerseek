use criterion::{criterion_group, criterion_main, Criterion};
use std::process::Command;
use tempfile::tempdir;

fn benchmark_rust_cli_performance(c: &mut Criterion) {
    let mut group = c.benchmark_group("Rust CLI Indexing Performance");
    group.sample_size(10);

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
            // Test different encodings
            for encoding in ["protein", "hp", "dayhoff"] {
                let benchmark_name = format!("{}_{}_k{}", name, encoding, ksize);

                group.bench_function(&benchmark_name, |b| {
                    b.iter(|| {
                        let temp_dir = tempdir().unwrap();
                        let output_path = temp_dir.path().join("test_output.db");

                        let status = Command::new("kmerseek-rust")
                            .args([
                                "index",
                                "--input",
                                file_path,
                                "--output",
                                output_path.to_str().unwrap(),
                                "--ksize",
                                &ksize.to_string(),
                                "--encoding",
                                encoding,
                            ])
                            .status()
                            .expect("Failed to execute kmerseek-rust");

                        assert!(status.success(), "kmerseek-rust failed");
                    });
                });
            }
        }
    }

    group.finish();
}

fn benchmark_memory_usage(c: &mut Criterion) {
    let mut group = c.benchmark_group("Memory Usage");
    group.sample_size(5);

    let test_file = "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";

    group.bench_function("rust_memory_usage", |b| {
        b.iter(|| {
            let temp_dir = tempdir().unwrap();
            let output_path = temp_dir.path().join("test_output.db");

            // Use /usr/bin/time to measure memory usage on Unix systems
            let output = Command::new("sh")
                .arg("-c")
                .arg(format!(
                    "/usr/bin/time -l kmerseek-rust index --input {} --output {} --ksize 10 --encoding hp 2>&1",
                    test_file,
                    output_path.to_str().unwrap()
                ))
                .output()
                .expect("Failed to execute command");

            // Parse memory usage from output (macOS format)
            let output_str = String::from_utf8_lossy(&output.stdout);
            if let Some(memory_line) = output_str.lines().find(|line| line.contains("maximum resident set size")) {
                println!("Rust memory usage: {}", memory_line);
            }
        });
    });

    group.finish();
}

fn benchmark_output_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("Output Size");
    group.sample_size(1);

    let test_file = "tests/testdata/fasta/ced9.fasta";

    group.bench_function("rust_output_size", |b| {
        b.iter(|| {
            let temp_dir = tempdir().unwrap();
            let output_path = temp_dir.path().join("test_output.db");

            let status = Command::new("kmerseek-rust")
                .args([
                    "index",
                    "--input",
                    test_file,
                    "--output",
                    output_path.to_str().unwrap(),
                    "--ksize",
                    "10",
                    "--encoding",
                    "protein",
                ])
                .status()
                .expect("Failed to execute kmerseek-rust");

            assert!(status.success(), "kmerseek-rust failed");

            // Measure output size
            let metadata = std::fs::metadata(&output_path).unwrap();
            let size_mb = metadata.len() as f64 / (1024.0 * 1024.0);
            println!("Rust output size: {:.2} MB", size_mb);
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
