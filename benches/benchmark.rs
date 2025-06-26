use criterion::{criterion_group, criterion_main, Criterion};
use kmerseek::encoding::{encode_kmer, encode_kmer_with_encoding_fn, get_encoding_fn_from_moltype};
use kmerseek::index::ProteomeIndex;
use kmerseek::signature::SEED;
use sourmash::signature::SigsTrait;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use tempfile::tempdir;

// Test protein sequence
const TEST_PROTEIN: &str = "PLANTYANDANIMALGENQMESCOFFEE";

fn setup_test_index(ksize: u32, moltype: &str) -> (ProteomeIndex, PathBuf) {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join(format!("proteins_{}_{}.fasta", moltype, ksize));
    let db_path = temp_dir.path().join(format!("db_{}_{}", moltype, ksize));

    // Create a simple test FASTA file
    let mut file = File::create(&fasta_path).unwrap();
    writeln!(file, ">test_protein\n{}", TEST_PROTEIN).unwrap();

    let index = ProteomeIndex::new(db_path, ksize, 1, moltype, SEED).unwrap();
    (index, fasta_path)
}

fn benchmark_proteome_index_encode_kmer(c: &mut Criterion) {
    for moltype in ["protein", "hp", "dayhoff"] {
        for ksize in [5, 10, 15, 20] {
            let (index, _) = setup_test_index(ksize, moltype);
            let encoding_fn = kmerseek::encoding::get_encoding_fn_from_moltype(moltype).unwrap();
            c.bench_function(&format!("proteome_index_encode_kmer_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    kmerseek::encoding::encode_kmer_with_encoding_fn(
                        &TEST_PROTEIN[..ksize as usize],
                        encoding_fn,
                    )
                })
            });
        }
    }
}

fn benchmark_encodings_encode_kmer(c: &mut Criterion) {
    for moltype in ["protein", "hp", "dayhoff"] {
        for ksize in [5, 10, 15, 20] {
            c.bench_function(&format!("encodings_encode_kmer_{}_{}", moltype, ksize), |b| {
                b.iter(|| encode_kmer(&TEST_PROTEIN[..ksize as usize], moltype))
            });
        }
    }
}

fn benchmark_encodings_encode_kmer_with_encoding_fn(c: &mut Criterion) {
    for moltype in ["protein", "hp", "dayhoff"] {
        let encoding_fn = get_encoding_fn_from_moltype(moltype).unwrap();
        for ksize in [5, 10, 15, 20] {
            c.bench_function(
                &format!("encodings_encode_kmer_with_encoding_fn_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        encode_kmer_with_encoding_fn(&TEST_PROTEIN[..ksize as usize], encoding_fn)
                    })
                },
            );
        }
    }
}

fn benchmark_process_protein_kmers(c: &mut Criterion) {
    for moltype in ["protein", "hp", "dayhoff"] {
        for ksize in [5, 10, 15, 20] {
            let (index, _) = setup_test_index(ksize, moltype);

            // Create a test protein signature
            let mut protein_sig = kmerseek::signature::ProteinSignature::new(
                "test_protein",
                ksize,
                1, // scaled
                moltype,
                kmerseek::signature::SEED,
            )
            .unwrap();

            // Add the protein sequence
            protein_sig.add_protein(TEST_PROTEIN.as_bytes()).unwrap();

            c.bench_function(&format!("process_protein_kmers_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    let mut sig = protein_sig.clone();
                    index.process_kmers(TEST_PROTEIN, &mut sig).unwrap();
                })
            });
        }
    }
}

criterion_group!(
    benches,
    benchmark_proteome_index_encode_kmer,
    benchmark_encodings_encode_kmer,
    benchmark_encodings_encode_kmer_with_encoding_fn,
    benchmark_process_protein_kmers,
);
criterion_main!(benches);
