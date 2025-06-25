use criterion::{criterion_group, criterion_main, Criterion};
use kmerseek::encoding::encode_kmer;
use kmerseek::index::ProteomeIndex;
use kmerseek::kmer_signature::SEED;
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
            c.bench_function(&format!("proteome_index_encode_kmer_{}_{}", moltype, ksize), |b| {
                b.iter(|| index.encode_kmer(&TEST_PROTEIN[..ksize as usize]))
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

fn benchmark_process_protein_kmers(c: &mut Criterion) {
    for moltype in ["protein", "hp", "dayhoff"] {
        for ksize in [5, 10, 15, 20] {
            let (index, _) = setup_test_index(ksize, moltype);

            // Create a test signature
            let mut minhash = sourmash::sketch::minhash::KmerMinHash::new(
                1, // scaled
                ksize,
                kmerseek::encoding::get_hash_function_from_moltype(moltype).unwrap(),
                SEED,
                true, // track_abundance
                0,    // num (use scaled instead)
            );
            minhash.add_sequence(TEST_PROTEIN.as_bytes(), true).unwrap();

            let small_sig = sourmash_plugin_branchwater::utils::multicollection::SmallSignature {
                location: "test".to_string(),
                name: "test".to_string(),
                md5sum: "test".to_string(),
                minhash,
            };

            c.bench_function(&format!("process_protein_kmers_{}_{}", moltype, ksize), |b| {
                b.iter(|| index.process_protein_kmers(TEST_PROTEIN, &small_sig))
            });
        }
    }
}

criterion_group!(
    benches,
    benchmark_proteome_index_encode_kmer,
    benchmark_encodings_encode_kmer,
    benchmark_process_protein_kmers
);
criterion_main!(benches);
