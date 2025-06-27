use criterion::{criterion_group, criterion_main, Criterion};
use kmerseek::encoding::{encode_kmer, encode_kmer_with_encoding_fn, get_encoding_fn_from_moltype};
use kmerseek::index::ProteomeIndex;
use kmerseek::signature::SEED;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use tempfile::tempdir;

// Test protein sequence
const TEST_PROTEIN: &str = "PLANTYANDANIMALGENQMESCOFFEE";

// Test protein sequences with different characteristics for benchmarking
const TEST_PROTEIN_WITH_AMBIGUOUS: &str = "PLANTYANDANIMALGENQMESCOFFEEBZJ";
const TEST_PROTEIN_WITH_SPECIAL: &str = "PLANTYANDANIMALGENQMESCOFFEEXUO";
const TEST_PROTEIN_WITH_STOP: &str = "PLANTYANDANIMALGENQMESCOFFEE*EXTRA";

const KSIZES: [u32; 3] = [5, 10, 20];

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

fn benchmark_create_protein_signature(c: &mut Criterion) {
    for moltype in ["protein", "hp", "dayhoff"] {
        for ksize in KSIZES {
            let (index, _) = setup_test_index(ksize, moltype);

            // Benchmark standard protein sequence
            c.bench_function(
                &format!("create_protein_signature_standard_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        index.create_protein_signature(TEST_PROTEIN, "test_protein").unwrap();
                    })
                },
            );

            // Benchmark protein sequence with ambiguous amino acids
            c.bench_function(
                &format!("create_protein_signature_ambiguous_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        index
                            .create_protein_signature(
                                TEST_PROTEIN_WITH_AMBIGUOUS,
                                "test_protein_ambiguous",
                            )
                            .unwrap();
                    })
                },
            );

            // Benchmark protein sequence with special amino acids (X, U, O)
            c.bench_function(
                &format!("create_protein_signature_special_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        index
                            .create_protein_signature(
                                TEST_PROTEIN_WITH_SPECIAL,
                                "test_protein_special",
                            )
                            .unwrap();
                    })
                },
            );

            // Benchmark protein sequence with stop codon
            c.bench_function(
                &format!("create_protein_signature_stop_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        index
                            .create_protein_signature(TEST_PROTEIN_WITH_STOP, "test_protein_stop")
                            .unwrap();
                    })
                },
            );
        }
    }
}

fn benchmark_proteome_index_encode_kmer(c: &mut Criterion) {
    for moltype in ["protein", "hp", "dayhoff"] {
        for ksize in KSIZES {
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
        for ksize in KSIZES {
            c.bench_function(&format!("encodings_encode_kmer_{}_{}", moltype, ksize), |b| {
                b.iter(|| encode_kmer(&TEST_PROTEIN[..ksize as usize], moltype))
            });
        }
    }
}

fn benchmark_encodings_encode_kmer_with_encoding_fn(c: &mut Criterion) {
    for moltype in ["protein", "hp", "dayhoff"] {
        let encoding_fn = get_encoding_fn_from_moltype(moltype).unwrap();
        for ksize in KSIZES {
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
        for ksize in KSIZES {
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

fn benchmark_process_fasta(c: &mut Criterion) {
    let temp_dir = tempdir().unwrap();
    // Create a temporary FASTA file for testing with distinct sequences
    let fasta_content = ">sp|O43236|SEPT4_HUMAN Septin-4 OS=Homo sapiens OX=9606 GN=SEPTIN4 PE=1 SV=1
MDRSLGWQGNSVPEDRTEAGIKRFLEDTTDDGELSKFVKDFSGNASCHPPEAKTWASRPQ
VPEPRPQAPDLYDDDLEFRPPSRPQSSDNQQYFCAPAPLSPSARPRSPWGKLDPYDSSED
DKEYVGFATLPNQVHRKSVKKGFDFTLMVAGESGLGKSTLVNSLFLTDLYRDRKLLGAEE
RIMQTVEITKHAVDIEEKGVRLRLTIVDTPGFGDAVNNTECWKPVAEYIDQQFEQYFRDE
SGLNRKNIQDNRVHCCLYFISPFGHGLRPLDVEFMKALHQRVNIVPILAKADTLTPPEVD
HKKRKIREEIEHFGIKIYQFPDCDSDEDEDFKLQDQALKESIPFAVIGSNTVVEARGRRV
RGRLYPWGIVEVENPGHCDFVKLRTMLVRTHMQDLKDVTRETHYENYRAQCIQSMTRLVV
KERNRNKLTRESGTDFPIPAVPPGTDPETEKLIREKDEELRRMQEMLHKIQKQMKENY
>sp|O43521|B2L11_HUMAN Bcl-2-like protein 11 OS=Homo sapiens OX=9606 GN=BCL2L11 PE=1 SV=1
MAKQPSDVSSECDREGRQLQPAERPPQLRPGAPTSLQTEPQGNPEGNHGGEGDSCPHGSP
QGPLAPPASPGPFATRSPLFIFMRRSSLLSRSSSGYFSFDTDRSPAPMSCDKSTQTPSPP
CQAFNHYLSAMASMRQAEPADMRPEIWIAQELRRIGDEFNAYYARRVFLNNYQAAEDHPR
MVILRLLRYIVRLVWRMH
>sp|O60238|BNI3L_HUMAN BCL2/adenovirus E1B 19 kDa protein-interacting protein 3-like OS=Homo sapiens OX=9606 GN=BNIP3L PE=1 SV=1
MSSHLVEPPPPLHNNNNNCEENEQSLPPPAGLNSSWVELPMNSSNGNDNGNGKNGGLEHV
PSSSSIHNGDMEKILLDAQHESGQSSSRGSSHCDSPSPQEDGQIMFDVEMHTSRDHSSQS
EEEVVEGEKEVEALKKSADWVSDWSSRPENIPPKEFHFRHPKRSVSLSMRKSGAMKKGGI
FSAEFLKVFIPSLFLSHVLALGLGIYIGKRLSTPSASTY";
    let fasta_path = temp_dir.path().join("test.fasta");
    std::fs::write(&fasta_path, fasta_content).unwrap();
    for moltype in ["protein", "hp", "dayhoff"] {
        for ksize in KSIZES {
            let (index, _) = setup_test_index(ksize, moltype);

            c.bench_function(&format!("process_fasta_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    index.process_fasta(&fasta_path, 0).unwrap();
                })
            });
        }
    }
}

criterion_group!(
    benches,
    benchmark_create_protein_signature,
    benchmark_proteome_index_encode_kmer,
    benchmark_encodings_encode_kmer,
    benchmark_encodings_encode_kmer_with_encoding_fn,
    benchmark_process_protein_kmers,
    benchmark_process_fasta,
);
criterion_main!(benches);
