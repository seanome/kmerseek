use criterion::{criterion_group, criterion_main, Criterion};
use kmerseek::encoding::{encode_by_moltype, encode_with_fn, get_encoding_fn_from_moltype};
use kmerseek::index::ProteomeIndex;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::mem;
use std::path::PathBuf;
use std::time::Instant;
use tempfile::tempdir;
use needletail;

// Test protein sequence
const TEST_PROTEIN: &str = "PLANTYANDANIMALGENQMESCOFFEE";

// Test protein sequences with different characteristics for benchmarking
const TEST_PROTEIN_WITH_AMBIGUOUS: &str = "PLANTYANDANIMALGENQMESCOFFEEBZJ";
const TEST_PROTEIN_WITH_SPECIAL: &str = "PLANTYANDANIMALGENQMESCOFFEEXUO";
const TEST_PROTEIN_WITH_STOP: &str = "PLANTYANDANIMALGENQMESCOFFEE*EXTRA";

const KSIZES: [u32; 3] = [5, 10, 20];
const MOLTYPES: [&str; 3] = ["protein", "hp", "dayhoff"];

fn setup_test_index(ksize: u32, moltype: &str) -> (ProteomeIndex, PathBuf) {
    let temp_dir = tempdir().unwrap();
    let fasta_path = temp_dir.path().join(format!("proteins_{}_{}.fasta", moltype, ksize));
    let db_path = temp_dir.path().join(format!("db_{}_{}", moltype, ksize));

    // Create a simple test FASTA file
    let mut file = File::create(&fasta_path).unwrap();
    writeln!(file, ">test_protein\n{}", TEST_PROTEIN).unwrap();

    let index = ProteomeIndex::new(db_path, ksize, 1, moltype, false).unwrap();
    (index, fasta_path)
}

fn benchmark_create_protein_signature(c: &mut Criterion) {
    for moltype in MOLTYPES {
        for ksize in KSIZES {
            let (index, _) = setup_test_index(ksize, moltype);

            // Benchmark standard protein sequence
            c.bench_function(
                &format!("create_protein_signature_standard_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        // Record start time for CPU measurement
                        let start_time = Instant::now();

                        // Create protein signature
                        let signature =
                            index.create_protein_signature(TEST_PROTEIN, "test_protein").unwrap();

                        // Record end time
                        let end_time = Instant::now();
                        let _cpu_time = end_time.duration_since(start_time);

                        // Measure basic memory usage (stack size)
                        let _memory_used = mem::size_of_val(&signature);
                    })
                },
            );

            // Benchmark protein sequence with ambiguous amino acids
            c.bench_function(
                &format!("create_protein_signature_ambiguous_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        // Record start time for CPU measurement
                        let start_time = Instant::now();

                        // Create protein signature
                        let signature = index
                            .create_protein_signature(
                                TEST_PROTEIN_WITH_AMBIGUOUS,
                                "test_protein_ambiguous",
                            )
                            .unwrap();

                        // Record end time
                        let end_time = Instant::now();
                        let _cpu_time = end_time.duration_since(start_time);

                        // Measure basic memory usage (stack size)
                        let _memory_used = mem::size_of_val(&signature);
                    })
                },
            );

            // Benchmark protein sequence with special amino acids (X, U, O)
            c.bench_function(
                &format!("create_protein_signature_special_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        // Record start time for CPU measurement
                        let start_time = Instant::now();

                        // Create protein signature
                        let signature = index
                            .create_protein_signature(
                                TEST_PROTEIN_WITH_SPECIAL,
                                "test_protein_special",
                            )
                            .unwrap();

                        // Record end time
                        let end_time = Instant::now();
                        let _cpu_time = end_time.duration_since(start_time);

                        // Measure basic memory usage (stack size)
                        let _memory_used = mem::size_of_val(&signature);
                    })
                },
            );

            // Benchmark protein sequence with stop codon
            c.bench_function(
                &format!("create_protein_signature_stop_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        // Record start time for CPU measurement
                        let start_time = Instant::now();

                        // Create protein signature
                        let signature = index
                            .create_protein_signature(TEST_PROTEIN_WITH_STOP, "test_protein_stop")
                            .unwrap();

                        // Record end time
                        let end_time = Instant::now();
                        let _cpu_time = end_time.duration_since(start_time);

                        // Measure basic memory usage (stack size)
                        let _memory_used = mem::size_of_val(&signature);
                    })
                },
            );
        }
    }
}

fn benchmark_proteome_index_encode_kmer(c: &mut Criterion) {
    for moltype in MOLTYPES {
        for ksize in KSIZES {
            let (_index, _) = setup_test_index(ksize, moltype);
            let encoding_fn = kmerseek::encoding::get_encoding_fn_from_moltype(moltype).unwrap();
            c.bench_function(&format!("proteome_index_encode_kmer_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    // Record start time for CPU measurement
                    let start_time = Instant::now();

                    // Encode kmer
                    let encoded = kmerseek::encoding::encode_with_fn(
                        &TEST_PROTEIN[..ksize as usize],
                        encoding_fn,
                    )
                    .unwrap();

                    // Record end time
                    let end_time = Instant::now();
                    let _cpu_time = end_time.duration_since(start_time);

                    // Measure basic memory usage (stack size)
                    let _memory_used = mem::size_of_val(&encoded);
                })
            });
        }
    }
}

fn benchmark_encodings_encode_kmer(c: &mut Criterion) {
    for moltype in MOLTYPES {
        for ksize in KSIZES {
            c.bench_function(&format!("encodings_encode_kmer_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    // Record start time for CPU measurement
                    let start_time = Instant::now();

                    // Encode kmer
                    let encoded =
                        encode_by_moltype(&TEST_PROTEIN[..ksize as usize], moltype).unwrap();

                    // Record end time
                    let end_time = Instant::now();
                    let _cpu_time = end_time.duration_since(start_time);

                    // Measure basic memory usage (stack size)
                    let _memory_used = mem::size_of_val(&encoded);
                })
            });
        }
    }
}

fn benchmark_encodings_encode_kmer_with_encoding_fn(c: &mut Criterion) {
    for moltype in MOLTYPES {
        let encoding_fn = get_encoding_fn_from_moltype(moltype).unwrap();
        for ksize in KSIZES {
            c.bench_function(
                &format!("encodings_encode_kmer_with_encoding_fn_{}_{}", moltype, ksize),
                |b| {
                    b.iter(|| {
                        // Record start time for CPU measurement
                        let start_time = Instant::now();

                        // Encode kmer
                        let encoded =
                            encode_with_fn(&TEST_PROTEIN[..ksize as usize], encoding_fn).unwrap();

                        // Record end time
                        let end_time = Instant::now();
                        let _cpu_time = end_time.duration_since(start_time);

                        // Measure basic memory usage (stack size)
                        let _memory_used = mem::size_of_val(&encoded);
                    })
                },
            );
        }
    }
}

fn benchmark_process_protein_kmers(c: &mut Criterion) {
    for moltype in MOLTYPES {
        for ksize in KSIZES {
            let (index, _) = setup_test_index(ksize, moltype);

            // Create a test protein signature
            let mut protein_sig = kmerseek::sketch::ProteinSketch::new(
                "test_protein",
                ksize,
                1, // scaled
                moltype,
            )
            .unwrap();

            // Add the protein sequence
            protein_sig.add_protein(TEST_PROTEIN, true).unwrap();

            c.bench_function(&format!("process_protein_kmers_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    // Record start time for CPU measurement
                    let start_time = Instant::now();

                    // Process kmers
                    let mut sig = protein_sig.clone();
                    index.process_kmers(TEST_PROTEIN, &mut sig).unwrap();

                    // Record end time
                    let end_time = Instant::now();
                    let _cpu_time = end_time.duration_since(start_time);

                    // Measure basic memory usage (stack size)
                    let _memory_used = mem::size_of_val(&sig);
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

    for moltype in MOLTYPES {
        for ksize in KSIZES {
            let temp_dir = tempdir().unwrap();
            let db_path = temp_dir.path().join(format!("db_{}_{}", moltype, ksize));

            // Create index with efficient storage disabled (store_raw_sequences = false)
            let index = ProteomeIndex::new(db_path.clone(), ksize, 1, moltype, false).unwrap();

            c.bench_function(&format!("process_fasta_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    // Record start time for CPU measurement
                    let start_time = Instant::now();

                    // Process the FASTA file
                    index.process_fasta(&fasta_path, 0, 1000).unwrap();

                    // Record end time
                    let end_time = Instant::now();
                    let cpu_time = end_time.duration_since(start_time);

                    // Calculate database file size
                    let mut total_size = 0u64;
                    if db_path.exists() {
                        if db_path.is_dir() {
                            // Sum up all files in the RocksDB directory
                            for entry in fs::read_dir(&db_path).unwrap().flatten() {
                                if let Ok(metadata) = entry.metadata() {
                                    total_size += metadata.len();
                                }
                            }
                        } else {
                            // Single file
                            if let Ok(metadata) = fs::metadata(&db_path) {
                                total_size = metadata.len();
                            }
                        }
                    }

                    // Print metrics (these will be captured by criterion)
                    println!("CPU time: {:?}", cpu_time);
                    println!("Database size: {} bytes", total_size);
                })
            });
        }
    }
}

fn benchmark_process_fasta_with_efficient_storage(c: &mut Criterion) {
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

    for moltype in MOLTYPES {
        for ksize in KSIZES {
            let temp_dir = tempdir().unwrap();
            let db_path = temp_dir.path().join(format!("db_efficient_{}_{}", moltype, ksize));

            // Create index with efficient storage enabled (store_raw_sequences = true)
            let index = ProteomeIndex::new(db_path.clone(), ksize, 1, moltype, true).unwrap();

            c.bench_function(&format!("process_fasta_efficient_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    // Record start time for CPU measurement
                    let start_time = Instant::now();

                    // Process the FASTA file
                    index.process_fasta(&fasta_path, 0, 1000).unwrap();

                    // Record end time
                    let end_time = Instant::now();
                    let cpu_time = end_time.duration_since(start_time);

                    // Calculate database file size
                    let mut total_size = 0u64;
                    if db_path.exists() {
                        if db_path.is_dir() {
                            // Sum up all files in the RocksDB directory
                            for entry in fs::read_dir(&db_path).unwrap().flatten() {
                                if let Ok(metadata) = entry.metadata() {
                                    total_size += metadata.len();
                                }
                            }
                        } else {
                            // Single file
                            if let Ok(metadata) = fs::metadata(&db_path) {
                                total_size = metadata.len();
                            }
                        }
                    }

                    // Print metrics (these will be captured by criterion)
                    println!("CPU time: {:?}", cpu_time);
                    println!("Database size: {} bytes", total_size);
                })
            });
        }
    }
}

/// Benchmark search throughput: queries/sec for search_one() with a pre-built index.
///
/// Uses bcl2_first25 as the target database and ced9 as the query protein.
/// The index is built once outside of timing; each iteration measures pure search cost.
///
/// Also benchmarks batch sizes to help tune the --batch-size CLI default.
fn benchmark_search_throughput(c: &mut Criterion) {
    use kmerseek::search::ProteinSearcher;
    use kmerseek::sketch::ProteinSketch;

    let target_fasta =
        "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";
    let query_fasta = "tests/testdata/fasta/ced9.fasta";

    // Read query sequence once
    let query_seq = {
        let mut reader = needletail::parse_fastx_file(query_fasta).unwrap();
        let record = reader.next().unwrap().unwrap();
        (
            std::str::from_utf8(record.id()).unwrap().to_string(),
            std::str::from_utf8(&record.seq()).unwrap().to_uppercase(),
        )
    };

    for moltype in ["hp", "protein", "dayhoff"] {
        for ksize in [10u32, 12] {
            let temp_dir = tempdir().unwrap();
            let db_path = temp_dir.path().join(format!("bench_search_{}_{}", moltype, ksize));

            // Build index once (not timed); drop it to release the RocksDB lock
            {
                let index = ProteomeIndex::new(db_path.clone(), ksize, 1, moltype, true).unwrap();
                index.process_fasta(target_fasta, 0, 1000).unwrap();
                index.save_state().unwrap();
            }

            // benchmark: ProteinSearcher::load() startup time
            {
                let bench_name = format!("searcher_load_{moltype}_k{ksize}");
                let db = db_path.clone();
                c.bench_function(&bench_name, |b| {
                    b.iter(|| ProteinSearcher::load(&db).unwrap())
                });
            }

            // Load searcher once for search benchmarks
            let searcher = ProteinSearcher::load(&db_path).unwrap();
            let mut query_sig =
                ProteinSketch::new(&query_seq.0, ksize, 1, moltype).unwrap();
            query_sig.add_protein(&query_seq.1, true).unwrap();

            // benchmark: single query search_one()
            {
                let bench_name = format!("search_one_{moltype}_k{ksize}");
                c.bench_function(&bench_name, |b| {
                    b.iter(|| searcher.search_one(&query_sig))
                });
            }

            // benchmark: batch search at several batch sizes (simulated by repeating query)
            for n_queries in [10usize, 100, 500] {
                let queries: Vec<ProteinSketch> = vec![query_sig.clone(); n_queries];
                let bench_name =
                    format!("search_batch{n_queries}_{moltype}_k{ksize}");
                c.bench_function(&bench_name, |b| {
                    b.iter(|| {
                        queries
                            .iter()
                            .map(|q| searcher.search_one(q))
                            .collect::<Vec<_>>()
                    })
                });
            }
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
    benchmark_process_fasta_with_efficient_storage,
    benchmark_search_throughput,
);
criterion_main!(benches);
