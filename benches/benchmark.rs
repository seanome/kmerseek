use criterion::{criterion_group, criterion_main, Criterion};
use kmerseek::hash_functions::{encode_by_alphabet, encode_with_fn, get_encoding_fn_from_moltype};
use kmerseek::index::ProteomeIndex;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::mem;
use std::path::PathBuf;
use std::time::Instant;
use tempfile::tempdir;

// Test protein sequence
const TEST_PROTEIN: &str = "PLANTYANDANIMALGENQMESCOFFEE";

// Test protein sequences with different characteristics for benchmarking
const TEST_PROTEIN_WITH_AMBIGUOUS: &str = "PLANTYANDANIMALGENQMESCOFFEEBZJ";
const TEST_PROTEIN_WITH_SPECIAL: &str = "PLANTYANDANIMALGENQMESCOFFEEXUO";
const TEST_PROTEIN_WITH_STOP: &str = "PLANTYANDANIMALGENQMESCOFFEE*EXTRA";

const KSIZES: [u32; 3] = [5, 10, 20];
const MOLTYPES: [&str; 3] = ["protein20", "hp_lehninger2", "dayhoff6"];

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
            let encoding_fn =
                kmerseek::hash_functions::get_encoding_fn_from_moltype(moltype).unwrap();
            c.bench_function(&format!("proteome_index_encode_kmer_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    // Record start time for CPU measurement
                    let start_time = Instant::now();

                    // Encode kmer
                    let encoded = kmerseek::hash_functions::encode_with_fn(
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
                        encode_by_alphabet(&TEST_PROTEIN[..ksize as usize], moltype).unwrap();

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
            c.bench_function(&format!("process_protein_kmers_{}_{}", moltype, ksize), |b| {
                b.iter(|| {
                    // Record start time for CPU measurement
                    let start_time = Instant::now();

                    // A fresh sketch each iteration: add_protein both sketches the
                    // sequence and maps k-mer positions, which is the work measured.
                    let mut sig = kmerseek::sketch::ProteinSketch::new(
                        "test_protein",
                        ksize,
                        1, // scaled
                        moltype,
                    )
                    .unwrap();
                    sig.add_protein(TEST_PROTEIN, true).unwrap();

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
    use kmerseek::search::{ProteinSearcher, SearchFilters};
    use kmerseek::sketch::ProteinSketch;

    let target_fasta =
        "tests/testdata/fasta/uniprotkb_protein_name_Uncharacterized_2025_04_15.fasta.gz";
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

    for moltype in ["hp_lehninger2", "protein20", "dayhoff6"] {
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
                c.bench_function(&bench_name, |b| b.iter(|| ProteinSearcher::load(&db).unwrap()));
            }

            // Load searcher once for search benchmarks
            let searcher = ProteinSearcher::load(&db_path).unwrap();
            let mut query_sig = ProteinSketch::new(&query_seq.0, ksize, 1, moltype).unwrap();
            query_sig.add_protein(&query_seq.1, true).unwrap();

            // benchmark: single query search_one()
            {
                let bench_name = format!("search_one_{moltype}_k{ksize}");
                c.bench_function(&bench_name, |b| {
                    b.iter(|| searcher.search_one(&query_sig, &SearchFilters::default(), 1))
                });
            }

            // benchmark: batch search at several batch sizes (simulated by repeating query)
            for n_queries in [10usize, 100, 500] {
                let queries: Vec<ProteinSketch> = vec![query_sig.clone(); n_queries];
                let bench_name = format!("search_batch{n_queries}_{moltype}_k{ksize}");
                c.bench_function(&bench_name, |b| {
                    b.iter(|| {
                        queries
                            .iter()
                            .map(|q| searcher.search_one(q, &SearchFilters::default(), n_queries))
                            .collect::<Vec<_>>()
                    })
                });
            }
        }
    }
}

/// Benchmark indexing and searching the 2.8k uncharacterized protein file with hp encoding
/// at large k-mer sizes (15, 20, 30) to see index build time, load time, and search selectivity.
fn benchmark_index_hp_large_k(c: &mut Criterion) {
    use kmerseek::search::{ProteinSearcher, SearchFilters};
    use kmerseek::sketch::ProteinSketch;

    let target_fasta =
        "tests/testdata/fasta/uniprotkb_protein_name_Uncharacterized_2025_04_15.fasta.gz";
    let query_fasta = "tests/testdata/fasta/ced9.fasta";

    let query_seq = {
        let mut reader = needletail::parse_fastx_file(query_fasta).unwrap();
        let record = reader.next().unwrap().unwrap();
        (
            std::str::from_utf8(record.id()).unwrap().to_string(),
            std::str::from_utf8(&record.seq()).unwrap().to_uppercase(),
        )
    };

    let mut group = c.benchmark_group("index_hp_large_k");
    group.sample_size(10);

    for ksize in [15u32, 20, 30] {
        // Benchmark index build time
        {
            let bench_name = format!("index_hp_k{ksize}");
            group.bench_function(&bench_name, |b| {
                b.iter(|| {
                    let temp_dir = tempdir().unwrap();
                    let db_path = temp_dir.path().join(format!("bench_idx_hp_{}", ksize));
                    let index =
                        ProteomeIndex::new(db_path.clone(), ksize, 1, "hp_lehninger2", true)
                            .unwrap();
                    index.process_fasta(target_fasta, 0, 1000).unwrap();
                    index.save_state().unwrap();
                });
            });
        }

        // Build index once for load + search benchmarks
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join(format!("bench_search_hp_{}", ksize));
        {
            let index =
                ProteomeIndex::new(db_path.clone(), ksize, 1, "hp_lehninger2", true).unwrap();
            index.process_fasta(target_fasta, 0, 1000).unwrap();
            index.save_state().unwrap();
        }

        // Benchmark searcher load time (proxy for inverted index size)
        {
            let load_name = format!("load_hp_k{ksize}");
            let db = db_path.clone();
            group.bench_function(&load_name, |b| b.iter(|| ProteinSearcher::load(&db).unwrap()));
        }

        // Benchmark search_one (shows selectivity improvement from larger k)
        {
            let searcher = ProteinSearcher::load(&db_path).unwrap();
            let mut query_sig =
                ProteinSketch::new(&query_seq.0, ksize, 1, "hp_lehninger2").unwrap();
            query_sig.add_protein(&query_seq.1, true).unwrap();
            let search_name = format!("search_one_hp_k{ksize}");
            group.bench_function(&search_name, |b| {
                b.iter(|| searcher.search_one(&query_sig, &SearchFilters::default(), 1))
            });
        }
    }
    group.finish();
}

/// Compare three approaches for k-mer position storage during indexing and search.
///
/// Approach 1 (current/as-is):
///   Index: HashMap<u64, KmerInfo{encoded_kmer, HashMap<String, Vec<usize>>}>
///   Search: O(1) HashMap lookup of pre-computed positions
///
/// Approach 2 (raw sequence only):
///   Index: just minhash.add_protein() + store raw AA string; no position pre-computation
///   Search: O(L) re-scan raw sequence against intersection hashes at search time
///
/// Approach 3 (positions only):
///   Index: minhash.add_protein() + flat HashMap<u64, Vec<usize>>
///   Search: O(1) HashMap lookup (same as approach 1 but much smaller data)
///
/// Serialized sizes are printed to stderr once per moltype/ksize combination.
fn benchmark_kmer_storage_approaches(c: &mut Criterion) {
    use kmerseek::hash_functions::get_hash_function_from_moltype;
    use kmerseek::search::find_matched_regions;
    use kmerseek::sketch::{ProteinSketch, PROTEIN_TO_MINHASH_RATIO};
    use kmerseek::SEED;
    use sourmash::_hash_murmur;
    use sourmash::signature::SigsTrait;
    use sourmash::sketch::minhash::KmerMinHash;
    use std::collections::{HashMap, HashSet};

    let query_fasta = "tests/testdata/fasta/ced9.fasta";
    let target_fasta = "tests/testdata/fasta/bcl2.fasta";

    let read_first_seq = |path: &str| -> (String, String) {
        let mut reader = needletail::parse_fastx_file(path).unwrap();
        let record = reader.next().unwrap().unwrap();
        (
            std::str::from_utf8(record.id()).unwrap().to_string(),
            std::str::from_utf8(&record.seq()).unwrap().to_uppercase(),
        )
    };

    let (query_name, query_seq) = read_first_seq(query_fasta);
    let (target_name, target_seq) = read_first_seq(target_fasta);

    let mut group = c.benchmark_group("kmer_storage");

    for moltype in ["hp_lehninger2", "protein20"] {
        for ksize in [10u32, 12] {
            let encoding_fn = get_encoding_fn_from_moltype(moltype).unwrap();
            let hash_fn = get_hash_function_from_moltype(moltype).unwrap();
            let minhash_ksize = ksize * PROTEIN_TO_MINHASH_RATIO;
            let k = ksize as usize;

            // ---- SERIALIZED SIZE (computed once per config, printed to stderr) ----
            {
                // Approach 1: full ProteinSketchStore (current)
                let sketch1 = ProteinSketch::from_protein_sequence(
                    &query_name,
                    &query_seq,
                    ksize,
                    1,
                    moltype,
                )
                .unwrap();
                let store1 = sketch1.to_efficient_data(true);
                let bytes1 = bincode::serialize(&store1).unwrap().len();
                let n_mins = store1.mins.len();

                // Approach 2: (name, mins, abunds, raw_sequence)
                let mut mh = KmerMinHash::new(1, minhash_ksize, hash_fn.clone(), SEED, true, 0);
                mh.add_protein(query_seq.as_bytes()).unwrap();
                let mins = mh.mins().to_vec();
                let abunds = mh.abunds().map(|a| a.to_vec());
                let bytes2 =
                    bincode::serialize(&(&query_name, &mins, &abunds, &query_seq)).unwrap().len();

                // Approach 3: (name, mins, abunds, HashMap<u64, Vec<usize>>)
                let hashvals: HashSet<u64> = mins.iter().copied().collect();
                let mut positions: HashMap<u64, Vec<usize>> = HashMap::new();
                for i in 0..query_seq.len().saturating_sub(k - 1) {
                    if let Ok(encoded) = encode_with_fn(&query_seq[i..i + k], encoding_fn) {
                        let hash = _hash_murmur(encoded.as_bytes(), SEED);
                        if hashvals.contains(&hash) {
                            positions.entry(hash).or_default().push(i);
                        }
                    }
                }
                let bytes3 =
                    bincode::serialize(&(&query_name, &mins, &abunds, &positions)).unwrap().len();

                eprintln!(
                    "[kmer_storage {moltype} k={ksize}] protein={}aa mins={n_mins} | \
                     approach1(kmer_infos)={bytes1}B  approach2(raw_seq)={bytes2}B  approach3(positions)={bytes3}B",
                    query_seq.len()
                );
            }

            // ---- INDEXING BENCHMARKS ----

            // Approach 1: current — add_protein builds full nested kmer_infos
            {
                let name = query_name.clone();
                let seq = query_seq.clone();
                group.bench_function(format!("index_approach1_current_{moltype}_k{ksize}"), |b| {
                    b.iter(|| {
                        let mut sketch = ProteinSketch::new(&name, ksize, 1, moltype).unwrap();
                        sketch.add_protein(&seq, true).unwrap();
                        std::hint::black_box(sketch)
                    })
                });
            }

            // Approach 2: raw sequence only — just minhash + clone string
            {
                let seq = query_seq.clone();
                let hfn2 = hash_fn.clone();
                group.bench_function(format!("index_approach2_raw_seq_{moltype}_k{ksize}"), |b| {
                    b.iter(|| {
                        let mut mh =
                            KmerMinHash::new(1, minhash_ksize, hfn2.clone(), SEED, true, 0);
                        mh.add_protein(seq.as_bytes()).unwrap();
                        let mins = mh.mins().to_vec();
                        let abunds = mh.abunds().map(|a| a.to_vec());
                        let raw = seq.clone();
                        std::hint::black_box((mins, abunds, raw))
                    })
                });
            }

            // Approach 3: positions only — minhash + flat HashMap<u64, Vec<usize>>
            {
                let seq = query_seq.clone();
                let hfn3 = hash_fn;
                group.bench_function(
                    format!("index_approach3_positions_{moltype}_k{ksize}"),
                    |b| {
                        b.iter(|| {
                            let mut mh =
                                KmerMinHash::new(1, minhash_ksize, hfn3.clone(), SEED, true, 0);
                            mh.add_protein(seq.as_bytes()).unwrap();
                            let hset: HashSet<u64> = mh.mins().iter().copied().collect();
                            let mut pos: HashMap<u64, Vec<usize>> = HashMap::new();
                            for i in 0..seq.len().saturating_sub(k - 1) {
                                if let Ok(enc) = encode_with_fn(&seq[i..i + k], encoding_fn) {
                                    let h = _hash_murmur(enc.as_bytes(), SEED);
                                    if hset.contains(&h) {
                                        pos.entry(h).or_default().push(i);
                                    }
                                }
                            }
                            let mins = mh.mins().to_vec();
                            let abunds = mh.abunds().map(|a| a.to_vec());
                            std::hint::black_box((mins, abunds, pos))
                        })
                    },
                );
            }

            // ---- SEARCH: FIND MATCHED REGIONS ----

            let q_sketch =
                ProteinSketch::from_protein_sequence(&query_name, &query_seq, ksize, 1, moltype)
                    .unwrap();
            let t_sketch =
                ProteinSketch::from_protein_sequence(&target_name, &target_seq, ksize, 1, moltype)
                    .unwrap();
            let intersection = q_sketch.intersect(&t_sketch);

            if !intersection.is_empty() {
                eprintln!(
                    "[kmer_storage {moltype} k={ksize}] intersection={} hashes between ced9 and bcl2",
                    intersection.len()
                );

                // Approach 1/3: full find_matched_regions with pre-computed kmer_infos
                group.bench_function(
                    format!("search_approach1_precomputed_{moltype}_k{ksize}"),
                    |b| {
                        b.iter(|| {
                            std::hint::black_box(find_matched_regions(
                                &q_sketch,
                                &t_sketch,
                                &intersection,
                            ))
                        })
                    },
                );

                // Approach 2: re-scan both sequences against intersection hashes
                // This is the extra cost approach 2 adds per matched pair at search time.
                {
                    let q = query_seq.clone();
                    let t = target_seq.clone();
                    let isect = intersection.clone();
                    group.bench_function(
                        format!("search_approach2_rescan_{moltype}_k{ksize}"),
                        |b| {
                            b.iter(|| {
                                let mut q_pos: HashMap<u64, Vec<usize>> = HashMap::new();
                                for i in 0..q.len().saturating_sub(k - 1) {
                                    if let Ok(enc) = encode_with_fn(&q[i..i + k], encoding_fn) {
                                        let h = _hash_murmur(enc.as_bytes(), SEED);
                                        if isect.contains(&h) {
                                            q_pos.entry(h).or_default().push(i);
                                        }
                                    }
                                }
                                let mut t_pos: HashMap<u64, Vec<usize>> = HashMap::new();
                                for i in 0..t.len().saturating_sub(k - 1) {
                                    if let Ok(enc) = encode_with_fn(&t[i..i + k], encoding_fn) {
                                        let h = _hash_murmur(enc.as_bytes(), SEED);
                                        if isect.contains(&h) {
                                            t_pos.entry(h).or_default().push(i);
                                        }
                                    }
                                }
                                std::hint::black_box((q_pos, t_pos))
                            })
                        },
                    );
                }
            }
        }
    }

    group.finish();
}

/// Benchmark rebuild_combined_minhash() in isolation at different k-mer sizes.
///
/// Builds the index once (not timed), then repeatedly times only the rebuild step.
/// Tests k=10, 20, 30 to show that the O(N log N) sort+add scales well as the
/// number of unique k-mers grows with larger k.
fn benchmark_rebuild_combined_minhash(c: &mut Criterion) {
    let target_fasta =
        "tests/testdata/fasta/uniprotkb_protein_name_Uncharacterized_2025_04_15.fasta.gz";

    let mut group = c.benchmark_group("rebuild_combined_minhash");
    group.sample_size(20);

    for ksize in [10u32, 20, 30] {
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join(format!("bench_rebuild_k{ksize}"));

        // Build index once so signatures are in memory (not timed)
        let index = ProteomeIndex::new(db_path, ksize, 1, "hp_lehninger2", false).unwrap();
        index.process_fasta(target_fasta, 0, 1000).unwrap();

        let n_unique = index.combined_minhash_size();
        let bench_name = format!("rebuild_hp_k{ksize}_{n_unique}_unique_kmers");
        group.bench_function(&bench_name, |b| b.iter(|| index.rebuild_combined_minhash().unwrap()));
    }
    group.finish();
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
    benchmark_index_hp_large_k,
    benchmark_kmer_storage_approaches,
    benchmark_rebuild_combined_minhash,
);
criterion_main!(benches);
