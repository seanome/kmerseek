use clap::{Parser, Subcommand, ValueEnum};
use kmerseek::alphabets::HpAlphabet;
use kmerseek::errors::IndexResult;
use kmerseek::types::MolType;
use kmerseek::{search::ProteinSearcher, ProteomeIndex};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "kmerseek")]
#[command(about = "Efficient protein domain annotation search with reduced amino acid k-mers")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Index a FASTA file (supports all compression formats)
    Index {
        /// Input FASTA file path (supports gzip, bzip2, xz, zstd, and uncompressed)
        #[arg(short, long)]
        input: PathBuf,

        /// Output database path (optional - will auto-generate if not provided)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// K-mer size for indexing
        #[arg(short, long, default_value = "10")]
        ksize: u32,

        /// Reduced amino acid alphabet to index with
        #[arg(short = 'a', long, default_value = "protein20")]
        alphabet: ProteinAlphabet,

        /// Seed for hp_random_control2 (1-10). Produces alphabet hp_random_control2_N.
        #[arg(long)]
        random_seed: Option<u64>,

        /// Progress notification interval (number of sequences between progress reports)
        #[arg(short, long, default_value = "10000")]
        progress_interval: u32,

        /// Write the k-mer frequency spectrum to this CSV path for plotting across alphabets
        /// and k-sizes. Gzip-compressed when the path ends in .gz. Columns:
        /// moltype, ksize, occurrences, n_kmers.
        #[arg(long, value_name = "PATH")]
        kmer_stats_out: Option<PathBuf>,

        /// Skip persisting a searchable index -- compute and write --kmer-stats-out only,
        /// with no RocksDB writes at all. Requires --kmer-stats-out. Use this when you only
        /// want the frequency spectrum, not a database to search later: it avoids both the
        /// SearchCache's RocksDB single-value size limit (~4 GiB) and the filesystem I/O load
        /// of chunked signature storage, neither of which stats-only output needs.
        #[arg(long, requires = "kmer_stats_out")]
        stats_only: bool,

        /// Remove low-complexity (homopolymer) k-mers from the index: raw
        /// amino-acid runs (e.g. "AAAAA") for any encoding, plus all-h or all-p
        /// runs for HP-family encodings (hp, hp_lehninger, hp_thomas_dill, etc.).
        /// The setting is stored in the index and reused automatically at search
        /// time, so you do not repeat it when searching.
        #[arg(long, default_value = "false")]
        remove_low_complexity: bool,
    },
    /// Search query sequences against a protein database
    Search {
        /// Query FASTA file path
        #[arg(short, long)]
        query: PathBuf,

        /// Target database path
        #[arg(short, long)]
        target: PathBuf,

        /// Output CSV file path (optional - will output to stdout if not provided)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// K-mer size (must match the database; if not provided, will use database value)
        #[arg(short, long)]
        ksize: Option<u32>,

        /// Reduced amino acid alphabet (must match the database)
        #[arg(short = 'a', long, default_value = "protein20")]
        alphabet: ProteinAlphabet,

        /// Seed for hp_random_control2 (1-10). Must match the seed used during indexing.
        #[arg(long)]
        random_seed: Option<u64>,

        /// Minimum containment threshold (0.0 = show all matches)
        #[arg(long, default_value = "0.0")]
        threshold: f64,

        /// Minimum number of shared k-mers required to report a match
        #[arg(long, default_value = "2")]
        min_shared_kmers: usize,

        /// Maximum uncorrected whole-query Poisson p-value required to report a match.
        /// A match is reported if either this or --min-region-score passes.
        #[arg(long, default_value = "0.05")]
        max_query_pvalue: f64,

        /// Minimum region-scoped score required to report a match, applied to the
        /// best-scoring region. Bigger means more surprising: the score is -log10 of the
        /// region's Poisson tail probability, so a p-value of 0.05 is a score of about 1.3,
        /// and a p-value of 0.0007 is a score of about 3.16. This is a heuristic ranking
        /// cutoff, not a statistically calibrated significance threshold (see the region
        /// scoring notes in the docs). A match is reported if either this or
        /// --max-query-pvalue passes, so a strong sub-protein domain hit survives even when
        /// the whole-query p-value is unimpressive. Defaults to about 1.3 (p=0.05).
        #[arg(long)]
        min_region_score: Option<f64>,

        /// Deprecated: use --max-query-pvalue (whole protein) or --min-region-score (per
        /// matched region). Kept as an alias that applies whole-query filtering only.
        #[arg(long)]
        max_pvalue: Option<f64>,

        /// Remove low-complexity (homopolymer) k-mers from query sketches.
        /// Omit this to follow whatever the target index was built with, which is
        /// almost always what you want. Pass it (or `--remove-low-complexity
        /// false`) only to override deliberately; a value that disagrees with the
        /// index is reported as a warning, because the two sides must match for
        /// containment to be comparable.
        #[arg(long, value_name = "BOOL", num_args = 0..=1, default_missing_value = "true")]
        remove_low_complexity: Option<bool>,

        /// Whether to output detailed match info to stderr (always extracts k-mers)
        #[arg(long, default_value = "false")]
        verbose: bool,

        /// Whether to treat query as a pre-indexed database instead of FASTA file
        #[arg(long, default_value = "false")]
        query_is_index: bool,

        /// Number of queries to process per parallel batch.
        /// Larger values use more memory but improve CPU utilization on many-core machines.
        /// Set to 1 to process queries one at a time (maximum streaming, minimum memory).
        #[arg(long, default_value = "500")]
        batch_size: usize,
    },
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq)]
enum ProteinAlphabet {
    /// The full 20-letter amino acid alphabet: no reduction
    #[value(name = "protein20")]
    Protein,
    /// Dayhoff, 6 classes
    #[value(name = "dayhoff6")]
    Dayhoff,
    /// HP Lehninger, 2 classes. Also accepts the older spellings `hp_lehninger` and `hp`;
    /// note that indexes built with the old sourmash-backed `hp` must be rebuilt, since it
    /// hashed the same partition differently
    #[value(name = "hp_lehninger2", alias = "hp-lehninger2")]
    HpLehninger,
    /// HP Thomas-Dill 1996 (C=h, G=p, P=p)
    #[value(name = "hp_thomas_dill2", alias = "hp-thomas-dill2")]
    HpThomasDill,
    /// HP Kyte-Doolittle 1982 binarized at hydropathy > 0 (W=p, Y=p)
    #[value(name = "hp_kyte_doolittle2", alias = "hp-kyte-doolittle2")]
    HpKyteDoolittle,
    /// HP Thomas-Dill with C reassigned to polar (isolation variant)
    #[value(name = "hp_thomas_dill_no_c2", alias = "hp-thomas-dill-no-c2")]
    HpThomasDillNoC,
    /// HP Lehninger with C reassigned to hydrophobic (isolation variant)
    #[value(name = "hp_lehninger_c_nonpolar2", alias = "hp-lehninger-c-nonpolar2")]
    HpLehningerCNonpolar,
    /// HPC Lehninger 3-letter: hydrophobic/polar/cystine, C split into its own class
    #[value(name = "hp_lehninger_hpc3", alias = "hp-lehninger-hpc3")]
    HpLehningerHpc,
    /// HP Physical Biology of the Cell 1st ed (Phillips et al. 2008)
    #[value(name = "hp_pbotc_1st_ed2", alias = "hp-pbotc-1st-ed2")]
    HpPBotC1stEd,
    /// HP negative control, 2 classes: the h/p split is randomized, scrambling the
    /// hydrophobicity signal. Use --random-seed for independent replicates
    #[value(name = "hp_random_control2", alias = "hp-random-control2")]
    HpRandomControl,
    /// GBMR4, 4 classes (Solis & Rackovsky 2000; best recall in Peterson et al. 2009)
    #[value(name = "gbmr4")]
    ReducedGbmr4,
    /// WWMJ5, 5 classes (Wang & Wang 1999, Miyazawa-Jernigan contact potentials)
    #[value(name = "wwmj5")]
    ReducedWwmj5,
    /// GBMR7, 7 classes (Solis & Rackovsky 2000)
    #[value(name = "gbmr7")]
    ReducedGbmr7,
    /// SDM12, 12 classes (Prlic et al. 2000; best AUC in Peterson et al. 2009)
    #[value(name = "sdm12")]
    ReducedSdm12,
    /// MMSEQS12, 12 classes (Steinegger & Soding 2018)
    #[value(name = "mmseqs12")]
    ReducedMmseqs12,
    /// WASS14, 14 classes, hydrophobicity-clustered (Ieremie et al. 2024)
    #[value(name = "wass14")]
    ReducedWass14,
    /// HSDM17, 17 classes (Prlic et al. 2000; best precision in Peterson et al. 2009)
    #[value(name = "hsdm17")]
    ReducedHsdm17,
    /// UNIPROT18, 18 classes, learned by a protein language model (Ieremie et al. 2024)
    #[value(name = "uniprot18")]
    ReducedUniprot18,
}

impl From<ProteinAlphabet> for &'static str {
    fn from(encoding: ProteinAlphabet) -> Self {
        match encoding {
            ProteinAlphabet::Protein => "protein20",
            ProteinAlphabet::Dayhoff => "dayhoff6",
            ProteinAlphabet::HpLehninger => "hp_lehninger2",
            ProteinAlphabet::HpThomasDill => "hp_thomas_dill2",
            ProteinAlphabet::HpKyteDoolittle => "hp_kyte_doolittle2",
            ProteinAlphabet::HpThomasDillNoC => "hp_thomas_dill_no_c2",
            ProteinAlphabet::HpLehningerCNonpolar => "hp_lehninger_c_nonpolar2",
            ProteinAlphabet::HpLehningerHpc => "hp_lehninger_hpc3",
            ProteinAlphabet::HpPBotC1stEd => "hp_pbotc_1st_ed2",
            ProteinAlphabet::HpRandomControl => "hp_random_control2",
            ProteinAlphabet::ReducedGbmr4 => "gbmr4",
            ProteinAlphabet::ReducedWwmj5 => "wwmj5",
            ProteinAlphabet::ReducedGbmr7 => "gbmr7",
            ProteinAlphabet::ReducedSdm12 => "sdm12",
            ProteinAlphabet::ReducedMmseqs12 => "mmseqs12",
            ProteinAlphabet::ReducedWass14 => "wass14",
            ProteinAlphabet::ReducedHsdm17 => "hsdm17",
            ProteinAlphabet::ReducedUniprot18 => "uniprot18",
        }
    }
}

fn main() -> IndexResult<()> {
    let cli = Cli::parse();

    eprintln!("kmerseek {}", env!("CARGO_PKG_VERSION"));

    match cli.command {
        Commands::Index {
            input,
            output,
            ksize,
            alphabet,
            random_seed,
            progress_interval,
            kmer_stats_out,
            stats_only,
            remove_low_complexity,
        } => {
            eprintln!("Indexing FASTA file: {}", input.display());

            // Scaled factor is always 1 (captures all k-mers)
            let scaled: u32 = 1;

            // Resolve effective moltype: seeded shuffled control -> "hp_shuffled_control_N".
            let effective_moltype: String = match (alphabet, random_seed) {
                (ProteinAlphabet::HpRandomControl, Some(seed)) => {
                    assert!((1..=10).contains(&seed), "--random-seed must be 1-10, got {seed}");
                    HpAlphabet::Random(seed).to_moltype()
                }
                _ => {
                    let s: &'static str = alphabet.into();
                    s.to_string()
                }
            };

            // Determine output path
            let output_path = if let Some(output) = output {
                eprintln!("Output database: {}", output.display());
                output
            } else {
                // Auto-generate filename based on input file
                let base_name =
                    input.file_name().and_then(|name| name.to_str()).unwrap_or("unknown");

                // Create a temporary index to generate the filename
                let mut temp_index = ProteomeIndex::new_with_auto_filename(
                    &input,
                    ksize,
                    scaled,
                    &effective_moltype,
                    true, // Always store raw sequences
                )?;
                // Must match the real index, so the generated name carries the
                // suffix and can't collide with a build that kept these k-mers.
                temp_index.set_remove_low_complexity(remove_low_complexity);

                let generated_filename = temp_index.generate_filename(base_name);
                let output_path = input
                    .parent()
                    .unwrap_or_else(|| std::path::Path::new("."))
                    .join(generated_filename);

                eprintln!("Auto-generated output database: {}", output_path.display());
                output_path
            };

            eprintln!("\n-------\nK-mer size: {}", ksize);
            eprintln!("Scaled: {}", scaled);
            eprintln!("Alphabet: {}", effective_moltype);
            eprintln!("Progress interval: {}", progress_interval);
            eprintln!("Remove low-complexity k-mers: {}", remove_low_complexity);
            eprintln!("-------\n");

            // Create the index
            let mut index = ProteomeIndex::new(
                &output_path,
                ksize,
                scaled,
                &effective_moltype,
                true, // Always store raw sequences
            )?;
            index.set_remove_low_complexity(remove_low_complexity);

            // Process the FASTA file
            eprintln!("Processing FASTA file...");
            index.process_fasta(&input, progress_interval, 1000)?;

            // Report what removal actually did, so its effect is visible without
            // having to rebuild and diff two indexes. Applies regardless of
            // --stats-only -- removal happens during process_fasta either way.
            if remove_low_complexity {
                let (examined, skipped) = index.low_complexity_counts();
                let percent =
                    if examined == 0 { 0.0 } else { 100.0 * skipped as f64 / examined as f64 };
                eprintln!(
                    "Removed {} of {} k-mer windows as low-complexity ({:.2}%)",
                    skipped, examined, percent
                );
            }

            if stats_only {
                // --stats-only: skip persisting a searchable index entirely (see
                // save_kmer_stats_only's doc comment for why). kmer_stats_out is
                // guaranteed Some here -- clap's `requires = "kmer_stats_out"` enforces it.
                let kmer_stats_out =
                    kmer_stats_out.expect("clap requires kmer_stats_out with stats_only");
                index.save_kmer_stats_only(&kmer_stats_out)?;
                eprintln!("Stats-only run completed successfully (no index persisted).");
            } else {
                // Enable compactions for better read performance
                eprintln!("Optimizing database for read operations...");
                index.enable_compactions()?;

                // Save the index state for loading
                index.save_state_with_kmer_stats(kmer_stats_out.as_deref())?;

                eprintln!("Indexing completed successfully!");
                eprintln!("Database saved to: {}", output_path.display());
            }
        }
        Commands::Search {
            query,
            target,
            output,
            ksize,
            alphabet,
            random_seed: _,
            threshold,
            min_shared_kmers,
            max_query_pvalue,
            min_region_score,
            max_pvalue,
            remove_low_complexity: remove_low_complexity_arg,
            verbose,
            query_is_index,
            batch_size,
        } => {
            eprintln!("Searching query sequences against target database");
            eprintln!("Query: {}", query.display());
            eprintln!("Target: {}", target.display());

            // Autodetect parameters from the target database
            eprintln!("Autodetecting parameters from target database...");
            let (detected_ksize, detected_scaled, detected_moltype) =
                ProteomeIndex::get_index_parameters(&target)?;

            // Validate and assign all parameters
            // WHY: This method centralizes parameter validation logic, making the main search
            // command handler much easier to read. It validates that user-provided parameters
            // match the database, or uses detected values if not provided. This is idiomatic
            // Rust - we extract complex logic into well-named methods for clarity.
            let (final_ksize, final_scaled, final_alphabet) = validate_and_assign_parameters(
                ksize,
                alphabet,
                detected_ksize,
                detected_scaled,
                &detected_moltype,
            )?;

            eprintln!("\n---\nUsing parameters:");
            eprintln!("  K-mer size: {} (detected: {})", final_ksize, detected_ksize);
            eprintln!("  Scaled: {} (detected: {})", final_scaled, detected_scaled);
            eprintln!("  Alphabet: {:?} (detected: {})", final_alphabet, detected_moltype);
            // --max-pvalue predates region scoring, so honour it as whole-query filtering only:
            // a region floor of infinity can never be cleared (the check is a strict >),
            // leaving the query scope as the only decider, the same as before region scoring
            // existed.
            let (max_query_pvalue, min_region_score) = match max_pvalue {
                Some(deprecated) => {
                    if let Some(ignored) = min_region_score {
                        eprintln!(
                            "WARNING: --min-region-score {ignored} is ignored because \
                             --max-pvalue was also passed; the region scope is forced to \
                             infinity (never passes) to reproduce pre-region-scoring \
                             behaviour."
                        );
                    }
                    eprintln!(
                        "WARNING: --max-pvalue is deprecated; it now applies whole-query \
                         filtering only.\n         Use --max-query-pvalue {deprecated} for the \
                         same behaviour, or --min-region-score to\n         keep sub-protein \
                         domain hits whose whole-query p-value is unimpressive."
                    );
                    (deprecated, f64::INFINITY)
                }
                // -log10(0.05): the score-scale equivalent of the same 0.05 default this flag
                // used before the -log10 transform.
                None => (max_query_pvalue, min_region_score.unwrap_or(-0.05_f64.log10())),
            };

            eprintln!("  Threshold: {}", threshold);
            eprintln!("  Minimum shared k-mers: {}", min_shared_kmers);
            eprintln!("  Maximum query p-value: {}", max_query_pvalue);
            eprintln!("  Minimum region score: {}", min_region_score);
            eprintln!("  Verbose output: {}", verbose);
            eprintln!("  Query is pre-indexed: {}\n---", query_is_index);

            use kmerseek::search::SearchFilters;
            let filters =
                SearchFilters { threshold, min_shared_kmers, max_query_pvalue, min_region_score };

            // Check if query and target are the same database (all-vs-all search)
            // WHY: RocksDB doesn't allow the same database to be opened twice by the same process.
            // When doing an all-vs-all search (query == target), we need to reuse the same database
            // instance instead of opening it twice. This prevents "No locks available" errors.
            let is_all_vs_all = if query_is_index {
                // Compare paths using canonicalize to handle symlinks and relative paths
                let query_path = query.canonicalize().ok().unwrap_or_else(|| query.clone());
                let target_path = target.canonicalize().ok().unwrap_or_else(|| target.clone());
                query_path == target_path
            } else {
                false
            };

            // Load the target database
            eprintln!("Loading target database...");
            let mut searcher = ProteinSearcher::load(&target)?;

            // Build query sketches the same way the target index was built.
            // WHY: if the index dropped low-complexity k-mers but queries keep them,
            // those k-mers match nothing yet still count toward the query cardinality,
            // deflating containment (intersection / query_size) for exactly the queries
            // that contain low-complexity regions.
            let index_removed = searcher.index().remove_low_complexity();
            let remove_low_complexity = remove_low_complexity_arg.unwrap_or(index_removed);

            eprintln!(
                "  Index: low-complexity k-mers were {} when it was built",
                if index_removed { "REMOVED" } else { "KEPT" }
            );
            eprintln!(
                "  This search: low-complexity k-mers are {} from query sketches ({})",
                if remove_low_complexity { "REMOVED" } else { "KEPT" },
                if remove_low_complexity_arg.is_some() {
                    "--remove-low-complexity"
                } else {
                    "matching the index"
                }
            );
            if remove_low_complexity != index_removed {
                eprintln!(
                    "  WARNING: this disagrees with the index. Containment is \
                     intersection / query_size, so k-mers present on only one side \
                     still count toward the denominator and skew scores."
                );
            }
            eprintln!(
                "  Index built by kmerseek: {}",
                searcher.index().kmerseek_version().unwrap_or("unknown (pre-versioning index)")
            );

            // Perform search - use optimized all-vs-all method if query == target
            let search_results = if is_all_vs_all {
                // Use optimized all-vs-all search that avoids cloning signatures
                // WHY: When query == target, we can use a specialized method that works directly
                // with references from the index, avoiding expensive clones. This is much more
                // memory-efficient for large databases and automatically skips self-matches.
                eprintln!(
                    "Detected all-vs-all search (query == target), using optimized search method..."
                );
                eprintln!("Skipping self-matches (comparing MD5 sums)...");
                searcher.search_all_vs_all(&filters)?
            } else if query_is_index {
                // Load pre-indexed query database
                eprintln!("Loading pre-indexed query database...");
                let query_index = ProteomeIndex::load(&query)?;
                let query_signatures: Vec<_> = query_index
                    .get_signatures()
                    .iter()
                    .map(|entry| entry.value().clone())
                    .collect();

                if query_signatures.is_empty() {
                    eprintln!("No query signatures found!");
                    return Ok(());
                }

                eprintln!("Found {} query signatures", query_signatures.len());
                eprintln!("Performing comprehensive search...");
                searcher.search(&query_signatures, &filters)?
            } else {
                // Stream queries from FASTA, writing CSV results as we go
                eprintln!("Streaming query sequences from FASTA...");
                use kmerseek::search::SearchResultCsv;
                use kmerseek::sketch::ProteinSketch;
                use needletail::parse_fastx_file;

                // First pass: build query-proteome k-mer frequencies for joint_kmer_freq. Also
                // counts the total number of queries up front. That count is attached to each
                // result as run_n_queries (see SearchResult::run_n_queries) and is not used in
                // any correction.
                eprintln!("First pass: scanning query proteome for k-mer frequencies...");
                let mut total_queries: usize = 0;
                {
                    use std::collections::HashMap;
                    let mut qfreqs: HashMap<u64, usize> = HashMap::new();
                    let mut freq_reader = parse_fastx_file(&query)
                        .map_err(|e| anyhow::anyhow!("Failed to parse query FASTA: {}", e))?;
                    while let Some(record) = freq_reader.next() {
                        let record =
                            record.map_err(|e| anyhow::anyhow!("FASTA parse error: {}", e))?;
                        let sequence = std::str::from_utf8(&record.seq())
                            .map_err(|e| anyhow::anyhow!("Invalid UTF-8: {}", e))?
                            .to_uppercase();
                        let name = std::str::from_utf8(record.id())
                            .map_err(|e| anyhow::anyhow!("Invalid UTF-8: {}", e))?;
                        let mut sig = ProteinSketch::new(
                            name,
                            final_ksize,
                            final_scaled,
                            final_alphabet.into(),
                        )?;
                        sig.set_remove_low_complexity(remove_low_complexity);
                        sig.add_protein(&sequence, true)?;
                        for min in sig.signature().minhash.mins() {
                            *qfreqs.entry(min).or_insert(0) += 1;
                        }
                        total_queries += 1;
                    }
                    eprintln!(
                        "First pass complete: {} query sequences, {} unique k-mers",
                        total_queries,
                        qfreqs.len()
                    );
                    searcher.set_query_frequencies(qfreqs, total_queries);
                }

                let mut reader = parse_fastx_file(&query)
                    .map_err(|e| anyhow::anyhow!("Failed to parse query FASTA: {}", e))?;

                // Create CSV writer up front so we stream rows as they're found
                let mut csv_writer: Box<dyn std::io::Write> = if let Some(ref output_path) = output
                {
                    eprintln!("Streaming results to: {}", output_path.display());
                    Box::new(std::io::BufWriter::new(std::fs::File::create(output_path)?))
                } else {
                    Box::new(std::io::BufWriter::new(std::io::stdout()))
                };
                let mut writer = csv::Writer::from_writer(&mut csv_writer);

                let mut query_count = 0u64;
                let mut match_count = 0u64;
                let mut row_count = 0u64;

                let progress = indicatif::ProgressBar::new_spinner();
                progress.set_style(
                    indicatif::ProgressStyle::with_template(
                        "{spinner:.green} [{elapsed_precise}] {msg}",
                    )
                    .unwrap()
                    .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ "),
                );
                progress.enable_steady_tick(std::time::Duration::from_millis(250));

                // Process queries in parallel batches: `batch_size` queries searched in parallel
                // (par_iter), then results written to CSV sequentially.
                // Larger batches = better CPU utilization; smaller = lower peak memory.
                use rayon::prelude::*;
                let mut batch: Vec<ProteinSketch> = Vec::with_capacity(batch_size);

                // Helper closure: process one batch and write results to CSV
                let process_batch = |batch: &[ProteinSketch],
                                     writer: &mut csv::Writer<&mut Box<dyn std::io::Write>>,
                                     match_count: &mut u64,
                                     row_count: &mut u64|
                 -> anyhow::Result<()> {
                    // Search all queries in this batch in parallel. Results failing `filters`
                    // are never included (see SearchFilters), so no post-hoc filtering needed here.
                    let batch_results: Vec<Vec<kmerseek::search::SearchResult>> = batch
                        .par_iter()
                        .map(|q| searcher.search_one(q, &filters, total_queries))
                        .collect();

                    // Write results sequentially (preserves per-query ordering within batch)
                    for results in &batch_results {
                        for result in results {
                            *match_count += 1;
                            for region in &result.matched_regions {
                                let csv_row = SearchResultCsv::from_result_and_region(
                                    result,
                                    region,
                                    remove_low_complexity,
                                );
                                writer.serialize(&csv_row)?;
                                *row_count += 1;
                            }
                        }
                    }
                    writer.flush()?;
                    Ok(())
                };

                while let Some(record) = reader.next() {
                    let record = record.map_err(|e| anyhow::anyhow!("FASTA parse error: {}", e))?;
                    let sequence = std::str::from_utf8(&record.seq())
                        .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence: {}", e))?
                        .to_uppercase();
                    let name = std::str::from_utf8(record.id())
                        .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in name: {}", e))?;

                    let mut query_sig =
                        ProteinSketch::new(name, final_ksize, final_scaled, final_alphabet.into())?;
                    query_sig.set_remove_low_complexity(remove_low_complexity);
                    query_sig.add_protein(&sequence, true)?;
                    batch.push(query_sig);

                    if batch.len() >= batch_size {
                        process_batch(&batch, &mut writer, &mut match_count, &mut row_count)?;
                        query_count += batch.len() as u64;
                        batch.clear();
                        progress.set_message(format!(
                            "{} queries | {} matches | {} rows written | {:.1} queries/sec",
                            query_count,
                            match_count,
                            row_count,
                            query_count as f64 / progress.elapsed().as_secs_f64(),
                        ));
                    }
                }

                // Process final partial batch
                if !batch.is_empty() {
                    process_batch(&batch, &mut writer, &mut match_count, &mut row_count)?;
                    query_count += batch.len() as u64;
                }

                writer.flush()?;
                drop(writer);
                drop(csv_writer);

                progress.finish_with_message(format!(
                    "Done! {} queries | {} matches | {} rows | {:.1} queries/sec",
                    query_count,
                    match_count,
                    row_count,
                    query_count as f64 / progress.elapsed().as_secs_f64(),
                ));

                eprintln!("\n=== Search Summary ===");
                eprintln!("Total queries: {}", query_count);
                eprintln!("Total matches: {}", match_count);
                eprintln!("Total CSV rows: {}", row_count);
                return Ok(());
            };

            // search() / search_all_vs_all() already applied `filters` internally, so
            // search_results only contains matches that passed threshold/min_shared_kmers and
            // cleared the query p-value or the region score.
            let filtered_results = search_results;

            eprintln!(
                "Found {} matches above threshold {} with at least {} shared k-mers and \
                 query p-value < {} or region score > {}",
                filtered_results.len(),
                threshold,
                min_shared_kmers,
                max_query_pvalue,
                min_region_score
            );

            use kmerseek::search::SearchResultCsv;
            if let Some(output_path) = output {
                eprintln!("Writing results to: {}", output_path.display());
                let mut writer = csv::Writer::from_path(output_path)?;

                for result in &filtered_results {
                    for region in &result.matched_regions {
                        let csv_row = SearchResultCsv::from_result_and_region(
                            result,
                            region,
                            remove_low_complexity,
                        );
                        writer.serialize(&csv_row)?;
                    }
                }

                writer.flush()?;
            } else {
                let mut writer = csv::Writer::from_writer(std::io::stdout());

                for result in &filtered_results {
                    for region in &result.matched_regions {
                        let csv_row = SearchResultCsv::from_result_and_region(
                            result,
                            region,
                            remove_low_complexity,
                        );
                        writer.serialize(&csv_row)?;
                    }
                }

                writer.flush()?;
            }

            eprintln!("\n=== Search Summary ===");
            eprintln!("Total matches found: {}", filtered_results.len());
            if !filtered_results.is_empty() {
                let avg_containment: f64 =
                    filtered_results.iter().map(|r| r.containment).sum::<f64>()
                        / filtered_results.len() as f64;
                let avg_tfidf: f64 = filtered_results.iter().map(|r| r.query_tfidf).sum::<f64>()
                    / filtered_results.len() as f64;
                let avg_database_kmer_freq: f64 =
                    filtered_results.iter().map(|r| r.mean_matched_kmer_freq).sum::<f64>()
                        / filtered_results.len() as f64;

                eprintln!("Average containment: {:.6}", avg_containment);
                eprintln!("Average TF-IDF: {:.6}", avg_tfidf);
                eprintln!("Average database k-mer frequency: {:.6}", avg_database_kmer_freq);
            }
        }
    }

    Ok(())
}

fn assign_encoding(
    alphabet: ProteinAlphabet,
    detected_moltype: &str,
) -> kmerseek::errors::IndexResult<ProteinAlphabet> {
    // Convert detected moltype string to enum
    // WHY: We need to compare the user-provided encoding with the detected encoding.
    // The detected encoding comes from the database as a string, so we convert it to
    // the enum type for comparison.
    // Indexes built before class counts were added to the HP names store the old
    // "hp_<name>" moltype. HpAlphabet::from_moltype() still parses those, so normalizing
    // here lets the match below deal only in current names.
    let canonical = MolType::new(detected_moltype)
        .map(|moltype| moltype.get().to_string())
        .unwrap_or_else(|_| detected_moltype.to_string());

    let detected_alphabet = match canonical.as_str() {
        "protein20" => ProteinAlphabet::Protein,
        "dayhoff6" => ProteinAlphabet::Dayhoff,
        "hp_lehninger2" => ProteinAlphabet::HpLehninger,
        "hp_thomas_dill2" => ProteinAlphabet::HpThomasDill,
        "hp_kyte_doolittle2" => ProteinAlphabet::HpKyteDoolittle,
        "hp_thomas_dill_no_c2" => ProteinAlphabet::HpThomasDillNoC,
        "hp_lehninger_c_nonpolar2" => ProteinAlphabet::HpLehningerCNonpolar,
        "hp_lehninger_hpc3" => ProteinAlphabet::HpLehningerHpc,
        "hp_pbotc_1st_ed2" => ProteinAlphabet::HpPBotC1stEd,
        "hp_random_control2" => ProteinAlphabet::HpRandomControl,
        // Seeded shuffled controls carry the seed in the moltype; map them back to
        // HpShuffledControl so the encoding path picks them up via
        // HpAlphabet::from_moltype(), which parses the numeric suffix.
        s if s.starts_with("hp_random_control2_") => ProteinAlphabet::HpRandomControl,
        "gbmr4" => ProteinAlphabet::ReducedGbmr4,
        "wwmj5" => ProteinAlphabet::ReducedWwmj5,
        "gbmr7" => ProteinAlphabet::ReducedGbmr7,
        "sdm12" => ProteinAlphabet::ReducedSdm12,
        "mmseqs12" => ProteinAlphabet::ReducedMmseqs12,
        "wass14" => ProteinAlphabet::ReducedWass14,
        "hsdm17" => ProteinAlphabet::ReducedHsdm17,
        "uniprot18" => ProteinAlphabet::ReducedUniprot18,
        _ => {
            return Err(kmerseek::errors::IndexError::ValidationError {
                message: format!(
                    "Unknown alphabet in database: {}. Expected one of: protein20, dayhoff6, \
                     hp_lehninger2, hp_thomas_dill2, hp_kyte_doolittle2, hp_thomas_dill_no_c2, \
                     hp_lehninger_c_nonpolar2, hp_lehninger_hpc3, hp_pbotc_1st_ed2, \
                     hp_random_control2 (or hp_random_control2_N for seeded variants), \
                     gbmr4, wwmj5, gbmr7, sdm12, mmseqs12, wass14, hsdm17, uniprot18",
                    detected_moltype
                ),
            });
        }
    };

    // Validate encoding: if user provided encoding doesn't match database, error
    // WHY: The database encoding is authoritative. If the user explicitly provides
    // an encoding that doesn't match, that's an error. This prevents silent failures
    // where searches would produce incorrect results. However, since encoding has
    // a default value, we can't distinguish "user specified" from "using default",
    // so we only error if it's clearly wrong (not the default and doesn't match).
    // In practice, users should not specify --encoding and let it autodetect.
    if alphabet != detected_alphabet && alphabet != ProteinAlphabet::Protein {
        // User explicitly provided a non-default encoding that doesn't match
        return Err(kmerseek::errors::IndexError::ValidationError {
            message: format!(
                "Alphabet mismatch: database was built with {}, but you specified \
                 --alphabet={:?}.\n\
                 The alphabet must match the database. Remove --alphabet to use the database \
                 value ({:?}).",
                detected_moltype, alphabet, detected_alphabet
            ),
        });
    }

    Ok(detected_alphabet)
}

/// Validate and assign search parameters from user input and database detection
///
/// WHY: This function centralizes the parameter validation logic, making the main search
/// command handler easier to read. It validates that user-provided parameters match the
/// database (which is authoritative), or uses detected values if not provided. This follows
/// idiomatic Rust patterns: extract complex logic into well-named functions, validate
/// preconditions, and provide clear error messages.
///
/// # Arguments
/// * `user_ksize` - User-provided ksize (None if not specified)
/// * `user_encoding` - User-provided encoding (may be default value)
/// * `detected_ksize` - Ksize detected from database
/// * `detected_scaled` - Scaled detected from database
/// * `detected_moltype` - Moltype detected from database (as string)
///
/// # Returns
/// Tuple of (final_ksize, final_scaled, final_alphabet) or ValidationError if mismatch
fn validate_and_assign_parameters(
    user_ksize: Option<u32>,
    user_encoding: ProteinAlphabet,
    detected_ksize: u32,
    detected_scaled: u32,
    detected_moltype: &str,
) -> IndexResult<(u32, u32, ProteinAlphabet)> {
    // Validate and assign ksize: use detected if not provided, error if mismatch
    // WHY: The database parameters are authoritative. If the user explicitly provides
    // a ksize that doesn't match, that's an error (they're trying to search with wrong
    // parameters). If they don't provide ksize, we use the detected value. This is
    // idiomatic Rust - we validate preconditions and fail fast with clear error messages.
    let final_ksize = match user_ksize {
        Some(ksize) if ksize != detected_ksize => {
            return Err(kmerseek::errors::IndexError::ValidationError {
                message: format!(
                    "K-mer size mismatch: database has ksize={}, but you specified --ksize={}.\n\
                    The ksize must match the database. Remove --ksize to use the database value ({}).",
                    detected_ksize, ksize, detected_ksize
                ),
            });
        }
        Some(ksize) => {
            // User provided ksize and it matches - use it (though it's the same as detected)
            ksize
        }
        None => {
            // User didn't provide ksize - use detected value
            eprintln!(
                "Using detected ksize: {} (not specified, using database value)",
                detected_ksize
            );
            detected_ksize
        }
    };

    // Scaled is always 1 (captures all k-mers). The database is authoritative, so
    // error out if it was built with a different scaled factor.
    if detected_scaled != 1 {
        return Err(kmerseek::errors::IndexError::ValidationError {
            message: format!(
                "Scaled factor mismatch: database has scaled={}, but kmerseek only supports scaled=1.",
                detected_scaled
            ),
        });
    }
    let final_scaled = 1;

    // Validate and assign encoding
    let final_alphabet = assign_encoding(user_encoding, detected_moltype)?;

    Ok((final_ksize, final_scaled, final_alphabet))
}
