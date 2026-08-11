use clap::{Parser, Subcommand, ValueEnum};
use kmerseek::errors::IndexResult;
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

        /// Protein encoding method
        #[arg(short, long, default_value = "protein")]
        encoding: ProteinEncoding,

        /// Seed for hp_shuffled_control (1-10). Produces moltype hp_shuffled_control_N.
        #[arg(long)]
        shuffled_seed: Option<u64>,

        /// Progress notification interval (number of sequences between progress reports)
        #[arg(short, long, default_value = "10000")]
        progress_interval: u32,

        /// Write the k-mer frequency spectrum to this CSV path for plotting across alphabets
        /// and k-sizes. Gzip-compressed when the path ends in .gz. Columns:
        /// moltype, ksize, occurrences, n_kmers.
        #[arg(long, value_name = "PATH")]
        kmer_stats_out: Option<PathBuf>,

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

        /// Protein encoding method (must match the database)
        #[arg(short, long, default_value = "protein")]
        encoding: ProteinEncoding,

        /// Seed for hp_shuffled_control (1-10). Must match the seed used during indexing.
        #[arg(long)]
        shuffled_seed: Option<u64>,

        /// Minimum containment threshold (0.0 = show all matches)
        #[arg(long, default_value = "0.0")]
        threshold: f64,

        /// Minimum number of shared k-mers required to report a match
        #[arg(long, default_value = "2")]
        min_shared_kmers: usize,

        /// Maximum uncorrected Poisson p-value required to report a match
        #[arg(long, default_value = "0.05")]
        max_pvalue: f64,

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
enum ProteinEncoding {
    /// Raw protein encoding (20 amino acids)
    Protein,
    /// Dayhoff encoding (6 groups)
    Dayhoff,
    /// HP encoding — sourmash built-in Lehninger classification (backward-compatible)
    Hp,
    /// HP Lehninger (explicit; identical hashes to hp)
    #[value(alias = "hp_lehninger")]
    HpLehninger,
    /// HP Thomas-Dill 1996 (C=h, G=p, P=p)
    #[value(alias = "hp_thomas_dill")]
    HpThomasDill,
    /// HP Kyte-Doolittle 1982 binarized at hydropathy > 0 (W=p, Y=p)
    #[value(alias = "hp_kyte_doolittle")]
    HpKyteDoolittle,
    /// HP Thomas-Dill with C reassigned to polar (isolation variant)
    #[value(alias = "hp_thomas_dill_no_c")]
    HpThomasDillNoC,
    /// HP Lehninger with C reassigned to hydrophobic (isolation variant)
    #[value(alias = "hp_lehninger_plus_c")]
    HpLehningerPlusC,
    /// HP Physical Biology of the Cell 1st ed (Phillips et al. 2008)
    #[value(name = "hp-pbotc-1st-ed", alias = "hp_pbotc_1st_ed")]
    HpPBotC1stEd,
    /// HP shuffled negative control (scrambled hydrophobicity signal)
    #[value(alias = "hp_shuffled_control")]
    HpShuffledControl,
}

impl From<ProteinEncoding> for &'static str {
    fn from(encoding: ProteinEncoding) -> Self {
        match encoding {
            ProteinEncoding::Protein => "protein",
            ProteinEncoding::Dayhoff => "dayhoff",
            ProteinEncoding::Hp => "hp",
            ProteinEncoding::HpLehninger => "hp_lehninger",
            ProteinEncoding::HpThomasDill => "hp_thomas_dill",
            ProteinEncoding::HpKyteDoolittle => "hp_kyte_doolittle",
            ProteinEncoding::HpThomasDillNoC => "hp_thomas_dill_no_c",
            ProteinEncoding::HpLehningerPlusC => "hp_lehninger_plus_c",
            ProteinEncoding::HpPBotC1stEd => "hp_pbotc_1st_ed",
            ProteinEncoding::HpShuffledControl => "hp_shuffled_control",
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
            encoding,
            shuffled_seed,
            progress_interval,
            kmer_stats_out,
            remove_low_complexity,
        } => {
            eprintln!("Indexing FASTA file: {}", input.display());

            // Scaled factor is always 1 (captures all k-mers)
            let scaled: u32 = 1;

            // Resolve effective moltype: seeded shuffled control -> "hp_shuffled_control_N".
            let effective_moltype: String = match (encoding, shuffled_seed) {
                (ProteinEncoding::HpShuffledControl, Some(seed)) => {
                    assert!((1..=10).contains(&seed), "--shuffled-seed must be 1-10, got {seed}");
                    format!("hp_shuffled_control_{seed}")
                }
                _ => {
                    let s: &'static str = encoding.into();
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
            eprintln!("Encoding: {}", effective_moltype);
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
            // having to rebuild and diff two indexes.
            if remove_low_complexity {
                let (examined, skipped) = index.low_complexity_counts();
                let percent =
                    if examined == 0 { 0.0 } else { 100.0 * skipped as f64 / examined as f64 };
                eprintln!(
                    "Removed {} of {} k-mer windows as low-complexity ({:.2}%)",
                    skipped, examined, percent
                );
            }

            // Enable compactions for better read performance
            eprintln!("Optimizing database for read operations...");
            index.enable_compactions()?;

            // Save the index state for loading
            index.save_state_with_kmer_stats(kmer_stats_out.as_deref())?;

            eprintln!("Indexing completed successfully!");
            eprintln!("Database saved to: {}", output_path.display());
        }
        Commands::Search {
            query,
            target,
            output,
            ksize,
            encoding,
            shuffled_seed: _,
            threshold,
            min_shared_kmers,
            max_pvalue,
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
            let (final_ksize, final_scaled, final_encoding) = validate_and_assign_parameters(
                ksize,
                encoding,
                detected_ksize,
                detected_scaled,
                &detected_moltype,
            )?;

            eprintln!("\n---\nUsing parameters:");
            eprintln!("  K-mer size: {} (detected: {})", final_ksize, detected_ksize);
            eprintln!("  Scaled: {} (detected: {})", final_scaled, detected_scaled);
            eprintln!("  Encoding: {:?} (detected: {})", final_encoding, detected_moltype);
            eprintln!("  Threshold: {}", threshold);
            eprintln!("  Minimum shared k-mers: {}", min_shared_kmers);
            eprintln!("  Maximum p-value: {}", max_pvalue);
            eprintln!("  Verbose output: {}", verbose);
            eprintln!("  Query is pre-indexed: {}\n---", query_is_index);

            use kmerseek::search::SearchFilters;
            let filters = SearchFilters { threshold, min_shared_kmers, max_pvalue };

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
            let remove_low_complexity = searcher.index().remove_low_complexity();
            eprintln!("  Remove low-complexity k-mers: {} (from index)", remove_low_complexity);
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

                // First pass: build query-proteome k-mer frequencies for joint_kmer_freq.
                eprintln!("First pass: scanning query proteome for k-mer frequencies...");
                {
                    use std::collections::HashMap;
                    let mut qfreqs: HashMap<u64, usize> = HashMap::new();
                    let mut total_queries: usize = 0;
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
                            final_encoding.into(),
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
                    let batch_results: Vec<Vec<kmerseek::search::SearchResult>> =
                        batch.par_iter().map(|q| searcher.search_one(q, &filters)).collect();

                    // Write results sequentially (preserves per-query ordering within batch)
                    for results in &batch_results {
                        for result in results {
                            *match_count += 1;
                            for region in &result.matched_regions {
                                let csv_row =
                                    SearchResultCsv::from_result_and_region(result, region);
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
                        ProteinSketch::new(name, final_ksize, final_scaled, final_encoding.into())?;
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
            // search_results only contains matches that passed threshold/min_shared_kmers/max_pvalue.
            let filtered_results = search_results;

            eprintln!(
                "Found {} matches above threshold {} with at least {} shared k-mers and p-value < {}",
                filtered_results.len(),
                threshold,
                min_shared_kmers,
                max_pvalue
            );

            use kmerseek::search::SearchResultCsv;
            if let Some(output_path) = output {
                eprintln!("Writing results to: {}", output_path.display());
                let mut writer = csv::Writer::from_path(output_path)?;

                for result in &filtered_results {
                    for region in &result.matched_regions {
                        let csv_row = SearchResultCsv::from_result_and_region(result, region);
                        writer.serialize(&csv_row)?;
                    }
                }

                writer.flush()?;
            } else {
                let mut writer = csv::Writer::from_writer(std::io::stdout());

                for result in &filtered_results {
                    for region in &result.matched_regions {
                        let csv_row = SearchResultCsv::from_result_and_region(result, region);
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
    encoding: ProteinEncoding,
    detected_moltype: &str,
) -> kmerseek::errors::IndexResult<ProteinEncoding> {
    // Convert detected moltype string to enum
    // WHY: We need to compare the user-provided encoding with the detected encoding.
    // The detected encoding comes from the database as a string, so we convert it to
    // the enum type for comparison.
    let detected_encoding = match detected_moltype {
        "protein" => ProteinEncoding::Protein,
        "dayhoff" => ProteinEncoding::Dayhoff,
        "hp" => ProteinEncoding::Hp,
        "hp_lehninger" => ProteinEncoding::HpLehninger,
        "hp_thomas_dill" => ProteinEncoding::HpThomasDill,
        "hp_kyte_doolittle" => ProteinEncoding::HpKyteDoolittle,
        "hp_thomas_dill_no_c" => ProteinEncoding::HpThomasDillNoC,
        "hp_lehninger_plus_c" => ProteinEncoding::HpLehningerPlusC,
        "hp_pbotc_1st_ed" => ProteinEncoding::HpPBotC1stEd,
        "hp_shuffled_control" => ProteinEncoding::HpShuffledControl,
        // Seeded shuffled controls (hp_shuffled_control_N) are stored with the seed
        // in the moltype; map them back to HpShuffledControl so the encoding path
        // picks them up via HpAlphabet::from_moltype() which parses the numeric suffix.
        s if s.starts_with("hp_shuffled_control_") => ProteinEncoding::HpShuffledControl,
        _ => {
            return Err(kmerseek::errors::IndexError::ValidationError {
                message: format!(
                    "Unknown encoding in database: {}. Expected one of: protein, dayhoff, hp, \
                     hp_lehninger, hp_thomas_dill, hp_kyte_doolittle, \
                     hp_thomas_dill_no_c, hp_lehninger_plus_c, hp_pbotc_1st_ed, \
                     hp_shuffled_control (or hp_shuffled_control_N for seeded variants)",
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
    if encoding != detected_encoding && encoding != ProteinEncoding::Protein {
        // User explicitly provided a non-default encoding that doesn't match
        return Err(kmerseek::errors::IndexError::ValidationError {
            message: format!(
                "Encoding mismatch: database has encoding={}, but you specified --encoding={:?}.\n\
                The encoding must match the database. Remove --encoding to use the database value ({:?}).",
                detected_moltype, encoding, detected_encoding
            ),
        });
    }

    Ok(detected_encoding)
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
/// Tuple of (final_ksize, final_scaled, final_encoding) or ValidationError if mismatch
fn validate_and_assign_parameters(
    user_ksize: Option<u32>,
    user_encoding: ProteinEncoding,
    detected_ksize: u32,
    detected_scaled: u32,
    detected_moltype: &str,
) -> IndexResult<(u32, u32, ProteinEncoding)> {
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
    let final_encoding = assign_encoding(user_encoding, detected_moltype)?;

    Ok((final_ksize, final_scaled, final_encoding))
}
