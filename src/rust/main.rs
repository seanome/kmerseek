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

        /// Scaled factor for minhash (1 = capture all k-mers)
        #[arg(short, long, default_value = "1")]
        scaled: u32,

        /// Protein encoding method
        #[arg(short, long, default_value = "protein")]
        encoding: ProteinEncoding,

        /// Progress notification interval (number of sequences between progress reports)
        #[arg(short, long, default_value = "10000")]
        progress_interval: u32,
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

        /// Scaled factor (must match the database; if not provided, will use database value)
        #[arg(short, long)]
        scaled: Option<u32>,

        /// Protein encoding method (must match the database)
        #[arg(short, long, default_value = "protein")]
        encoding: ProteinEncoding,

        /// Minimum containment threshold (0.0 = show all matches)
        #[arg(long, default_value = "0.0")]
        threshold: f64,

        /// Whether to output detailed match info to stderr (always extracts k-mers)
        #[arg(long, default_value = "false")]
        verbose: bool,

        /// Whether to treat query as a pre-indexed database instead of FASTA file
        #[arg(long, default_value = "false")]
        query_is_index: bool,
    },
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq)]
enum ProteinEncoding {
    /// Raw protein encoding (20 amino acids)
    Protein,
    /// Dayhoff encoding (6 groups)
    Dayhoff,
    /// HP encoding (hydrophobic/polar)
    Hp,
}

impl From<ProteinEncoding> for &'static str {
    fn from(encoding: ProteinEncoding) -> Self {
        match encoding {
            ProteinEncoding::Protein => "protein",
            ProteinEncoding::Dayhoff => "dayhoff",
            ProteinEncoding::Hp => "hp",
        }
    }
}

fn main() -> IndexResult<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Index { input, output, ksize, scaled, encoding, progress_interval } => {
            eprintln!("Indexing FASTA file: {}", input.display());

            // Determine output path
            let output_path = if let Some(output) = output {
                eprintln!("Output database: {}", output.display());
                output
            } else {
                // Auto-generate filename based on input file
                let base_name =
                    input.file_name().and_then(|name| name.to_str()).unwrap_or("unknown");

                // Create a temporary index to generate the filename
                let temp_index = ProteomeIndex::new_with_auto_filename(
                    &input,
                    ksize,
                    scaled,
                    encoding.into(),
                    true, // Always store raw sequences
                )?;

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
            eprintln!("Encoding: {:?}", encoding);
            eprintln!("Progress interval: {}", progress_interval);
            eprintln!("-------\n");

            // Create the index
            let index = ProteomeIndex::new(
                &output_path,
                ksize,
                scaled,
                encoding.into(),
                true, // Always store raw sequences
            )?;

            // Process the FASTA file
            eprintln!("Processing FASTA file...");
            index.process_fasta(&input, progress_interval, 1000)?;

            // Enable compactions for better read performance
            eprintln!("Optimizing database for read operations...");
            index.enable_compactions()?;

            // Save the index state for loading
            index.save_state()?;

            eprintln!("Indexing completed successfully!");
            eprintln!("Database saved to: {}", output_path.display());
        }
        Commands::Search {
            query,
            target,
            output,
            ksize,
            scaled,
            encoding,
            threshold,
            verbose,
            query_is_index,
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
                scaled,
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
            eprintln!("  Verbose output: {}", verbose);
            eprintln!("  Query is pre-indexed: {}\n---", query_is_index);

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
            let searcher = ProteinSearcher::load(&target)?;

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
                searcher.search_all_vs_all()?
            } else if query_is_index {
                // Load pre-indexed query database
                eprintln!("Loading pre-indexed query database...");
                let query_index = ProteomeIndex::load(&query)?;
                let query_signatures: Vec<_> = query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect();

                if query_signatures.is_empty() {
                    eprintln!("No query signatures found!");
                    return Ok(());
                }

                eprintln!("Found {} query signatures", query_signatures.len());
                eprintln!("Performing comprehensive search...");
                searcher.search(&query_signatures)?
            } else {
                // Stream queries from FASTA one at a time to minimize memory usage
                // WHY: Loading all query signatures into memory can require 50+ GB for large
                // proteomes. Instead, we process one sequence at a time: create signature,
                // search against target, collect results, discard signature.
                eprintln!("Streaming query sequences from FASTA...");
                use needletail::parse_fastx_file;
                use kmerseek::sketch::ProteinSketch;

                let mut reader = parse_fastx_file(&query)
                    .map_err(|e| anyhow::anyhow!("Failed to parse query FASTA: {}", e))?;

                let mut all_results = Vec::new();
                let mut query_count = 0u64;

                let progress = indicatif::ProgressBar::new_spinner();
                progress.set_style(
                    indicatif::ProgressStyle::with_template(
                        "{spinner:.green} [{elapsed_precise}] {msg}",
                    )
                    .unwrap()
                    .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ "),
                );

                while let Some(record) = reader.next() {
                    let record = record.map_err(|e| anyhow::anyhow!("FASTA parse error: {}", e))?;
                    let sequence = std::str::from_utf8(&record.seq())
                        .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in sequence: {}", e))?
                        .to_uppercase();
                    let name = std::str::from_utf8(record.id())
                        .map_err(|e| anyhow::anyhow!("Invalid UTF-8 in name: {}", e))?;

                    // Create a single query signature with raw sequences stored
                    let mut query_sig = ProteinSketch::new(
                        name, final_ksize, final_scaled, final_encoding.into(),
                    )?;
                    query_sig.add_protein(&sequence, true)?;

                    // Search this one query against all targets
                    let results = searcher.search_one(&query_sig);
                    all_results.extend(results);

                    query_count += 1;
                    if query_count % 1000 == 0 {
                        progress.set_message(format!("Searched {} queries, {} matches so far", query_count, all_results.len()));
                        progress.tick();
                    }
                }

                progress.finish_with_message(format!("Searched {} queries, {} matches total", query_count, all_results.len()));

                // Sort by containment
                all_results.sort_by(|a, b| {
                    b.containment.partial_cmp(&a.containment).unwrap_or(std::cmp::Ordering::Equal)
                });

                all_results
            };

            // Filter results by threshold
            let filtered_results: Vec<_> = search_results
                .into_iter()
                .filter(|result| result.containment >= threshold)
                .collect();

            eprintln!("Found {} matches above threshold {}", filtered_results.len(), threshold);

            // Output CSV to stdout or file
            // WHY: We expand each SearchResult into multiple rows - one per matched region.
            // Each row contains all the SearchResult similarity metrics plus the matched region
            // information. We only output SearchResults that have matched regions - if there are
            // no matched regions, the SearchResult is skipped. This ensures every CSV row has
            // complete matched region information.
            use kmerseek::search::SearchResultCsv;
            if let Some(output_path) = output {
                eprintln!("Writing results to: {}", output_path.display());
                let mut writer = csv::Writer::from_path(output_path)?;

                for result in &filtered_results {
                    // Output one row per matched region
                    // WHY: Every CSV row must have matched region data. Each SearchResult produces
                    // multiple CSV rows (one per matched region), with all similarity metrics
                    // repeated for each region.
                    for region in &result.matched_regions {
                        let csv_row =
                            SearchResultCsv::from_result_and_region(result, region);
                        writer.serialize(&csv_row)?;
                    }
                }

                writer.flush()?;
            } else {
                // Output to stdout
                let mut writer = csv::Writer::from_writer(std::io::stdout());

                for result in &filtered_results {
                    // Output one row per matched region
                    // WHY: Every CSV row must have matched region data. Each SearchResult produces
                    // multiple CSV rows (one per matched region), with all similarity metrics
                    // repeated for each region.
                    for region in &result.matched_regions {
                        let csv_row =
                            SearchResultCsv::from_result_and_region(result, region);
                        writer.serialize(&csv_row)?;
                    }
                }

                writer.flush()?;
            }

            // Display summary statistics (TF-IDF and overlap probabilities are now included in results)
            // WHY: Summary statistics go to stderr so they don't interfere with CSV output to stdout.
            // This allows users to pipe CSV data to other tools while still seeing progress and summary info.
            eprintln!("\n=== Search Summary ===");
            eprintln!("Total matches found: {}", filtered_results.len());
            if !filtered_results.is_empty() {
                let avg_containment: f64 =
                    filtered_results.iter().map(|r| r.containment).sum::<f64>()
                        / filtered_results.len() as f64;
                let avg_tfidf: f64 = filtered_results.iter().map(|r| r.tfidf).sum::<f64>()
                    / filtered_results.len() as f64;
                let avg_database_kmer_freq: f64 =
                    filtered_results.iter().map(|r| r.average_database_kmer_frequency).sum::<f64>()
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
        _ => {
            return Err(kmerseek::errors::IndexError::ValidationError {
                message: format!(
                    "Unknown encoding in database: {}. Expected one of: protein, dayhoff, hp",
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
/// * `user_scaled` - User-provided scaled (None if not specified)
/// * `user_encoding` - User-provided encoding (may be default value)
/// * `detected_ksize` - Ksize detected from database
/// * `detected_scaled` - Scaled detected from database
/// * `detected_moltype` - Moltype detected from database (as string)
///
/// # Returns
/// Tuple of (final_ksize, final_scaled, final_encoding) or ValidationError if mismatch
fn validate_and_assign_parameters(
    user_ksize: Option<u32>,
    user_scaled: Option<u32>,
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

    // Validate and assign scaled: use detected if not provided, error if mismatch
    // WHY: Same logic as ksize - database parameters are authoritative.
    let final_scaled = match user_scaled {
        Some(scaled) if scaled != detected_scaled => {
            return Err(kmerseek::errors::IndexError::ValidationError {
                message: format!(
                    "Scaled factor mismatch: database has scaled={}, but you specified --scaled={}.\n\
                    The scaled factor must match the database. Remove --scaled to use the database value ({}).",
                    detected_scaled, scaled, detected_scaled
                ),
            });
        }
        Some(scaled) => {
            // User provided scaled and it matches - use it
            scaled
        }
        None => {
            // User didn't provide scaled - use detected value
            eprintln!(
                "Using detected scaled: {} (not specified, using database value)",
                detected_scaled
            );
            detected_scaled
        }
    };

    // Validate and assign encoding
    let final_encoding = assign_encoding(user_encoding, detected_moltype)?;

    Ok((final_ksize, final_scaled, final_encoding))
}
