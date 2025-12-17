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

        /// K-mer size (must match the database)
        #[arg(short, long, default_value = "10")]
        ksize: u32,

        /// Scaled factor (must match the database)
        #[arg(short, long, default_value = "1")]
        scaled: u32,

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
            println!("Indexing FASTA file: {}", input.display());

            // Determine output path
            let output_path = if let Some(output) = output {
                println!("Output database: {}", output.display());
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

                println!("Auto-generated output database: {}", output_path.display());
                output_path
            };

            println!("\n-------\nK-mer size: {}", ksize);
            println!("Scaled: {}", scaled);
            println!("Encoding: {:?}", encoding);
            println!("Progress interval: {}", progress_interval);
            println!("-------\n");

            // Create the index
            let index = ProteomeIndex::new(
                &output_path,
                ksize,
                scaled,
                encoding.into(),
                true, // Always store raw sequences
            )?;

            // Process the FASTA file
            println!("Processing FASTA file...");
            index.process_fasta(&input, progress_interval, 1000)?;

            // Enable compactions for better read performance
            println!("Optimizing database for read operations...");
            index.enable_compactions()?;

            // Save the index state for loading
            index.save_state()?;

            println!("Indexing completed successfully!");
            println!("Database saved to: {}", output_path.display());
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
            println!("Searching query sequences against target database");
            println!("Query: {}", query.display());
            println!("Target: {}", target.display());

            // Autodetect parameters from the target database
            println!("Autodetecting parameters from target database...");
            let (detected_ksize, detected_scaled, detected_moltype) =
                ProteomeIndex::get_index_parameters(&target)?;

            let final_ksize = assign_with_warning(ksize, detected_ksize, "ksize");
            let final_scaled = assign_with_warning(scaled, detected_scaled, "scaled");
            let final_encoding = assign_encoding(encoding, &detected_moltype);

            println!("\n---\nUsing parameters:");
            println!("  K-mer size: {} (detected: {})", final_ksize, detected_ksize);
            println!("  Scaled: {} (detected: {})", final_scaled, detected_scaled);
            println!("  Encoding: {:?} (detected: {})", final_encoding, detected_moltype);
            println!("  Threshold: {}", threshold);
            println!("  Verbose output: {}", verbose);
            println!("  Query is pre-indexed: {}\n---", query_is_index);

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
            println!("Loading target database...");
            let searcher = ProteinSearcher::load(&target)?;

            // Perform search - use optimized all-vs-all method if query == target
            let search_results = if is_all_vs_all {
                // Use optimized all-vs-all search that avoids cloning signatures
                // WHY: When query == target, we can use a specialized method that works directly
                // with references from the index, avoiding expensive clones. This is much more
                // memory-efficient for large databases and automatically skips self-matches.
                println!(
                    "Detected all-vs-all search (query == target), using optimized search method..."
                );
                println!("Skipping self-matches (comparing MD5 sums)...");
                searcher.search_all_vs_all()?
            } else {
                // Get query signatures using the detected parameters
                let query_signatures: Vec<_> = if query_is_index {
                    // Load pre-indexed query database
                    println!("Loading pre-indexed query database...");
                    let query_index = ProteomeIndex::load(&query)?;
                    query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect()
                } else {
                    // Process query sequences and create signatures using detected parameters
                    println!("Processing query sequences with detected parameters...");
                    // WHY: We must store raw sequences for query signatures so that matched regions
                    // can be found. The find_matched_regions function requires raw sequences to extract
                    // subsequences. Without stored sequences, matched_regions will be empty and those
                    // SearchResults won't be included in the CSV output.
                    let query_index = ProteomeIndex::new_with_auto_filename(
                        &query,
                        final_ksize,
                        final_scaled,
                        final_encoding.into(),
                        true, // Store raw sequences so matched regions can be found
                    )?;

                    query_index.process_fasta(&query, 1000, 1000)?;
                    query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect()
                };

                if query_signatures.is_empty() {
                    eprintln!("No query signatures found!");
                    return Ok(());
                }

                println!("Found {} query signatures", query_signatures.len());

                // Perform comprehensive search (includes TF-IDF and overlap probability calculations)
                println!("Performing comprehensive search...");
                searcher.search(&query_signatures)?
            };

            // Filter results by threshold
            let filtered_results: Vec<_> = search_results
                .into_iter()
                .filter(|result| result.containment >= threshold)
                .collect();

            println!("Found {} matches above threshold {}", filtered_results.len(), threshold);

            // Output CSV to stdout or file
            // WHY: We expand each SearchResult into multiple rows - one per matched region.
            // Each row contains all the SearchResult similarity metrics plus the matched region
            // information. We only output SearchResults that have matched regions - if there are
            // no matched regions, the SearchResult is skipped. This ensures every CSV row has
            // complete matched region information.
            use kmerseek::search::SearchResultWithRegionCsv;
            if let Some(output_path) = output {
                println!("Writing results to: {}", output_path.display());
                let mut writer = csv::Writer::from_path(output_path)?;

                for result in &filtered_results {
                    // Output one row per matched region
                    // WHY: Every CSV row must have matched region data. Each SearchResult produces
                    // multiple CSV rows (one per matched region), with all similarity metrics
                    // repeated for each region.
                    for region in &result.matched_regions {
                        let csv_row =
                            SearchResultWithRegionCsv::from_result_and_region(result, region);
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
                            SearchResultWithRegionCsv::from_result_and_region(result, region);
                        writer.serialize(&csv_row)?;
                    }
                }

                writer.flush()?;
            }

            // Display summary statistics (TF-IDF and overlap probabilities are now included in results)
            println!("\n=== Search Summary ===");
            println!("Total matches found: {}", filtered_results.len());
            if !filtered_results.is_empty() {
                let avg_containment: f64 =
                    filtered_results.iter().map(|r| r.containment).sum::<f64>()
                        / filtered_results.len() as f64;
                let avg_tfidf: f64 = filtered_results.iter().map(|r| r.tfidf).sum::<f64>()
                    / filtered_results.len() as f64;
                let avg_overlap_prob: f64 =
                    filtered_results.iter().map(|r| r.overlap_probability).sum::<f64>()
                        / filtered_results.len() as f64;

                println!("Average containment: {:.6}", avg_containment);
                println!("Average TF-IDF: {:.6}", avg_tfidf);
                println!("Average overlap probability: {:.6}", avg_overlap_prob);
            }
        }
    }

    Ok(())
}

fn assign_encoding(encoding: ProteinEncoding, detected_moltype: &str) -> ProteinEncoding {
    let final_encoding = if encoding != ProteinEncoding::Protein {
        // If user specified non-default encoding
        eprintln!(
            "Warning: Overriding detected encoding {detected_moltype} with user-specified {encoding:?}",
        );
        encoding
    } else {
        // Convert detected moltype string to enum
        match detected_moltype {
            "protein" => ProteinEncoding::Protein,
            "dayhoff" => ProteinEncoding::Dayhoff,
            "hp" => ProteinEncoding::Hp,
            _ => {
                eprintln!("Warning: Unknown detected encoding {detected_moltype}, using protein",);
                ProteinEncoding::Protein
            }
        }
    };
    final_encoding
}

fn assign_with_warning(value: u32, detected_value: u32, value_name: &str) -> u32 {
    // Use detected parameters, but allow user overrides
    let final_value = if value != detected_value {
        // If user specified non-default ksize
        println!(
            "Warning: Overriding detected {value_name} {detected_value} with user-specified {value}",
        );
        value
    } else {
        detected_value
    };
    final_value
}
