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

            // Use detected parameters, but allow user overrides
            let final_ksize = if ksize != 10 {
                // If user specified non-default ksize
                println!(
                    "Warning: Overriding detected ksize {} with user-specified {}",
                    detected_ksize, ksize
                );
                ksize
            } else {
                detected_ksize
            };

            let final_scaled = if scaled != 1 {
                // If user specified non-default scaled
                println!(
                    "Warning: Overriding detected scaled {} with user-specified {}",
                    detected_scaled, scaled
                );
                scaled
            } else {
                detected_scaled
            };

            let final_encoding = if encoding != ProteinEncoding::Protein {
                // If user specified non-default encoding
                println!(
                    "Warning: Overriding detected encoding '{}' with user-specified {:?}",
                    detected_moltype, encoding
                );
                encoding
            } else {
                // Convert detected moltype string to enum
                match detected_moltype.as_str() {
                    "protein" => ProteinEncoding::Protein,
                    "dayhoff" => ProteinEncoding::Dayhoff,
                    "hp" => ProteinEncoding::Hp,
                    _ => {
                        println!(
                            "Warning: Unknown detected encoding '{}', using protein",
                            detected_moltype
                        );
                        ProteinEncoding::Protein
                    }
                }
            };

            println!("\n---\nUsing parameters:");
            println!("  K-mer size: {} (detected: {})", final_ksize, detected_ksize);
            println!("  Scaled: {} (detected: {})", final_scaled, detected_scaled);
            println!("  Encoding: {:?} (detected: {})", final_encoding, detected_moltype);
            println!("  Threshold: {}", threshold);
            println!("  Verbose output: {}", verbose);
            println!("  Query is pre-indexed: {}\n---", query_is_index);

            // Load the target database
            println!("Loading target database...");
            let searcher = ProteinSearcher::load(&target)?;

            // Get query signatures using the detected parameters
            let query_signatures: Vec<_> = if query_is_index {
                // Load pre-indexed query database
                println!("Loading pre-indexed query database...");
                let query_index = ProteomeIndex::load(&query)?;
                query_index.get_signatures().iter().map(|entry| entry.value().clone()).collect()
            } else {
                // Process query sequences and create signatures using detected parameters
                println!("Processing query sequences with detected parameters...");
                let query_index = ProteomeIndex::new_with_auto_filename(
                    &query,
                    final_ksize,
                    final_scaled,
                    final_encoding.into(),
                    false, // Don't store raw sequences for query
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
            let search_results = searcher.search(&query_signatures)?;

            // Filter results by threshold
            let filtered_results: Vec<_> = search_results
                .into_iter()
                .filter(|result| result.containment >= threshold)
                .collect();

            println!("Found {} matches above threshold {}", filtered_results.len(), threshold);

            // Output detailed k-mer information to stderr if verbose
            if verbose {
                let detailed_results = searcher.search_with_kmer_extraction(&query_signatures)?;
                for result in &detailed_results {
                    eprintln!("{}", result.to_print);
                }
            }

            // Output CSV to stdout or file
            if let Some(output_path) = output {
                println!("Writing results to: {}", output_path.display());
                let mut writer = csv::Writer::from_path(output_path)?;

                for result in &filtered_results {
                    writer.serialize(result)?;
                }

                writer.flush()?;
            } else {
                // Output to stdout
                let mut writer = csv::Writer::from_writer(std::io::stdout());

                for result in &filtered_results {
                    writer.serialize(result)?;
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
