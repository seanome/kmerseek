use clap::{Parser, Subcommand, ValueEnum};
use kmerseek::{ProteomeIndex, SEED};
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
    /// Index a FASTA file (supports .gz compression)
    Index {
        /// Input FASTA file path (supports .gz compression)
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

        /// Whether to store raw protein sequences (increases storage size)
        #[arg(long, default_value = "false")]
        store_raw_sequences: bool,
    },
}

#[derive(ValueEnum, Clone, Copy, Debug)]
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

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Index {
            input,
            output,
            ksize,
            scaled,
            encoding,
            progress_interval,
            store_raw_sequences,
        } => {
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
                    SEED,
                    store_raw_sequences,
                )?;

                let generated_filename = temp_index.generate_filename(base_name);
                let output_path = input
                    .parent()
                    .unwrap_or_else(|| std::path::Path::new("."))
                    .join(generated_filename);

                println!("Auto-generated output database: {}", output_path.display());
                output_path
            };

            println!("K-mer size: {}", ksize);
            println!("Scaled: {}", scaled);
            println!("Encoding: {:?}", encoding);
            println!("Progress interval: {}", progress_interval);
            println!("Store raw sequences: {}", store_raw_sequences);

            // Create the index
            let index = ProteomeIndex::new(
                &output_path,
                ksize,
                scaled,
                encoding.into(),
                SEED,
                store_raw_sequences,
            )?;

            // Process the FASTA file
            println!("Processing FASTA file...");
            index.process_fasta(&input, progress_interval)?;

            println!("Indexing completed successfully!");
            println!("Database saved to: {}", output_path.display());
        }
    }

    Ok(())
}
