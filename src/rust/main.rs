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

        /// Output database path
        #[arg(short, long)]
        output: PathBuf,

        /// K-mer size for indexing
        #[arg(short, long, default_value = "10")]
        ksize: u32,

        /// Scaled factor for minhash (1 = capture all k-mers)
        #[arg(short, long, default_value = "1")]
        scaled: u32,

        /// Protein encoding method
        #[arg(short, long, default_value = "protein")]
        encoding: ProteinEncoding,
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
        Commands::Index { input, output, ksize, scaled, encoding } => {
            println!("Indexing FASTA file: {}", input.display());
            println!("Output database: {}", output.display());
            println!("K-mer size: {}", ksize);
            println!("Scaled: {}", scaled);
            println!("Encoding: {:?}", encoding);

            // Create the index
            let index = ProteomeIndex::new(&output, ksize, scaled, encoding.into(), SEED)?;

            // Process the FASTA file
            println!("Processing FASTA file...");
            index.process_fasta(&input)?;

            println!("Indexing completed successfully!");
            println!("Database saved to: {}", output.display());
        }
    }

    Ok(())
}
