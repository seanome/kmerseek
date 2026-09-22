use clap::{Parser, Subcommand, ValueEnum};
use kmerseek::errors::{IndexError, IndexResult};
use kmerseek::karlin_altschul::{DecoyNull, KaCalibration, BIN_WIDTH};
use kmerseek::search::{KaCalibrationReport, KaCalibrationSettings, KaSource};
use kmerseek::types::{MolType, Scaled};
use kmerseek::{pair, search::ProteinSearcher, ProteomeIndex};
use std::path::{Path, PathBuf};

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

        /// Keep only k-mers whose hash falls in the lowest 1/scaled of the hash space
        /// (FracMinHash). 1 keeps every k-mer. Indexing memory and index size fall
        /// almost linearly with this value. A match is found from any one kept k-mer
        /// and reported at its full length; a match none of whose k-mers were kept is
        /// missed, which is likelier the shorter it is. Stored in the index; search
        /// reads it back. Maximum 10.
        #[arg(short, long, default_value = "1")]
        scaled: u32,

        /// Reduced amino acid alphabet to index with
        #[arg(short = 'a', long, default_value = "protein20")]
        alphabet: ProteinAlphabet,

        /// Progress notification interval (number of sequences between progress reports)
        #[arg(short, long, default_value = "10000")]
        progress_interval: u32,

        /// Write the k-mer frequency spectrum to this CSV path for plotting across alphabets
        /// and k-sizes. Gzip-compressed when the path ends in .gz. Columns:
        /// moltype, ksize, occurrences, n_kmers.
        #[arg(long, value_name = "PATH")]
        kmer_stats_out: Option<PathBuf>,

        /// Do not keep a searchable index: build it in a scratch directory that is removed
        /// on exit, and write only --kmer-stats-out. Requires --kmer-stats-out. Use this
        /// when you only want the frequency spectrum, not a database to search later.
        #[arg(long, requires = "kmer_stats_out")]
        stats_only: bool,

        /// Remove low-complexity (homopolymer) k-mers from the index: raw
        /// amino-acid runs (e.g. "AAAAA") for any encoding, plus all-h or all-p
        /// runs for HP-family encodings (hp, hp_lehninger, hp_thomas_dill, etc.).
        /// The setting is stored in the index and reused automatically at search
        /// time, so you do not repeat it when searching.
        #[arg(long, default_value = "false")]
        remove_low_complexity: bool,

        /// Fit the Karlin-Altschul K that `kmerseek search` uses for E-values, for this
        /// mismatch penalty (the `--extend-mismatch-penalty` a search will pass). K depends
        /// on the alphabet, the seed length, the penalty, the X-drop and the database, so it
        /// is fitted here, on this index, and stored in it. A search with a different
        /// penalty or X-drop refits on the fly.
        #[arg(long, default_value = "2.0")]
        extend_mismatch_penalty: f64,

        /// The `--extend-xdrop` the K fit assumes.
        #[arg(long, default_value = "8.0")]
        extend_xdrop: f64,

        /// How many database sequences to search against the index to fit lambda and K.
        /// ln(regions with score >= S) is a straight line in S whose slope is -lambda and
        /// whose intercept gives K; related pairs bend it upward and are cut off. 0 skips
        /// the fit, and a search then has to fit its own or be given --ka-k.
        #[arg(long, default_value = "200")]
        ka_queries: usize,

        /// Seed for picking the calibration sequences, so the fit is reproducible.
        #[arg(long, default_value = "1")]
        ka_seed: u64,

        /// What the calibration queries are. `database`: the sequences as they are, with
        /// the homolog bend cut off. `shuffled`: residues shuffled, the independent-letter
        /// model, which loses the hydrophobic runs and periodicity real proteins have.
        /// `reversed`: read back to front, which in a hydrophobic/polar alphabet still
        /// matches the forward helices and strands.
        #[arg(long, value_enum, default_value_t = DecoyNull::Database)]
        ka_null: DecoyNull,

        /// For `--ka-null database`: how the reference queries that decide where the fit
        /// stops are made. `shuffled-dipeptide` keeps each pair of neighbouring residues as
        /// often as in the original, so hydrophobic runs survive and only relatives and
        /// periodicity lift the real curve above it; `shuffled` keeps composition only.
        #[arg(long, value_enum, default_value_t = DecoyNull::ShuffledDipeptide)]
        ka_reference: DecoyNull,

        /// Write the survival curve the fit was read from (score, regions with score >= it,
        /// and the fit) to this CSV, for plotting with scripts/plot_ka_survival.py. Written
        /// even when the fit is refused, with the fit columns empty and `fitted` false.
        #[arg(long, value_name = "PATH")]
        ka_survival_out: Option<PathBuf>,
    },
    /// Fit and store the Karlin-Altschul lambda and K of an existing index for one more
    /// mismatch penalty and X-drop. `kmerseek index` fits one pair at build time; a search
    /// with a different `--extend-mismatch-penalty` then either finds a fit stored here or
    /// has to fit its own on every run.
    Calibrate {
        /// Index directory, as written by `kmerseek index`. Opened read-write, so no
        /// search may have it open at the same time.
        #[arg(short, long)]
        target: PathBuf,

        /// The `--extend-mismatch-penalty` the fit is for.
        #[arg(long, default_value = "2.0")]
        extend_mismatch_penalty: f64,

        /// The `--extend-xdrop` the fit assumes.
        #[arg(long, default_value = "8.0")]
        extend_xdrop: f64,

        /// How many database sequences to search against the index to fit lambda and K.
        #[arg(long, default_value = "200")]
        ka_queries: usize,

        /// Seed for picking the calibration sequences, so the fit is reproducible.
        #[arg(long, default_value = "1")]
        ka_seed: u64,

        /// What the calibration queries are; see `kmerseek index --help`.
        #[arg(long, value_enum, default_value_t = DecoyNull::Database)]
        ka_null: DecoyNull,

        /// For `--ka-null database`: how the reference queries are made; see
        /// `kmerseek index --help`.
        #[arg(long, value_enum, default_value_t = DecoyNull::ShuffledDipeptide)]
        ka_reference: DecoyNull,

        /// Write the survival curve the fit was read from to this CSV, for plotting with
        /// scripts/plot_ka_survival.py. Written even when the fit is refused, with the fit
        /// columns empty and `fitted` false.
        #[arg(long, value_name = "PATH")]
        ka_survival_out: Option<PathBuf>,
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

        /// Grow each matched region past its exact k-mer run, charging this much per
        /// encoded position where query and target disagree (+1 per agreeing position).
        /// 0 keeps regions exact, the default. A remote homolog conserves the HP pattern per
        /// position far better than it conserves any 23-residue stretch of it exactly, so an
        /// exact run is treated as a seed and extended with X-drop (see --extend-xdrop).
        /// The region's shared k-mer count and Poisson score still count exact k-mers only;
        /// `region_n_mismatches` reports how many positions inside the region disagree.
        #[arg(long, default_value = "0.0")]
        extend_mismatch_penalty: f64,

        /// Stop extending once the running score has fallen this far below its best.
        #[arg(long, default_value = "8.0")]
        extend_xdrop: f64,

        /// Karlin-Altschul K for `region_evalue` and `region_ka_bits` on extended regions,
        /// with the closed-form lambda per pair. Normally left unset: the lambda and K
        /// fitted when the index was built (for its penalty and X-drop) are used, or, for
        /// another penalty or X-drop, a fit on --ka-queries database sequences runs before
        /// the search. Used only with --extend-mismatch-penalty.
        #[arg(long)]
        ka_k: Option<f64>,

        /// Calibration queries to fit lambda and K on when the index has no fit for this
        /// penalty and X-drop and --ka-k is unset. 0 refuses to search without a fit.
        #[arg(long, default_value = "200")]
        ka_queries: usize,

        /// Seed for picking the calibration sequences.
        #[arg(long, default_value = "1")]
        ka_seed: u64,

        /// What the calibration queries are when a fit runs here; see `kmerseek index --help`.
        #[arg(long, value_enum, default_value_t = DecoyNull::Database)]
        ka_null: DecoyNull,

        /// The reference for `--ka-null database` when a fit runs here; see `kmerseek index --help`.
        #[arg(long, value_enum, default_value_t = DecoyNull::ShuffledDipeptide)]
        ka_reference: DecoyNull,

        /// Chain extended regions on one diagonal at most this many residues apart into one
        /// region scored with Karlin-Altschul sum statistics (Karlin & Altschul 1993). A
        /// domain that a single gapless run cannot cover becomes one call. 0 (default) keeps
        /// every region separate. Used only with --extend-mismatch-penalty.
        #[arg(long, default_value = "0")]
        chain_max_gap: u32,

        /// Largest diagonal shift (net indel) between chained regions. 0 chains only along
        /// one diagonal. Used with --chain-max-gap.
        #[arg(long, default_value = "0")]
        chain_max_shift: u32,

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
    /// Compare one query sequence with one target sequence: list every shared k-mer with
    /// its position in both, and the matched regions they chain into, as JSON. Plot the
    /// output with scripts/visualize_pair.py.
    Pair {
        /// Query FASTA file path; the first record is used unless --query-name is given
        #[arg(short, long)]
        query: PathBuf,

        /// Target FASTA file path; the first record is used unless --target-name is given
        #[arg(short, long)]
        target: PathBuf,

        /// Header of the query record to use, either the whole header or its first token
        /// (e.g. sp|P10415|BCL2_HUMAN)
        #[arg(long)]
        query_name: Option<String>,

        /// Header of the target record to use, either the whole header or its first token
        #[arg(long)]
        target_name: Option<String>,

        /// Output JSON path (optional - will output to stdout if not provided)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// K-mer size
        #[arg(short, long, default_value = "10")]
        ksize: u32,

        /// Reduced amino acid alphabet
        #[arg(short = 'a', long, default_value = "protein20")]
        alphabet: ProteinAlphabet,
    },
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq)]
enum ProteinAlphabet {
    /// The full 20-letter amino acid alphabet: no reduction. Also accepts sourmash's
    /// name for it, `protein`
    #[value(name = "protein20", alias = "protein")]
    Protein,
    /// Dayhoff, 6 classes. Also accepts sourmash's name for it, `dayhoff`
    #[value(name = "dayhoff6", alias = "dayhoff")]
    Dayhoff,
    /// HP Lehninger, 2 classes. Also accepts sourmash's name for it, `hp`, which is the
    /// same partition hashed the same way
    #[value(name = "hp_lehninger2", aliases = ["hp-lehninger2", "hp"])]
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
    /// GBMR4, 4 classes (Solis & Rackovsky 2000; best recall in Peterson et al. 2009)
    #[value(name = "gbmr4")]
    ReducedGbmr4,
    /// POLARITY4, 4 classes on polarity and charge (Ball, Hill & Scott 2014)
    #[value(name = "polarity4")]
    ReducedPolarity4,
    /// WWMJ5, 5 classes (Wang & Wang 1999, Miyazawa-Jernigan contact potentials)
    #[value(name = "wwmj5")]
    ReducedWwmj5,
    /// GBMR7, 7 classes (Solis & Rackovsky 2000)
    #[value(name = "gbmr7")]
    ReducedGbmr7,
    /// FUNCGROUPS8, 8 classes, one per side-chain functional group (Jain et al. 2014)
    #[value(name = "funcgroups8")]
    ReducedFuncGroups8,
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
            ProteinAlphabet::ReducedGbmr4 => "gbmr4",
            ProteinAlphabet::ReducedPolarity4 => "polarity4",
            ProteinAlphabet::ReducedWwmj5 => "wwmj5",
            ProteinAlphabet::ReducedGbmr7 => "gbmr7",
            ProteinAlphabet::ReducedFuncGroups8 => "funcgroups8",
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
            scaled,
            alphabet,
            progress_interval,
            kmer_stats_out,
            stats_only,
            remove_low_complexity,
            extend_mismatch_penalty,
            extend_xdrop,
            ka_queries,
            ka_seed,
            ka_null,
            ka_reference,
            ka_survival_out,
        } => {
            eprintln!("Indexing FASTA file: {}", input.display());

            // Fail on a bad value here, before any database is created.
            let scaled = Scaled::new(scaled)
                .map_err(|message| IndexError::ConfigurationError {
                    field: "scaled".to_string(),
                    message,
                })?
                .get();

            let effective_moltype: &'static str = alphabet.into();

            // --stats-only builds in a scratch directory that is removed on exit. The
            // index is still written there: the streaming indexer sorts k-mers on disk,
            // so there is no in-memory path that could produce the spectrum without it.
            let scratch = if stats_only { Some(tempfile::tempdir()?) } else { None };

            // Determine output path
            let output_path = if let Some(scratch) = &scratch {
                scratch.path().join("stats-only.kmerseek.rocksdb")
            } else if let Some(output) = output {
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
                    effective_moltype,
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
                effective_moltype,
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
                // kmer_stats_out is guaranteed Some here -- clap's
                // `requires = "kmer_stats_out"` enforces it.
                let kmer_stats_out =
                    kmer_stats_out.expect("clap requires kmer_stats_out with stats_only");
                index.save_state_with_kmer_stats(Some(&kmer_stats_out))?;
                drop(index);
                drop(scratch);
                eprintln!("Stats-only run completed successfully (no index persisted).");
            } else {
                // Enable compactions for better read performance
                eprintln!("Optimizing database for read operations...");
                index.enable_compactions()?;

                // Save the index state for loading
                index.save_state_with_kmer_stats(kmer_stats_out.as_deref())?;

                if ka_queries > 0 {
                    let settings = KaCalibrationSettings {
                        mismatch_penalty: extend_mismatch_penalty,
                        xdrop: extend_xdrop,
                        null: ka_null,
                        reference: ka_reference,
                        n_queries: ka_queries,
                        seed: ka_seed,
                    };
                    calibrate_index(index, settings, ka_survival_out.as_deref())?;
                } else {
                    eprintln!(
                        "Skipping the Karlin-Altschul fit (--ka-queries 0); a search will \
                         have to fit lambda and K itself or be given --ka-k."
                    );
                }

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
            threshold,
            min_shared_kmers,
            max_query_pvalue,
            min_region_score,
            max_pvalue,
            remove_low_complexity: remove_low_complexity_arg,
            extend_mismatch_penalty,
            extend_xdrop,
            ka_k,
            ka_queries,
            ka_seed,
            ka_null,
            ka_reference,
            chain_max_gap,
            chain_max_shift,
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
            if extend_mismatch_penalty > 0.0 {
                eprintln!(
                    "  Seed extension: mismatch penalty {}, X-drop {}, chain gap {} shift {}",
                    extend_mismatch_penalty, extend_xdrop, chain_max_gap, chain_max_shift
                );
            } else {
                eprintln!("  Seed extension: off (regions are exact runs)");
            }
            eprintln!("  Verbose output: {}", verbose);
            eprintln!("  Query is pre-indexed: {}\n---", query_is_index);

            use kmerseek::search::SearchFilters;
            let filters = SearchFilters {
                threshold,
                min_shared_kmers,
                max_query_pvalue,
                min_region_score,
                skip_self_matches: false,
            };

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
            if extend_mismatch_penalty > 0.0 {
                use kmerseek::search::ExtensionParams;
                let settings = KaCalibrationSettings {
                    mismatch_penalty: extend_mismatch_penalty,
                    xdrop: extend_xdrop,
                    null: ka_null,
                    reference: ka_reference,
                    n_queries: ka_queries,
                    seed: ka_seed,
                };
                let (ka, source) = searcher.resolve_ka(ka_k, settings)?;
                eprintln!(
                    "  Karlin-Altschul: K {:.4}, lambda scale {:.3} ({source})",
                    ka.k, ka.lambda_scale
                );
                if let KaSource::Index(fit) | KaSource::Fitted(fit) = &source {
                    warn_on_short_fit(fit);
                }
                searcher.set_extension(Some(ExtensionParams {
                    mismatch_penalty: extend_mismatch_penalty,
                    xdrop: extend_xdrop,
                    ka_k: ka.k,
                    ka_lambda_scale: ka.lambda_scale,
                    chain_max_gap,
                    chain_max_shift,
                }));
            }

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
                // Sketches from two indexes are compared as they are, so the two must agree
                // on every sketch parameter. A scaled mismatch would otherwise reach
                // find_matched_regions, which asserts it.
                let query_params =
                    (query_index.ksize(), query_index.scaled(), query_index.moltype());
                let target_params = (final_ksize, final_scaled, detected_moltype.as_str());
                if query_params != target_params {
                    return Err(IndexError::ValidationError {
                        message: format!(
                            "Query index (ksize={}, scaled={}, alphabet={}) was not built with \
                             the target's parameters (ksize={}, scaled={}, alphabet={})",
                            query_params.0,
                            query_params.1,
                            query_params.2,
                            target_params.0,
                            target_params.1,
                            target_params.2,
                        ),
                    });
                }
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
        Commands::Calibrate {
            target,
            extend_mismatch_penalty,
            extend_xdrop,
            ka_queries,
            ka_seed,
            ka_null,
            ka_reference,
            ka_survival_out,
        } => {
            if ka_queries == 0 {
                return Err(
                    anyhow::anyhow!("--ka-queries must be at least 1 to fit anything").into()
                );
            }
            eprintln!("Opening index {} read-write", target.display());
            let index = ProteomeIndex::open_for_calibration(&target)?;
            let settings = KaCalibrationSettings {
                mismatch_penalty: extend_mismatch_penalty,
                xdrop: extend_xdrop,
                null: ka_null,
                reference: ka_reference,
                n_queries: ka_queries,
                seed: ka_seed,
            };
            calibrate_index(index, settings, ka_survival_out.as_deref())?;
        }
        Commands::Pair { query, target, query_name, target_name, output, ksize, alphabet } => {
            run_pair(&query, &target, query_name, target_name, output, ksize, alphabet.into())?;
        }
    }

    Ok(())
}

fn run_pair(
    query: &Path,
    target: &Path,
    query_name: Option<String>,
    target_name: Option<String>,
    output: Option<PathBuf>,
    ksize: u32,
    moltype: &str,
) -> IndexResult<()> {
    let query = pair::read_record(query, query_name.as_deref())?;
    let target = pair::read_record(target, target_name.as_deref())?;
    let report = pair::compare_pair(&query, &target, ksize, moltype)?;
    eprintln!(
        "{} shared {}-mers in {} matched regions ({})",
        report.shared_kmers.len(),
        report.ksize,
        report.regions.len(),
        report.moltype
    );
    let json = report.to_json()?;
    match output {
        Some(path) => std::fs::write(&path, json)?,
        None => println!("{json}"),
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
    // "hp_<name>" moltype. Alphabet::from_moltype() still parses those, so normalizing
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
        "gbmr4" => ProteinAlphabet::ReducedGbmr4,
        "polarity4" => ProteinAlphabet::ReducedPolarity4,
        "wwmj5" => ProteinAlphabet::ReducedWwmj5,
        "gbmr7" => ProteinAlphabet::ReducedGbmr7,
        "funcgroups8" => ProteinAlphabet::ReducedFuncGroups8,
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
                     gbmr4, polarity4, wwmj5, gbmr7, funcgroups8, sdm12, mmseqs12, \
                     wass14, hsdm17, uniprot18",
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
    // In practice, users should not specify --alphabet and let it autodetect.
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

/// Say so when the homolog excess left the fit fewer bins than `FIT_WINDOW`: the slope is
/// then read from the seed end of the curve, where the seed requirement still shapes it.
fn warn_on_short_fit(fit: &KaCalibration) {
    if fit.n_fit_points() < kmerseek::karlin_altschul::FIT_WINDOW {
        eprintln!(
            "  WARNING: the fit has only {} bins (x {:.1}..{:.1}) below the relatives at x {}. \
             Related sequences are dense in this database; the slope is read close to the \
             seed. More --ka-queries gives a second opinion.",
            fit.n_fit_points(),
            fit.x_range().0,
            fit.x_range().1,
            fit.bend_score
                .map_or("none".to_string(), |b| format!("{:.1}", b as f64 * fit.bin_width))
        );
    }
}

/// Fit lambda and K on `n_queries` calibration queries of the index just built and store
/// the fit in the index. Also prints the closed-form lambda and K for the database's own
/// composition, so the effect of the seed requirement and of real sequence structure on
/// each is visible.
fn calibrate_index(
    index: ProteomeIndex,
    settings: KaCalibrationSettings,
    survival_out: Option<&std::path::Path>,
) -> IndexResult<()> {
    use kmerseek::karlin_altschul::karlin_altschul_k_theory;
    let KaCalibrationSettings { mismatch_penalty, xdrop, null, reference, n_queries, .. } =
        settings;
    eprintln!(
        "Fitting Karlin-Altschul lambda and K on {n_queries} {null} sequences (penalty {mismatch_penalty}, X-drop {xdrop}){}...",
        if null == DecoyNull::Database {
            format!(", censored against the same sequences {reference}")
        } else {
            String::new()
        }
    );
    let mut searcher = ProteinSearcher::new(index)?;
    let report = searcher.calibrate_ka(settings)?;
    let theory_k = karlin_altschul_k_theory(report.match_probability, mismatch_penalty)
        .map_or("none".to_string(), |k| format!("{k:.4}"));
    eprintln!(
        "  Closed form at the database's own match probability {:.3}: K {theory_k} (independent positions, no seed, one lambda for every pair)",
        report.match_probability
    );
    match &report.fitted {
        Some(fit) => {
            eprintln!("  {}", KaSource::Fitted(fit.clone()));
            warn_on_short_fit(fit);
            searcher.index().put_ka_calibration(fit)?;
            eprintln!(
                "  Stored in the index for --extend-mismatch-penalty {mismatch_penalty} --extend-xdrop {xdrop}"
            );
        }
        None => eprintln!(
            "  {} queries gave {} regions but fewer than {} score bins above the peak with {} \
             regions each, too few to fit; nothing stored. A search will have to fit its own \
             lambda and K (--ka-queries) or be given --ka-k.",
            report.n_queries,
            report.n_regions,
            kmerseek::karlin_altschul::MIN_FIT_POINTS,
            kmerseek::karlin_altschul::MIN_BIN_COUNT
        ),
    }
    // Written after the verdict, and whether or not there was a fit: a refused fit is the
    // one whose histogram most needs looking at.
    if let Some(path) = survival_out {
        write_survival_csv(path, &report, &settings)?;
        eprintln!("  Survival curve written to {}", path.display());
    }
    Ok(())
}

fn write_survival_csv(
    path: &std::path::Path,
    report: &KaCalibrationReport,
    settings: &KaCalibrationSettings,
) -> IndexResult<()> {
    let mut w = csv::Writer::from_path(path)?;
    w.write_record([
        "x",
        "n_regions_at_least",
        "fitted_n_regions_at_least",
        "reference_n_regions_at_least",
        "in_fit",
        "slope",
        "k",
        "lambda_analytic",
        "match_probability",
        "null",
        "reference",
        "mismatch_penalty",
        "xdrop",
        "n_queries",
        "query_residues",
        "database_kmers",
        "bin_width",
        "fitted",
    ])?;
    // The fit is a line through ln(regions in the bin at x); its survival is the same line
    // divided by (1 - e^(-slope w)). Without a fit those three columns stay empty and the
    // curve itself is still written.
    let line = report.fitted.as_ref().map(|fit| {
        let per_bin = 1.0 - (-fit.slope * fit.bin_width).exp();
        let ln_intercept =
            (fit.k * fit.query_residues as f64 * fit.database_kmers as f64 * per_bin).ln();
        (fit, per_bin, ln_intercept)
    });
    let reference: std::collections::HashMap<i64, u64> =
        report.reference_survival.iter().copied().collect();
    let reference_null = (settings.null == DecoyNull::Database).then_some(settings.reference);
    for &(bin, count) in &report.survival {
        let x = bin as f64 * BIN_WIDTH;
        let (fitted, in_fit, slope, k) = match line {
            Some((fit, per_bin, ln_intercept)) => (
                format!("{:.3}", (ln_intercept - fit.slope * x).exp() / per_bin),
                (fit.score_lo <= bin && bin <= fit.score_hi).to_string(),
                fit.slope.to_string(),
                fit.k.to_string(),
            ),
            None => (String::new(), String::new(), String::new(), String::new()),
        };
        w.write_record([
            format!("{x:.3}"),
            count.to_string(),
            fitted,
            reference.get(&bin).map_or(String::new(), |r| r.to_string()),
            in_fit,
            slope,
            k,
            report.lambda_analytic.to_string(),
            report.match_probability.to_string(),
            settings.null.to_string(),
            reference_null.map_or(String::new(), |r| r.to_string()),
            settings.mismatch_penalty.to_string(),
            settings.xdrop.to_string(),
            report.n_queries.to_string(),
            report.query_residues.to_string(),
            report.database_kmers.to_string(),
            BIN_WIDTH.to_string(),
            report.fitted.is_some().to_string(),
        ])?;
    }
    w.flush()?;
    Ok(())
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

    // Query sketches must use the database's scaled factor, or the FracMinHash cutoffs
    // disagree and shared k-mers go unseen. There is no --scaled on search, so the
    // detected value is the only one.
    let final_scaled = detected_scaled;

    // Validate and assign encoding
    let final_alphabet = assign_encoding(user_encoding, detected_moltype)?;

    Ok((final_ksize, final_scaled, final_alphabet))
}
