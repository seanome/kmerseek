pub const TEST_FASTA_GZ: &str =
    "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";
pub const TEST_CED9_FASTA: &str = "tests/testdata/fasta/ced9.fasta";
pub const TEST_BLC2_FASTA: &str = "tests/testdata/fasta/bcl2.fasta";
/// Two decoys from notebook 344 (seanome/2024-kmerseek-analysis PR 174): dipeptide shuffles,
/// replicate 2, of DSPP and of elastin, rebuilt from the QfO 2020_04 human proteome with that
/// notebook's seeds. The search there called a chance region between each and a real protein
/// significant: DSPP's against TBC1D10B's KE repeat, elastin's against FMN2's proline-rich
/// stretch.
pub const NB344_DECOYS_FASTA: &str =
    "tests/testdata/fasta/nb344_decoys_dspp_eln_dipeptide_r2.fasta";
/// TBC1D10B and FMN2 from the same proteome: the two targets of `NB344_DECOYS_FASTA`.
pub const NB344_TARGETS_FASTA: &str = "tests/testdata/fasta/nb344_targets_tbc1d10b_fmn2.fasta";
/// The 25 BCL-2-like proteins of TEST_FASTA_GZ, each shuffled 20 times with its 2-mer
/// (dipeptide) counts kept, so a decoy has the composition of a real protein and nothing
/// else in common with any query. Made with
/// `python shuffle_fasta_2mer.py <TEST_FASTA_GZ> -n 20 --seed 1 | gzip`; the Python test
/// `scripts/test_shuffle_fasta_2mer.py` regenerates it and checks it is unchanged.
pub const TEST_DECOYS_2MER_GZ: &str =
    "tests/testdata/fasta/bcl2_25_shuffled_2mer_20x_seed1.fasta.gz";
/// BCL2_HUMAN (UniProt P10415) residues 1-32, the first 32 residues of the protein.
pub const BCL2_1_32: &str = "MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDA";
/// k-mers in `BCL2_1_32` at k=5: 32 - 5 + 1 windows, all distinct.
pub const BCL2_1_32_KMERS_K5: usize = 28;
pub const TEST_FASTA_ZST: &str = "tests/testdata/fasta/test_compression.fasta.zst";

pub const TEST_FASTA_CONTENT: &str =
    ">test_protein1\nPLANTANDANIMALGENQMES\n>test_protein2\nLIVINGALIVE";
pub const TEST_KMER: &str = "LIVINGALIVE";
pub const TEST_PROTEIN: &str = "PLANTANDANIMALGENQMES";

// Contains invalid character '1' which is not a valid amino acid
pub const TEST_PROTEIN_INVALID: &str = "PLANTANDANIMALGEN1MES";

// Mixed case sequences for testing uppercasing functionality
pub const TEST_FASTA_MIXED_CASE_CONTENT: &str =
    ">test_protein_mixed1\nmAaGgCcTt\n>test_protein_mixed2\nmAaGgCcTtNnRrSsVvWwYyHhKkDdEeFfPpQqIiLl";
