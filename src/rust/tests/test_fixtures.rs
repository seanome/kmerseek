pub const TEST_FASTA_GZ: &str =
    "tests/testdata/fasta/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";
pub const TEST_CED9_FASTA: &str = "tests/testdata/fasta/ced9.fasta";
pub const TEST_BLC2_FASTA: &str = "tests/testdata/fasta/bcl2.fasta";
/// The 25 BCL-2-like proteins of TEST_FASTA_GZ, each shuffled 20 times with its 2-mer
/// (dipeptide) counts kept, so a decoy has the composition of a real protein and nothing
/// else in common with any query. `python shuffle_fasta_2mer.py <TEST_FASTA_GZ> -n 20 --seed 1`.
pub const TEST_DECOYS_2MER_GZ: &str =
    "tests/testdata/fasta/bcl2_25_shuffled_2mer_20x_seed1.fasta.gz";
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
