pub const TEST_FASTA_GZ: &str =
    "tests/testdata/index/bcl2_first25_uniprotkb_accession_O43236_OR_accession_2025_02_06.fasta.gz";
pub const TEST_FASTA: &str = "tests/testdata/fasta/ced9.fasta";

pub const TEST_FASTA_CONTENT: &str =
    ">test_protein1\nPLANTANDANIMALGENQMES\n>test_protein2\nLIVINGALIVE";
pub const TEST_KMER: &str = "LIVINGALIVE";
pub const TEST_PROTEIN: &str = "PLANTANDANIMALGENQMES";

// Contains invalid character '1' which is not a valid amino acid
pub const TEST_PROTEIN_INVALID: &str = "PLANTANDANIMALGEN1MES";
