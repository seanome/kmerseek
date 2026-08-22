//! Multi-letter reduced amino acid alphabets from the fold-recognition literature.
//!
//! Unlike [`crate::hp_alphabets`], which splits the 20 residues two or three ways,
//! these alphabets keep 4 to 18 classes. Peterson et al. (2009) benchmarked over 150
//! published clustering schemes against DALI fold assignments and reported GBMR4,
//! SDM12 and HSDM17 as the best performers on recall, AUC and mean pooled precision
//! respectively. Ieremie et al. (2024) reused those three and added GBMR7, WWMJ5,
//! MMSEQS12, WASS14 and UNIPROT18 when testing how alphabet reduction affects protein
//! language models. Where the two papers overlap, their partitions agree.
//!
//! References:
//!   Peterson, E. L., Kondev, J., Theriot, J. A. & Phillips, R. (2009). Reduced amino
//!   acid alphabets exhibit an improved sensitivity and selectivity in fold assignment.
//!   Bioinformatics 25(11):1356-1362. doi:10.1093/bioinformatics/btp164. Table 2.
//!
//!   Ieremie, I., Ewing, R. M. & Niranjan, M. (2024). Protein language models meet
//!   reduced amino acid alphabets. Bioinformatics 40(2):btae061.
//!   doi:10.1093/bioinformatics/btae061. Table 1.

use std::collections::HashMap;
use std::sync::LazyLock;

/// The 20 canonical residues, sorted. Every alphabet below must cover all of these and
/// nothing else.
const CANONICAL_AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

/// Reduced alphabets with more than the two or three classes of an HP table.
///
/// The numeric suffix in each name is the number of classes, following the convention
/// of both source papers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReducedAlphabet {
    /// Solis & Rackovsky (2000), 4 classes. Top recall at 0.01 EPQ in Peterson et al.
    Gbmr4,
    /// Wang & Wang (1999), 5 classes, derived from the Miyazawa-Jernigan matrix.
    Wwmj5,
    /// Solis & Rackovsky (2000), 7 classes.
    Gbmr7,
    /// Prlic et al. (2000) structure-derived substitution matrix, 12 classes.
    /// Top AUC in Peterson et al.
    Sdm12,
    /// Steinegger & Soding (2018), 12 classes, from MMseqs2 sequence clustering.
    Mmseqs12,
    /// Ieremie et al. (2024), 14 classes, clustered on hydrophobicity.
    Wass14,
    /// Prlic et al. (2000) homologous structure-derived matrix, 17 classes.
    /// Top mean pooled precision in Peterson et al.
    Hsdm17,
    /// Ieremie et al. (2024), 18 classes, from clusters learned by a protein language
    /// model trained on the full alphabet.
    Uniprot18,
}

impl ReducedAlphabet {
    pub fn table(&self) -> &'static HashMap<u8, u8> {
        match self {
            Self::Gbmr4 => &GBMR4,
            Self::Wwmj5 => &WWMJ5,
            Self::Gbmr7 => &GBMR7,
            Self::Sdm12 => &SDM12,
            Self::Mmseqs12 => &MMSEQS12,
            Self::Wass14 => &WASS14,
            Self::Hsdm17 => &HSDM17,
            Self::Uniprot18 => &UNIPROT18,
        }
    }

    /// The residue clusters, in the order given by the source paper.
    pub fn clusters(&self) -> &'static [&'static str] {
        match self {
            Self::Gbmr4 => GBMR4_CLUSTERS,
            Self::Wwmj5 => WWMJ5_CLUSTERS,
            Self::Gbmr7 => GBMR7_CLUSTERS,
            Self::Sdm12 => SDM12_CLUSTERS,
            Self::Mmseqs12 => MMSEQS12_CLUSTERS,
            Self::Wass14 => WASS14_CLUSTERS,
            Self::Hsdm17 => HSDM17_CLUSTERS,
            Self::Uniprot18 => UNIPROT18_CLUSTERS,
        }
    }

    /// Number of classes, matching the numeric suffix in the alphabet's name.
    pub fn size(&self) -> usize {
        self.clusters().len()
    }

    /// Short identifier used in output filenames and Nextflow traces.
    pub fn name(&self) -> &'static str {
        match self {
            Self::Gbmr4 => "gbmr4",
            Self::Wwmj5 => "wwmj5",
            Self::Gbmr7 => "gbmr7",
            Self::Sdm12 => "sdm12",
            Self::Mmseqs12 => "mmseqs12",
            Self::Wass14 => "wass14",
            Self::Hsdm17 => "hsdm17",
            Self::Uniprot18 => "uniprot18",
        }
    }

    pub fn all() -> &'static [ReducedAlphabet] {
        &[
            ReducedAlphabet::Gbmr4,
            ReducedAlphabet::Wwmj5,
            ReducedAlphabet::Gbmr7,
            ReducedAlphabet::Sdm12,
            ReducedAlphabet::Mmseqs12,
            ReducedAlphabet::Wass14,
            ReducedAlphabet::Hsdm17,
            ReducedAlphabet::Uniprot18,
        ]
    }

    /// Moltype string stored in the index (e.g. `"sdm12"`).
    ///
    /// These are the names the source papers use, digit included, so a moltype is directly
    /// citeable: searching for SDM12 or GBMR4 finds the paper it came from.
    pub fn to_moltype(&self) -> String {
        self.name().to_string()
    }

    /// Parse from a moltype string. Returns `None` for unrecognized strings.
    pub fn from_moltype(s: &str) -> Option<ReducedAlphabet> {
        Self::all().iter().find(|a| a.name() == s).copied()
    }
}

/// Build a residue-to-symbol table from a list of clusters.
///
/// WHY the first residue as the symbol: it keeps the encoded sequence readable, so
/// SDM12's `LIVM` class shows up as `l` rather than an opaque index, and the encoded
/// string can be eyeballed against the source residues. The assertions below enforce
/// that this is unambiguous. Two clusters starting with the same residue would merge
/// into one class without any error.
fn build_reduced(clusters: &[&str]) -> HashMap<u8, u8> {
    let mut m = HashMap::with_capacity(21);
    for cluster in clusters {
        let symbol = cluster.as_bytes()[0].to_ascii_lowercase();
        for &residue in cluster.as_bytes() {
            m.insert(residue, symbol);
        }
    }

    let mut covered: Vec<u8> = m.keys().copied().collect();
    covered.sort_unstable();
    assert_eq!(
        covered, CANONICAL_AA,
        "reduced alphabet must cover each of the 20 canonical amino acids exactly once, got {:?}",
        clusters
    );

    let mut symbols: Vec<u8> = m.values().copied().collect();
    symbols.sort_unstable();
    symbols.dedup();
    assert_eq!(
        symbols.len(),
        clusters.len(),
        "two clusters share a leading residue, so they would collapse into one class: {:?}",
        clusters
    );

    m.insert(b'*', b'*'); // stop codon pass-through
    m
}

// -----------------------------------------------------------------------------
// GBMR4 — Solis & Rackovsky (2000), 4 classes.
//
// Glycine and proline are singled out as structurally unlike anything else; the
// remaining two classes are a hydrophobic (YFLIVMCWH) and a polar (ADKERNTSQ) group,
// which makes GBMR4 a refinement of the plain HP split. Best recall at 0.01 errors
// per query of every alphabet Peterson et al. tested.
// -----------------------------------------------------------------------------
const GBMR4_CLUSTERS: &[&str] = &["ADKERNTSQ", "YFLIVMCWH", "G", "P"];
static GBMR4: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| build_reduced(GBMR4_CLUSTERS));

// -----------------------------------------------------------------------------
// WWMJ5 — Wang & Wang (1999), 5 classes.
//
// Derived by minimizing the mismatch within the Miyazawa-Jernigan contact potential
// matrix. Close to the five-letter alphabet (I, K, E, A, G) that Riddle et al. (1997)
// found sufficient to fold the SH3 domain.
// -----------------------------------------------------------------------------
const WWMJ5_CLUSTERS: &[&str] = &["CMFILVWY", "ATH", "GP", "DE", "SNQRK"];
static WWMJ5: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| build_reduced(WWMJ5_CLUSTERS));

// -----------------------------------------------------------------------------
// GBMR7 — Solis & Rackovsky (2000), 7 classes.
//
// Same maximum-mutual-information construction as GBMR4 at a finer level. Note that
// it is not a refinement of GBMR4: the large AEFIKLMQRVWY class mixes residues that
// GBMR4 keeps on opposite sides of its hydrophobic/polar split.
// -----------------------------------------------------------------------------
const GBMR7_CLUSTERS: &[&str] = &["DN", "AEFIKLMQRVWY", "CH", "T", "S", "G", "P"];
static GBMR7: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| build_reduced(GBMR7_CLUSTERS));

// -----------------------------------------------------------------------------
// SDM12 — Prlic et al. (2000) structure-derived matrix, 12 classes.
//
// Best AUC in Peterson et al., and derivable from GBMR4 by splitting classes without
// moving any residue across an existing boundary. Keeps acidic/basic (KER), polar
// (TSQ), aromatic (YF) and aliphatic (LIVM) classes. Two groupings run against
// intuition: aspartic acid sits outside the KER class, and methionine sits inside the
// otherwise aliphatic LIVM class.
// -----------------------------------------------------------------------------
const SDM12_CLUSTERS: &[&str] =
    &["A", "D", "KER", "N", "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"];
static SDM12: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| build_reduced(SDM12_CLUSTERS));

// -----------------------------------------------------------------------------
// MMSEQS12 — Steinegger & Soding (2018), 12 classes.
//
// Built for fast clustering of large sequence sets in MMseqs2 rather than for fold
// recognition, so it is the natural size-matched comparison for SDM12.
// -----------------------------------------------------------------------------
const MMSEQS12_CLUSTERS: &[&str] =
    &["AST", "LM", "IV", "KR", "EQ", "ND", "FY", "C", "G", "H", "P", "W"];
static MMSEQS12: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| build_reduced(MMSEQS12_CLUSTERS));

// -----------------------------------------------------------------------------
// WASS14 — Ieremie et al. (2024), 14 classes.
//
// Clusters residues with similar solvent-accessibility-derived hydrophobicity. Pairs
// residues that most other schemes separate (W with M, D with I), which is why it
// distinguishes point mutations better than same-size alphabets but carries less
// evolutionary signal.
// -----------------------------------------------------------------------------
const WASS14_CLUSTERS: &[&str] =
    &["WM", "DI", "P", "C", "AV", "K", "T", "RE", "G", "L", "Y", "SH", "F", "NQ"];
static WASS14: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| build_reduced(WASS14_CLUSTERS));

// -----------------------------------------------------------------------------
// HSDM17 — Prlic et al. (2000) homologous structure-derived matrix, 17 classes.
//
// Best mean pooled precision in Peterson et al. Refines SDM12 down to only the
// strongest associations: acidic/basic collapses to KE and aliphatic to LIV, with
// every other residue on its own.
// -----------------------------------------------------------------------------
const HSDM17_CLUSTERS: &[&str] =
    &["A", "D", "KE", "R", "N", "T", "S", "Q", "Y", "F", "LIV", "M", "C", "W", "H", "G", "P"];
static HSDM17: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| build_reduced(HSDM17_CLUSTERS));

// -----------------------------------------------------------------------------
// UNIPROT18 — Ieremie et al. (2024), 18 classes.
//
// Read off the clustering that a protein language model trained on the full alphabet
// arrived at on its own: only glutamate/proline and histidine/leucine merge. The
// mildest reduction here, and the reference point for how much a reduction costs.
// -----------------------------------------------------------------------------
const UNIPROT18_CLUSTERS: &[&str] =
    &["A", "R", "N", "D", "C", "Q", "EP", "G", "HL", "I", "K", "M", "F", "S", "T", "W", "Y", "V"];
static UNIPROT18: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| build_reduced(UNIPROT18_CLUSTERS));

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_canonical(alphabet: ReducedAlphabet) -> String {
        let table = alphabet.table();
        CANONICAL_AA.iter().map(|r| *table.get(r).expect("canonical residue") as char).collect()
    }

    #[test]
    fn test_gbmr4_partition() {
        // ADKERNTSQ -> a, YFLIVMCWH -> y, G -> g, P -> p
        assert_eq!(encode_canonical(ReducedAlphabet::Gbmr4), "ayaaygyyayyapaaaayyy");
    }

    #[test]
    fn test_wwmj5_partition() {
        // CMFILVWY -> c, ATH -> a, GP -> g, DE -> d, SNQRK -> s
        assert_eq!(encode_canonical(ReducedAlphabet::Wwmj5), "acddcgacsccsgsssaccc");
    }

    #[test]
    fn test_gbmr7_partition() {
        // DN -> d, AEFIKLMQRVWY -> a, CH -> c, T -> t, S -> s, G -> g, P -> p
        assert_eq!(encode_canonical(ReducedAlphabet::Gbmr7), "acdaagcaaaadpaastaaa");
    }

    #[test]
    fn test_sdm12_partition() {
        // A a | D d | KER k | N n | TSQ t | YF y | LIVM l | C c | W w | H h | G g | P p
        assert_eq!(encode_canonical(ReducedAlphabet::Sdm12), "acdkyghlkllnptkttlwy");
    }

    #[test]
    fn test_mmseqs12_partition() {
        // AST a | LM l | IV i | KR k | EQ e | ND n | FY f | C c | G g | H h | P p | W w
        assert_eq!(encode_canonical(ReducedAlphabet::Mmseqs12), "acnefghikllnpekaaiwf");
    }

    #[test]
    fn test_wass14_partition() {
        // WM w | DI d | P p | C c | AV a | K k | T t | RE r | G g | L l | Y y | SH s | F f | NQ n
        assert_eq!(encode_canonical(ReducedAlphabet::Wass14), "acdrfgsdklwnpnrstawy");
    }

    #[test]
    fn test_hsdm17_partition() {
        // KE k | LIV l | every other residue on its own
        assert_eq!(encode_canonical(ReducedAlphabet::Hsdm17), "acdkfghlklmnpqrstlwy");
    }

    #[test]
    fn test_uniprot18_partition() {
        // EP e | HL h | every other residue on its own
        assert_eq!(encode_canonical(ReducedAlphabet::Uniprot18), "acdefghikhmneqrstvwy");
    }

    /// `build_reduced` asserts full coverage and distinct symbols, so constructing every
    /// table is itself the check that no cluster list drifts out of shape.
    #[test]
    fn test_every_alphabet_builds_and_has_declared_size() {
        let expected_sizes = [4usize, 5, 7, 12, 12, 14, 17, 18];
        for (alphabet, expected) in ReducedAlphabet::all().iter().zip(expected_sizes) {
            let table = alphabet.table();
            // 20 canonical residues plus the stop-codon pass-through.
            assert_eq!(table.len(), 21, "{}", alphabet.name());
            assert_eq!(table.get(&b'*'), Some(&b'*'), "{}", alphabet.name());
            assert_eq!(alphabet.size(), expected, "{}", alphabet.name());
        }
    }

    /// The numeric suffix in the name has to be the class count, since that is how both
    /// papers and every downstream filename refer to these alphabets.
    #[test]
    fn test_name_suffix_matches_class_count() {
        for alphabet in ReducedAlphabet::all() {
            let digits: String = alphabet.name().chars().filter(|c| c.is_ascii_digit()).collect();
            assert_eq!(digits.parse::<usize>().unwrap(), alphabet.size(), "{}", alphabet.name());
        }
    }

    /// The moltype is the paper's own name, digit included.
    #[test]
    fn test_moltype_round_trip() {
        for alphabet in ReducedAlphabet::all() {
            let moltype = alphabet.to_moltype();
            assert_eq!(moltype, alphabet.name());
            assert_eq!(ReducedAlphabet::from_moltype(&moltype), Some(*alphabet));
        }
        assert_eq!(ReducedAlphabet::from_moltype("sdm12"), Some(ReducedAlphabet::Sdm12));
        assert_eq!(ReducedAlphabet::from_moltype("hsdm17"), Some(ReducedAlphabet::Hsdm17));
    }

    #[test]
    fn test_from_moltype_rejects_other_moltypes() {
        for moltype in ["protein20", "dayhoff6", "hp_lehninger2", "reduced_sdm12", "nope"] {
            assert_eq!(ReducedAlphabet::from_moltype(moltype), None, "{moltype}");
        }
    }

    /// GBMR4, SDM12 and HSDM17 appear in both source papers, and Peterson et al. note
    /// that each is a refinement of the previous one: splitting classes only, never
    /// moving a residue across an existing boundary. Two residues sharing a class in the
    /// finer alphabet must therefore share one in the coarser alphabet too.
    #[test]
    fn test_hsdm17_refines_sdm12_refines_gbmr4() {
        let chain = [
            (ReducedAlphabet::Hsdm17, ReducedAlphabet::Sdm12),
            (ReducedAlphabet::Sdm12, ReducedAlphabet::Gbmr4),
        ];
        for (finer, coarser) in chain {
            let (fine, coarse) = (finer.table(), coarser.table());
            for &a in CANONICAL_AA {
                for &b in CANONICAL_AA {
                    if fine[&a] == fine[&b] {
                        assert_eq!(
                            coarse[&a],
                            coarse[&b],
                            "{} groups {} with {} but {} splits them",
                            finer.name(),
                            a as char,
                            b as char,
                            coarser.name()
                        );
                    }
                }
            }
        }
    }
}
