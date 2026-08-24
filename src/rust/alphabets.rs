//! Every alphabet kmerseek can index with, as one [`Alphabet`] enum.
//!
//! An alphabet is a partition of the 20 canonical residues into classes, stored as a
//! residue-to-class table. `protein20` is the degenerate case that keeps all 20.
//!
//! Two families make up most of the list. The HP tables split the residues into
//! hydrophobic and polar, differing only on the borderline residues C, G, P, W and Y, plus
//! one three-class variant that gives cysteine its own symbol. The multi-letter alphabets
//! keep 4 to 18 classes and come from the fold-recognition literature.
//!
//! [`alphabet_table`] gives the table to pre-encode a sequence with. It returns `None` for
//! `protein20`, `dayhoff6` and `hp_lehninger2`, which sourmash encodes itself; those go
//! through [`crate::hash_functions`] instead.

use std::collections::HashMap;
use std::sync::LazyLock;

/// The 20 canonical residues, sorted. Every alphabet here must cover all of these and
/// nothing else, whichever family it belongs to.
const CANONICAL_AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

// ----------------------------------------------------------------------------
// HP alphabets: two classes, or three where cysteine is split out.
// HP alphabet mappings for the Kmerseek robustness sweep.
//
// Each alphabet partitions the 20 canonical amino acids into
// hydrophobic (h) and polar (p) classes. The alphabets differ on
// borderline residues — principally C, G, P, W, and Y.
// ----------------------------------------------------------------------------

/// Translate a sourmash moltype into the kmerseek name for the same alphabet.
///
/// sourmash calls an alphabet a moltype, and writes three: `protein`, `dayhoff` and `hp`.
/// Each names a partition kmerseek also has. kmerseek hashes all three through sourmash's
/// own hash function rather than pre-encoding them, so a sketch or index carrying one of
/// these names already holds hashes kmerseek can read.
///
/// Every other string passes through unchanged, including kmerseek's own earlier spellings
/// (`raw`, `hp_<name>` without a class count). Those were never sourmash names, so there is
/// nothing to stay compatible with and they are rejected further along.
pub fn canonical_moltype(moltype: &str) -> &str {
    match moltype {
        "protein" => "protein20",
        "dayhoff" => "dayhoff6",
        "hp" => "hp_lehninger2",
        other => other,
    }
}

/// Every alphabet kmerseek can index with.
///
/// Each variant is named for its moltype, so `Sdm12` is `sdm12` and `HpLehninger2` is
/// `hp_lehninger2`. The trailing digit is the class count, and [`Self::size`] returns it.
///
/// Three of these are encoded by sourmash rather than by a table here, which
/// [`Self::uses_sourmash_encoder`] reports: `protein20`, `dayhoff6` and `hp_lehninger2`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Alphabet {
    /// All 20 residues, no reduction.
    Protein20,
    /// Dayhoff's six groups.
    Dayhoff6,

    /// Lehninger's hydrophobic/polar split, which is also sourmash's `aa_to_hp`.
    HpLehninger2,
    /// Thomas & Dill 1996: C hydrophobic, G and P polar.
    HpThomasDill2,
    /// Kyte & Doolittle 1982, split at hydropathy > 0: W and Y polar.
    HpKyteDoolittle2,
    /// Thomas-Dill with C moved to polar.
    HpThomasDillNoC2,
    /// Lehninger with C moved to hydrophobic.
    HpLehningerCNonpolar2,
    /// Lehninger with C in a third class of its own.
    HpLehningerHpc3,
    /// Physical Biology of the Cell, 1st edition.
    HpPBotC1stEd2,
    /// Negative control: the h/p split is randomized, scrambling the hydrophobicity signal.
    HpRandomControl2,
    /// One of ten pre-computed random controls, for a null distribution. Seed must be 1-10.
    HpRandomControl2Seed(u64),

    /// Solis & Rackovsky 2000, 4 classes.
    Gbmr4,
    /// Wang & Wang 1999, 5 classes.
    Wwmj5,
    /// Solis & Rackovsky 2000, 7 classes.
    Gbmr7,
    /// Prlic et al. 2000 structure-derived matrix, 12 classes.
    Sdm12,
    /// Steinegger & Soding 2018, 12 classes.
    Mmseqs12,
    /// Ieremie et al. 2024, 14 classes, clustered on hydrophobicity.
    Wass14,
    /// Prlic et al. 2000 homologous structure-derived matrix, 17 classes.
    Hsdm17,
    /// Ieremie et al. 2024, 18 classes, learned by a protein language model.
    Uniprot18,
}

impl Alphabet {
    /// Every alphabet that has a fixed moltype, so not the seeded random controls.
    pub fn all() -> &'static [Alphabet] {
        &[
            Alphabet::Protein20,
            Alphabet::Dayhoff6,
            Alphabet::HpLehninger2,
            Alphabet::HpThomasDill2,
            Alphabet::HpKyteDoolittle2,
            Alphabet::HpThomasDillNoC2,
            Alphabet::HpLehningerCNonpolar2,
            Alphabet::HpLehningerHpc3,
            Alphabet::HpPBotC1stEd2,
            Alphabet::HpRandomControl2,
            Alphabet::Gbmr4,
            Alphabet::Wwmj5,
            Alphabet::Gbmr7,
            Alphabet::Sdm12,
            Alphabet::Mmseqs12,
            Alphabet::Wass14,
            Alphabet::Hsdm17,
            Alphabet::Uniprot18,
        ]
    }

    /// The moltype for alphabets that have a fixed one. `None` for a seeded random control,
    /// whose moltype carries its seed and so cannot be a `&'static str`.
    fn fixed_moltype(&self) -> Option<&'static str> {
        Some(match self {
            Self::Protein20 => "protein20",
            Self::Dayhoff6 => "dayhoff6",
            Self::HpLehninger2 => "hp_lehninger2",
            Self::HpThomasDill2 => "hp_thomas_dill2",
            Self::HpKyteDoolittle2 => "hp_kyte_doolittle2",
            Self::HpThomasDillNoC2 => "hp_thomas_dill_no_c2",
            Self::HpLehningerCNonpolar2 => "hp_lehninger_c_nonpolar2",
            Self::HpLehningerHpc3 => "hp_lehninger_hpc3",
            Self::HpPBotC1stEd2 => "hp_pbotc_1st_ed2",
            Self::HpRandomControl2 => "hp_random_control2",
            Self::Gbmr4 => "gbmr4",
            Self::Wwmj5 => "wwmj5",
            Self::Gbmr7 => "gbmr7",
            Self::Sdm12 => "sdm12",
            Self::Mmseqs12 => "mmseqs12",
            Self::Wass14 => "wass14",
            Self::Hsdm17 => "hsdm17",
            Self::Uniprot18 => "uniprot18",
            Self::HpRandomControl2Seed(_) => return None,
        })
    }

    /// The moltype stored in an index, e.g. `"sdm12"` or `"hp_thomas_dill2"`.
    ///
    /// A seeded control puts its seed after the class count, so `hp_random_control2_3` is
    /// seed 3 of a 2-class control rather than a 23-class alphabet.
    pub fn to_moltype(&self) -> String {
        match self {
            Self::HpRandomControl2Seed(seed) => format!("hp_random_control2_{seed}"),
            fixed => fixed.fixed_moltype().expect("only the seeded control lacks one").to_string(),
        }
    }

    /// Parse a moltype. Returns `None` for anything unrecognized.
    ///
    /// Takes kmerseek names only. sourmash's `protein`, `dayhoff` and `hp` arrive here
    /// through [`canonical_moltype`], which every entry point applies first.
    pub fn from_moltype(moltype: &str) -> Option<Alphabet> {
        if let Some(alphabet) = Self::all().iter().find(|a| a.fixed_moltype() == Some(moltype)) {
            return Some(*alphabet);
        }
        let seed = moltype.strip_prefix("hp_random_control2_")?.parse().ok()?;
        Some(Self::HpRandomControl2Seed(seed))
    }

    /// Number of classes this alphabet collapses the 20 residues into, which is the number
    /// its moltype ends in.
    pub fn size(&self) -> usize {
        match self {
            Self::Protein20 => 20,
            Self::Dayhoff6 => 6,
            Self::HpLehningerHpc3 => 3,
            Self::Gbmr4 => 4,
            Self::Wwmj5 => 5,
            Self::Gbmr7 => 7,
            Self::Sdm12 | Self::Mmseqs12 => 12,
            Self::Wass14 => 14,
            Self::Hsdm17 => 17,
            Self::Uniprot18 => 18,
            _ => 2,
        }
    }

    /// Whether sourmash encodes this alphabet rather than kmerseek pre-encoding it.
    ///
    /// True for `protein20` (no reduction), `dayhoff6`, and `hp_lehninger2`, which is
    /// sourmash's own `aa_to_hp` partition. Sketches for these are byte-identical to
    /// sourmash's, which is what makes sourmash's own moltypes readable. Everything else
    /// has no sourmash equivalent and goes through [`Self::partition`].
    pub fn uses_sourmash_encoder(&self) -> bool {
        matches!(self, Self::Protein20 | Self::Dayhoff6 | Self::HpLehninger2)
    }

    /// The residue-to-class table. `None` for `protein20`, which does not reduce, and for
    /// `dayhoff6`, whose table lives in sourmash.
    ///
    /// This is the partition itself, so `hp_lehninger2` has one even though sourmash does
    /// its encoding. Use [`alphabet_table`] for the table to pre-encode with.
    pub fn partition(&self) -> Option<&'static HashMap<u8, u8>> {
        Some(match self {
            Self::Protein20 | Self::Dayhoff6 => return None,
            Self::HpLehninger2 => &LEHNINGER_HP,
            Self::HpThomasDill2 => &THOMAS_DILL_HP,
            Self::HpKyteDoolittle2 => &KYTE_DOOLITTLE_HP,
            Self::HpThomasDillNoC2 => &THOMAS_DILL_NO_C_HP,
            Self::HpLehningerCNonpolar2 => &LEHNINGER_C_NONPOLAR_HP,
            Self::HpLehningerHpc3 => &LEHNINGER_HPC,
            Self::HpPBotC1stEd2 => &PBOTC_1ST_ED_HP,
            Self::HpRandomControl2 => &RANDOM_CONTROL_HP,
            Self::HpRandomControl2Seed(1) => &RANDOM_HP_1,
            Self::HpRandomControl2Seed(2) => &RANDOM_HP_2,
            Self::HpRandomControl2Seed(3) => &RANDOM_HP_3,
            Self::HpRandomControl2Seed(4) => &RANDOM_HP_4,
            Self::HpRandomControl2Seed(5) => &RANDOM_HP_5,
            Self::HpRandomControl2Seed(6) => &RANDOM_HP_6,
            Self::HpRandomControl2Seed(7) => &RANDOM_HP_7,
            Self::HpRandomControl2Seed(8) => &RANDOM_HP_8,
            Self::HpRandomControl2Seed(9) => &RANDOM_HP_9,
            Self::HpRandomControl2Seed(10) => &RANDOM_HP_10,
            Self::HpRandomControl2Seed(seed) => {
                panic!("random-control seed {seed} not pre-computed (only 1-10 supported)")
            }
            Self::Gbmr4 => &GBMR4,
            Self::Wwmj5 => &WWMJ5,
            Self::Gbmr7 => &GBMR7,
            Self::Sdm12 => &SDM12,
            Self::Mmseqs12 => &MMSEQS12,
            Self::Wass14 => &WASS14,
            Self::Hsdm17 => &HSDM17,
            Self::Uniprot18 => &UNIPROT18,
        })
    }

    /// The residue clusters in the order the source paper lists them. `None` for the HP
    /// family and the two sourmash-encoded alphabets, which are not defined that way.
    pub fn clusters(&self) -> Option<&'static [&'static str]> {
        Some(match self {
            Self::Gbmr4 => GBMR4_CLUSTERS,
            Self::Wwmj5 => WWMJ5_CLUSTERS,
            Self::Gbmr7 => GBMR7_CLUSTERS,
            Self::Sdm12 => SDM12_CLUSTERS,
            Self::Mmseqs12 => MMSEQS12_CLUSTERS,
            Self::Wass14 => WASS14_CLUSTERS,
            Self::Hsdm17 => HSDM17_CLUSTERS,
            Self::Uniprot18 => UNIPROT18_CLUSTERS,
            _ => return None,
        })
    }

    /// The alphabets defined by residue clusters, meaning the multi-letter ones.
    pub fn multi_letter() -> impl Iterator<Item = &'static Alphabet> {
        Self::all().iter().filter(|a| a.clusters().is_some())
    }

    /// The HP family: two classes, or three where cysteine is split out.
    pub fn hp_family() -> impl Iterator<Item = &'static Alphabet> {
        Self::all()
            .iter()
            .filter(|a| a.fixed_moltype().is_some_and(|moltype| moltype.starts_with("hp_")))
    }
}

fn build_hp(h_residues: &[u8], p_residues: &[u8]) -> HashMap<u8, u8> {
    assert_eq!(
        h_residues.len() + p_residues.len(),
        20,
        "HP table must cover all 20 canonical amino acids"
    );
    let mut m = HashMap::with_capacity(21);
    for &r in h_residues {
        m.insert(r, b'h');
    }
    for &r in p_residues {
        m.insert(r, b'p');
    }
    m.insert(b'*', b'*'); // stop codon pass-through
    m
}

fn build_hpc(h_residues: &[u8], p_residues: &[u8], c_residues: &[u8]) -> HashMap<u8, u8> {
    assert_eq!(
        h_residues.len() + p_residues.len() + c_residues.len(),
        20,
        "HPC table must cover all 20 canonical amino acids"
    );
    let mut m = HashMap::with_capacity(21);
    for &r in h_residues {
        m.insert(r, b'h');
    }
    for &r in p_residues {
        m.insert(r, b'p');
    }
    for &r in c_residues {
        m.insert(r, b'c');
    }
    m.insert(b'*', b'*'); // stop codon pass-through
    m
}

// -----------------------------------------------------------------------------
// Lehninger five-group collapse  (current Kmerseek default as of 2026-04).
//
//   h: A F G I L M P V W Y       (10 residues)
//        = nonpolar aliphatic {G, A, V, L, I, M, P}
//        + aromatic            {F, W, Y}
//   p: C D E H K N Q R S T       (10 residues)
//        = polar uncharged     {S, T, C, N, Q}
//        + positively charged  {K, R, H}
//        + negatively charged  {D, E}
//
// Reference:
//   Nelson, D. L. & Cox, M. M. (2021). Lehninger Principles of
//   Biochemistry, 8th ed. W. H. Freeman. Chapter 3.
// -----------------------------------------------------------------------------
static LEHNINGER_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"AFGILMPVWY", b"CDEHKNQRST"));

// -----------------------------------------------------------------------------
// Thomas-Dill ENERGI 2-class reduction.
// Also adopted by PBotC 2nd ed, Figure 8.28.
//
//   h: A C F I L M V W Y         (9 residues)
//   p: D E G H K N P Q R S T     (11 residues)
//
// References:
//   Thomas, P. D. & Dill, K. A. (1996). An iterative method for
//   extracting energy-like quantities from protein structures.
//   PNAS 93(21):11628-11633. doi:10.1073/pnas.93.21.11628
//   See Fig. 3 (2-class level, "VILMFWYAC vs EDRKGPSTHQN").
//
//   Phillips, R., Kondev, J., Theriot, J., Garcia, H. (2012).
//   Physical Biology of the Cell, 2nd ed. Garland Science.
//   Figure 8.28 (HP classification) and Figure 8.29 (full hierarchy,
//   adapted from Thomas & Dill 1996).
// -----------------------------------------------------------------------------
static THOMAS_DILL_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"ACFILMVWY", b"DEGHKNPQRST"));

// -----------------------------------------------------------------------------
// Kyte-Doolittle hydropathy scale, binarized at hydropathy > 0.
//
//   h: A C F I L M V             (7 residues; all with hydropathy > 0)
//   p: D E G H K N P Q R S T W Y (13 residues; hydropathy <= 0)
//
// Notably: W (-0.9) and Y (-1.3) are polar on this scale, despite
// being aromatics, because of their H-bond donor groups (indole NH,
// hydroxyl).
//
// Reference:
//   Kyte, J. & Doolittle, R. F. (1982). A simple method for displaying
//   the hydropathic character of a protein. J Mol Biol 157(1):105-132.
//   doi:10.1016/0022-2836(82)90515-0
// -----------------------------------------------------------------------------
static KYTE_DOOLITTLE_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"ACFILMV", b"DEGHKNPQRSTWY"));

// -----------------------------------------------------------------------------
// Thomas-Dill with cysteine reassigned to polar.
// Isolation variant: tests how much of any TD-vs-Lehninger performance
// gap is attributable to C placement specifically.
//
//   h: A F I L M V W Y           (8 residues)
//   p: C D E G H K N P Q R S T   (12 residues)
// -----------------------------------------------------------------------------
static THOMAS_DILL_NO_C_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"AFILMVWY", b"CDEGHKNPQRST"));

// -----------------------------------------------------------------------------
// Lehninger with cysteine added to H.
// Isolation variant: keeps Lehninger's G/P/Y placement but adopts the
// Kyte-Doolittle/Thomas-Dill treatment of C.
//
//   h: A C F G I L M P V W Y     (11 residues)
//   p: D E H K N Q R S T         (9 residues)
// -----------------------------------------------------------------------------
static LEHNINGER_C_NONPOLAR_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"ACFGILMPVWY", b"DEHKNQRST"));

// -----------------------------------------------------------------------------
// Lehninger HPC: 3-letter extension of LehningerCNonpolar.
//
// LehningerCNonpolar folds cysteine into the hydrophobic class because its thiol
// side chain is nonpolar. But cysteine's ability to form disulfide bonds is a
// distinct chemistry from ordinary hydrophobic packing, so this variant keeps
// Lehninger's H/P split for the other 19 residues and gives cysteine its own
// third symbol ('c', for cystine) instead of merging it into 'h'.
//
//   h: A F G I L M P V W Y       (10 residues; = Lehninger's h, C excluded)
//   p: D E H K N Q R S T         (9 residues; = Lehninger's p, C excluded)
//   c: C                         (1 residue)
// -----------------------------------------------------------------------------
static LEHNINGER_HPC: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hpc(b"AFGILMPVWY", b"DEHKNQRST", b"C"));

// -----------------------------------------------------------------------------
// Phillips et al. PBotC 1st ed, Figure 8.30.
// Cys and Tyr are drawn in parentheses in the original figure,
// indicating borderline H status. We follow the figure's grouping and
// place them in h. Pro is firmly h in this edition (moved to p in the
// 2nd ed, which adopted Thomas-Dill).
//
//   h: A C F I L M P V W Y       (10 residues)
//   p: D E G H K N Q R S T       (10 residues)
//
// Reference:
//   Phillips, R., Kondev, J., Theriot, J. (2008). Physical Biology of
//   the Cell, 1st ed. Garland Science. Figure 8.30.
// -----------------------------------------------------------------------------
static PBOTC_1ST_ED_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"ACFILMPVWY", b"DEGHKNQRST"));

// -----------------------------------------------------------------------------
// Random negative control.
//
// Partition generated by shuffling the 20 canonical amino acids and
// splitting at the midpoint. It scrambles the hydrophobicity signal
// while maintaining a 10/10 split. By design, it mixes strongly
// hydrophobic residues (I, V, F) into p and polar/charged residues
// (D, K, R, Q) into h.
//
//   h: A D G K L M Q R W Y       (10 residues)
//   p: C E F H I N P S T V       (10 residues)
//
// For the supplementary figure, consider running multiple independent
// shuffles (e.g. 10 different seeds) and reporting the distribution of
// AUCs rather than a single point — this makes the negative control
// statistically robust. See `random_hp()` below.
// -----------------------------------------------------------------------------
static RANDOM_CONTROL_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"ADGKLMQRWY", b"CEFHINPSTV"));

// Seeded random controls — 10 independent random partitions for null-distribution estimation.
static RANDOM_HP_1: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(1));
static RANDOM_HP_2: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(2));
static RANDOM_HP_3: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(3));
static RANDOM_HP_4: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(4));
static RANDOM_HP_5: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(5));
static RANDOM_HP_6: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(6));
static RANDOM_HP_7: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(7));
static RANDOM_HP_8: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(8));
static RANDOM_HP_9: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(9));
static RANDOM_HP_10: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| random_hp(10));

/// Generate a seeded random HP partition for negative-control runs.
/// Pass multiple seeds to characterize the null distribution.
pub fn random_hp(seed: u64) -> HashMap<u8, u8> {
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    let mut residues: Vec<u8> = b"ACDEFGHIKLMNPQRSTVWY".to_vec();
    let mut rng = StdRng::seed_from_u64(seed);
    residues.shuffle(&mut rng);
    build_hp(&residues[..10], &residues[10..])
}

/// The partition of an alphabet that has one, for tests that compare tables directly.
#[cfg(test)]
fn test_partition(alphabet: Alphabet) -> &'static HashMap<u8, u8> {
    alphabet.partition().expect("these alphabets all reduce")
}

#[cfg(test)]
mod hp_tests {
    use super::*;

    fn all_named_alphabets() -> Vec<Alphabet> {
        Alphabet::hp_family().copied().collect()
    }

    #[test]
    fn all_alphabets_cover_20_residues() {
        for alpha in all_named_alphabets() {
            let t = test_partition(alpha);
            // 20 amino acids + 1 stop codon
            assert_eq!(t.len(), 21, "alphabet {:?} has wrong table size", alpha);
            for &r in CANONICAL_AA {
                assert!(t.contains_key(&r), "alphabet {:?} missing residue {}", alpha, r as char);
            }
        }
    }

    #[test]
    fn all_mappings_are_h_p_or_c() {
        for alpha in all_named_alphabets() {
            for &r in CANONICAL_AA {
                let encoded = test_partition(alpha)[&r];
                assert!(
                    encoded == b'h' || encoded == b'p' || encoded == b'c',
                    "alphabet {:?} residue {} mapped to unexpected byte {}",
                    alpha,
                    r as char,
                    encoded as char
                );
            }
        }
    }

    #[test]
    fn stop_codon_passes_through() {
        for alpha in all_named_alphabets() {
            assert_eq!(
                test_partition(alpha)[&b'*'],
                b'*',
                "alphabet {:?} does not pass through stop codon",
                alpha
            );
        }
    }

    // Controversial residues: C, G, P, W, Y
    // Source: table in kmerseek paper supplementary (comparing all six alphabets).
    // Each assertion cites the reference that places the residue in that class.

    #[test]
    fn cysteine_placement() {
        // C = polar in Lehninger (polar uncharged group) and ThomasDillNoC (by construction).
        assert_eq!(
            test_partition(Alphabet::HpLehninger2)[&b'C'],
            b'p',
            "Lehninger: C polar (Nelson & Cox 2021, Ch. 3)"
        );
        assert_eq!(
            test_partition(Alphabet::HpThomasDillNoC2)[&b'C'],
            b'p',
            "ThomasDillNoC: C polar by construction"
        );

        // C = hydrophobic in ThomasDill, KyteDoolittle, LehningerCNonpolar, PBotC1stEd.
        assert_eq!(
            test_partition(Alphabet::HpThomasDill2)[&b'C'],
            b'h',
            "ThomasDill: C hydrophobic (Thomas & Dill 1996, Fig. 3)"
        );
        assert_eq!(
            test_partition(Alphabet::HpKyteDoolittle2)[&b'C'],
            b'h',
            "KyteDoolittle: C hydrophobic (hydropathy +2.5, Kyte & Doolittle 1982)"
        );
        assert_eq!(
            test_partition(Alphabet::HpLehningerCNonpolar2)[&b'C'],
            b'h',
            "LehningerCNonpolar: C hydrophobic by construction"
        );
        assert_eq!(
            test_partition(Alphabet::HpPBotC1stEd2)[&b'C'],
            b'h',
            "PBotC1stEd: C borderline-h (Phillips et al. 2008, Fig. 8.30)"
        );

        // C = its own third class (cystine) in LehningerHpc, distinct from h/p.
        assert_eq!(
            test_partition(Alphabet::HpLehningerHpc3)[&b'C'],
            b'c',
            "LehningerHpc: C is cystine, its own class distinct from h/p"
        );
    }

    #[test]
    fn glycine_placement() {
        // G = hydrophobic in Lehninger (nonpolar aliphatic) and LehningerCNonpolar.
        assert_eq!(
            test_partition(Alphabet::HpLehninger2)[&b'G'],
            b'h',
            "Lehninger: G nonpolar aliphatic (Nelson & Cox 2021, Ch. 3)"
        );
        assert_eq!(
            test_partition(Alphabet::HpLehningerCNonpolar2)[&b'G'],
            b'h',
            "LehningerCNonpolar: G inherits Lehninger placement"
        );
        assert_eq!(
            test_partition(Alphabet::HpLehningerHpc3)[&b'G'],
            b'h',
            "LehningerHpc: G inherits Lehninger placement"
        );

        // G = polar in ThomasDill, KyteDoolittle, ThomasDillNoC, PBotC1stEd.
        assert_eq!(
            test_partition(Alphabet::HpThomasDill2)[&b'G'],
            b'p',
            "ThomasDill: G polar (Thomas & Dill 1996, Fig. 3)"
        );
        assert_eq!(
            test_partition(Alphabet::HpKyteDoolittle2)[&b'G'],
            b'p',
            "KyteDoolittle: G polar (hydropathy -0.4, Kyte & Doolittle 1982)"
        );
        assert_eq!(
            test_partition(Alphabet::HpThomasDillNoC2)[&b'G'],
            b'p',
            "ThomasDillNoC: G polar (inherits ThomasDill)"
        );
        assert_eq!(
            test_partition(Alphabet::HpPBotC1stEd2)[&b'G'],
            b'p',
            "PBotC1stEd: G polar (Phillips et al. 2008, Fig. 8.30)"
        );
    }

    #[test]
    fn proline_placement() {
        // P = hydrophobic in Lehninger (nonpolar aliphatic), LehningerCNonpolar, PBotC1stEd.
        assert_eq!(
            test_partition(Alphabet::HpLehninger2)[&b'P'],
            b'h',
            "Lehninger: P nonpolar aliphatic (Nelson & Cox 2021, Ch. 3)"
        );
        assert_eq!(
            test_partition(Alphabet::HpLehningerCNonpolar2)[&b'P'],
            b'h',
            "LehningerCNonpolar: P inherits Lehninger placement"
        );
        assert_eq!(
            test_partition(Alphabet::HpPBotC1stEd2)[&b'P'],
            b'h',
            "PBotC1stEd: P hydrophobic (Phillips et al. 2008, Fig. 8.30; moved to p in 2nd ed)"
        );
        assert_eq!(
            test_partition(Alphabet::HpLehningerHpc3)[&b'P'],
            b'h',
            "LehningerHpc: P inherits Lehninger placement"
        );

        // P = polar in ThomasDill, KyteDoolittle, ThomasDillNoC.
        assert_eq!(
            test_partition(Alphabet::HpThomasDill2)[&b'P'],
            b'p',
            "ThomasDill: P polar (Thomas & Dill 1996, Fig. 3)"
        );
        assert_eq!(
            test_partition(Alphabet::HpKyteDoolittle2)[&b'P'],
            b'p',
            "KyteDoolittle: P polar (hydropathy -1.6, Kyte & Doolittle 1982)"
        );
        assert_eq!(
            test_partition(Alphabet::HpThomasDillNoC2)[&b'P'],
            b'p',
            "ThomasDillNoC: P polar (inherits ThomasDill)"
        );
    }

    #[test]
    fn tryptophan_placement() {
        // W = polar in KyteDoolittle only (hydropathy -0.9 due to indole NH).
        assert_eq!(
            test_partition(Alphabet::HpKyteDoolittle2)[&b'W'],
            b'p',
            "KyteDoolittle: W polar (hydropathy -0.9, Kyte & Doolittle 1982)"
        );

        // W = hydrophobic in all other named alphabets.
        for alpha in [
            Alphabet::HpLehninger2,
            Alphabet::HpThomasDill2,
            Alphabet::HpThomasDillNoC2,
            Alphabet::HpLehningerCNonpolar2,
            Alphabet::HpLehningerHpc3,
            Alphabet::HpPBotC1stEd2,
        ] {
            assert_eq!(test_partition(alpha)[&b'W'], b'h', "{:?}: W hydrophobic (aromatic)", alpha);
        }
    }

    #[test]
    fn tyrosine_placement() {
        // Y = polar in KyteDoolittle only (hydropathy -1.3 due to hydroxyl).
        assert_eq!(
            test_partition(Alphabet::HpKyteDoolittle2)[&b'Y'],
            b'p',
            "KyteDoolittle: Y polar (hydropathy -1.3, Kyte & Doolittle 1982)"
        );

        // Y = hydrophobic in all other named alphabets.
        for alpha in [
            Alphabet::HpLehninger2,
            Alphabet::HpThomasDill2,
            Alphabet::HpThomasDillNoC2,
            Alphabet::HpLehningerCNonpolar2,
            Alphabet::HpLehningerHpc3,
            Alphabet::HpPBotC1stEd2,
        ] {
            assert_eq!(test_partition(alpha)[&b'Y'], b'h', "{:?}: Y hydrophobic (aromatic)", alpha);
        }
    }

    #[test]
    fn random_hp_deterministic() {
        let a = random_hp(42);
        let b = random_hp(42);
        assert_eq!(a, b, "random_hp must be deterministic for the same seed");
    }

    #[test]
    fn random_hp_covers_all_residues() {
        let m = random_hp(0);
        assert_eq!(m.len(), 21);
        for &r in CANONICAL_AA {
            assert!(m.contains_key(&r), "random_hp missing residue {}", r as char);
            let v = m[&r];
            assert!(v == b'h' || v == b'p');
        }
    }

    #[test]
    fn random_hp_differs_by_seed() {
        let a = random_hp(1);
        let b = random_hp(2);
        // Two different seeds should produce different partitions
        // (astronomically unlikely to collide with 20 elements)
        assert_ne!(a, b, "random_hp with different seeds should differ");
    }

    #[test]
    fn random_control_static_matches_documented_partition() {
        // h: ADGKLMQRWY  p: CEFHINPSTV — frozen in static for reproducibility.
        let t = test_partition(Alphabet::HpRandomControl2);
        for &r in b"ADGKLMQRWY" {
            assert_eq!(t[&r], b'h', "RandomControl: {} should be h", r as char);
        }
        for &r in b"CEFHINPSTV" {
            assert_eq!(t[&r], b'p', "RandomControl: {} should be p", r as char);
        }
    }

    #[test]
    fn all_alphabet_names_are_unique() {
        let names: Vec<_> = Alphabet::hp_family().map(|a| a.to_moltype()).collect();
        let mut sorted = names.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(names.len(), sorted.len(), "duplicate alphabet names detected");
    }

    /// Every alphabet's moltype must end in its class count, so the HP names read the same
    /// way as `sdm12` and `gbmr4`.
    #[test]
    fn moltype_ends_in_class_count() {
        let expected = [
            (Alphabet::HpLehninger2, "hp_lehninger2"),
            (Alphabet::HpThomasDill2, "hp_thomas_dill2"),
            (Alphabet::HpKyteDoolittle2, "hp_kyte_doolittle2"),
            (Alphabet::HpThomasDillNoC2, "hp_thomas_dill_no_c2"),
            (Alphabet::HpLehningerCNonpolar2, "hp_lehninger_c_nonpolar2"),
            (Alphabet::HpLehningerHpc3, "hp_lehninger_hpc3"),
            (Alphabet::HpPBotC1stEd2, "hp_pbotc_1st_ed2"),
            (Alphabet::HpRandomControl2, "hp_random_control2"),
        ];

        for (alphabet, moltype) in expected {
            assert_eq!(alphabet.to_moltype(), moltype);
            assert_eq!(Alphabet::from_moltype(moltype), Some(alphabet));
        }
    }

    /// The count in the name has to be the number of distinct symbols the table emits,
    /// otherwise the name lies about what the index stores.
    #[test]
    fn class_count_matches_distinct_table_symbols() {
        for alphabet in Alphabet::hp_family() {
            let symbols: std::collections::HashSet<u8> = test_partition(*alphabet)
                .iter()
                .filter(|(residue, _)| **residue != b'*')
                .map(|(_, symbol)| *symbol)
                .collect();
            assert_eq!(symbols.len(), alphabet.size(), "{}", alphabet.to_moltype());
            assert!(
                alphabet.to_moltype().ends_with(&alphabet.size().to_string()),
                "{}",
                alphabet.to_moltype()
            );
        }
    }

    /// Seeded controls put the seed after the class count, so seed 3 of a 2-class control
    /// cannot be misread as a 23-class alphabet.
    #[test]
    fn seeded_random_control_moltype_separates_count_from_seed() {
        assert_eq!(Alphabet::HpRandomControl2Seed(3).to_moltype(), "hp_random_control2_3");
        assert_eq!(
            Alphabet::from_moltype("hp_random_control2_3"),
            Some(Alphabet::HpRandomControl2Seed(3))
        );
        assert_eq!(
            Alphabet::from_moltype("hp_random_control2_10"),
            Some(Alphabet::HpRandomControl2Seed(10))
        );
    }

    #[test]
    fn from_moltype_rejects_unknown_strings() {
        // `protein` and `hp` are sourmash spellings; they reach an Alphabet only through
        // canonical_moltype, not from_moltype. The rest are simply not alphabets.
        for moltype in ["protein", "hp", "dayhoff", "raw", "hp_", "hp_nope2", "reduced_sdm12"] {
            assert_eq!(Alphabet::from_moltype(moltype), None, "{moltype}");
        }
    }

    /// The claim the README's table rests on: every scheme agrees on 15 of the 20 residues
    /// and disagrees only about C, G, P, W and Y. A new alphabet that moved one of the other
    /// 15 would make the table wrong.
    #[test]
    fn schemes_differ_only_on_the_five_borderline_residues() {
        const ALWAYS_HYDROPHOBIC: &[u8] = b"AFILMV";
        const ALWAYS_POLAR: &[u8] = b"DEHKNQRST";
        const BORDERLINE: &[u8] = b"CGPWY";

        assert_eq!(
            ALWAYS_HYDROPHOBIC.len() + ALWAYS_POLAR.len() + BORDERLINE.len(),
            20,
            "the three groups must partition the canonical residues"
        );

        for alphabet in Alphabet::hp_family() {
            // The randomized control is a negative control, so it has no reason to agree.
            if matches!(alphabet, Alphabet::HpRandomControl2) {
                continue;
            }
            let table = test_partition(*alphabet);
            for &residue in ALWAYS_HYDROPHOBIC {
                assert_eq!(
                    table[&residue],
                    b'h',
                    "{}: {} should be hydrophobic in every scheme",
                    alphabet.to_moltype(),
                    residue as char
                );
            }
            for &residue in ALWAYS_POLAR {
                assert_eq!(
                    table[&residue],
                    b'p',
                    "{}: {} should be polar in every scheme",
                    alphabet.to_moltype(),
                    residue as char
                );
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Multi-letter alphabets: 4 to 18 classes.
// Multi-letter reduced amino acid alphabets from the fold-recognition literature.
//
// Unlike [`crate::hp_alphabets`], which splits the 20 residues two or three ways,
// these alphabets keep 4 to 18 classes. Peterson et al. (2009) benchmarked over 150
// published clustering schemes against DALI fold assignments and reported GBMR4,
// SDM12 and HSDM17 as the best performers on recall, AUC and mean pooled precision
// respectively. Ieremie et al. (2024) reused those three and added GBMR7, WWMJ5,
// MMSEQS12, WASS14 and UNIPROT18 when testing how alphabet reduction affects protein
// language models. Where the two papers overlap, their partitions agree.
//
// References:
//   Peterson, E. L., Kondev, J., Theriot, J. A. & Phillips, R. (2009). Reduced amino
//   acid alphabets exhibit an improved sensitivity and selectivity in fold assignment.
//   Bioinformatics 25(11):1356-1362. doi:10.1093/bioinformatics/btp164. Table 2.
//
//   Ieremie, I., Ewing, R. M. & Niranjan, M. (2024). Protein language models meet
//   reduced amino acid alphabets. Bioinformatics 40(2):btae061.
//   doi:10.1093/bioinformatics/btae061. Table 1.
// ----------------------------------------------------------------------------

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
// Same maximum-mutual-information construction as GBMR4 at a finer level, but not a
// refinement of it: the large AEFIKLMQRVWY class mixes residues that GBMR4 keeps on
// opposite sides of its hydrophobic/polar split.
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
mod multi_letter_tests {
    use super::*;

    fn encode_canonical(alphabet: Alphabet) -> String {
        let table = test_partition(alphabet);
        CANONICAL_AA.iter().map(|r| *table.get(r).expect("canonical residue") as char).collect()
    }

    #[test]
    fn test_gbmr4_partition() {
        // ADKERNTSQ -> a, YFLIVMCWH -> y, G -> g, P -> p
        assert_eq!(encode_canonical(Alphabet::Gbmr4), "ayaaygyyayyapaaaayyy");
    }

    #[test]
    fn test_wwmj5_partition() {
        // CMFILVWY -> c, ATH -> a, GP -> g, DE -> d, SNQRK -> s
        assert_eq!(encode_canonical(Alphabet::Wwmj5), "acddcgacsccsgsssaccc");
    }

    #[test]
    fn test_gbmr7_partition() {
        // DN -> d, AEFIKLMQRVWY -> a, CH -> c, T -> t, S -> s, G -> g, P -> p
        assert_eq!(encode_canonical(Alphabet::Gbmr7), "acdaagcaaaadpaastaaa");
    }

    #[test]
    fn test_sdm12_partition() {
        // A a | D d | KER k | N n | TSQ t | YF y | LIVM l | C c | W w | H h | G g | P p
        assert_eq!(encode_canonical(Alphabet::Sdm12), "acdkyghlkllnptkttlwy");
    }

    #[test]
    fn test_mmseqs12_partition() {
        // AST a | LM l | IV i | KR k | EQ e | ND n | FY f | C c | G g | H h | P p | W w
        assert_eq!(encode_canonical(Alphabet::Mmseqs12), "acnefghikllnpekaaiwf");
    }

    #[test]
    fn test_wass14_partition() {
        // WM w | DI d | P p | C c | AV a | K k | T t | RE r | G g | L l | Y y | SH s | F f | NQ n
        assert_eq!(encode_canonical(Alphabet::Wass14), "acdrfgsdklwnpnrstawy");
    }

    #[test]
    fn test_hsdm17_partition() {
        // KE k | LIV l | every other residue on its own
        assert_eq!(encode_canonical(Alphabet::Hsdm17), "acdkfghlklmnpqrstlwy");
    }

    #[test]
    fn test_uniprot18_partition() {
        // EP e | HL h | every other residue on its own
        assert_eq!(encode_canonical(Alphabet::Uniprot18), "acdefghikhmneqrstvwy");
    }

    /// `build_reduced` asserts full coverage and distinct symbols, so constructing every
    /// table is itself the check that no cluster list drifts out of shape.
    #[test]
    fn test_every_alphabet_builds_and_has_declared_size() {
        let expected_sizes = [4usize, 5, 7, 12, 12, 14, 17, 18];
        for (alphabet, expected) in Alphabet::multi_letter().zip(expected_sizes) {
            let table = test_partition(*alphabet);
            // 20 canonical residues plus the stop-codon pass-through.
            assert_eq!(table.len(), 21, "{}", alphabet.to_moltype());
            assert_eq!(table.get(&b'*'), Some(&b'*'), "{}", alphabet.to_moltype());
            assert_eq!(alphabet.size(), expected, "{}", alphabet.to_moltype());
        }
    }

    /// The numeric suffix in the name has to be the class count, since that is how both
    /// papers and every downstream filename refer to these alphabets.
    #[test]
    fn test_name_suffix_matches_class_count() {
        for alphabet in Alphabet::multi_letter() {
            let digits: String =
                alphabet.to_moltype().chars().filter(|c| c.is_ascii_digit()).collect();
            assert_eq!(
                digits.parse::<usize>().unwrap(),
                alphabet.size(),
                "{}",
                alphabet.to_moltype()
            );
        }
    }

    /// The moltype is the paper's own name, digit included.
    #[test]
    fn test_moltype_round_trip() {
        for alphabet in Alphabet::multi_letter() {
            let moltype = alphabet.to_moltype();
            assert_eq!(moltype, alphabet.to_moltype());
            assert_eq!(Alphabet::from_moltype(&moltype), Some(*alphabet));
        }
        assert_eq!(Alphabet::from_moltype("sdm12"), Some(Alphabet::Sdm12));
        assert_eq!(Alphabet::from_moltype("hsdm17"), Some(Alphabet::Hsdm17));
    }

    /// GBMR4, SDM12 and HSDM17 appear in both source papers, and Peterson et al. note
    /// that each is a refinement of the previous one: splitting classes only, never
    /// moving a residue across an existing boundary. Two residues sharing a class in the
    /// finer alphabet must therefore share one in the coarser alphabet too.
    #[test]
    fn test_hsdm17_refines_sdm12_refines_gbmr4() {
        let chain = [(Alphabet::Hsdm17, Alphabet::Sdm12), (Alphabet::Sdm12, Alphabet::Gbmr4)];
        for (finer, coarser) in chain {
            let (fine, coarse) = (test_partition(finer), test_partition(coarser));
            for &a in CANONICAL_AA {
                for &b in CANONICAL_AA {
                    if fine[&a] == fine[&b] {
                        assert_eq!(
                            coarse[&a],
                            coarse[&b],
                            "{} groups {} with {} but {} splits them",
                            finer.to_moltype(),
                            a as char,
                            b as char,
                            coarser.to_moltype()
                        );
                    }
                }
            }
        }
    }
}

/// The table to pre-encode a sequence with before hashing.
///
/// `None` for the three alphabets sourmash encodes itself, which [`crate::hash_functions`]
/// routes to sourmash instead.
pub fn alphabet_table(moltype: &str) -> Option<&'static HashMap<u8, u8>> {
    let alphabet = Alphabet::from_moltype(canonical_moltype(moltype))?;
    if alphabet.uses_sourmash_encoder() {
        return None;
    }
    alphabet.partition()
}

#[cfg(test)]
mod sourmash_compat_tests {
    use super::*;

    /// sourmash's three moltypes name partitions kmerseek also has, and kmerseek hashes all
    /// three through sourmash's own hash function, so data carrying these names is readable.
    #[test]
    fn sourmash_moltypes_map_to_the_same_alphabet() {
        assert_eq!(canonical_moltype("protein"), "protein20");
        assert_eq!(canonical_moltype("dayhoff"), "dayhoff6");
        assert_eq!(canonical_moltype("hp"), "hp_lehninger2");
    }

    /// kmerseek's own names pass through untouched.
    #[test]
    fn kmerseek_names_pass_through() {
        for moltype in
            ["protein20", "dayhoff6", "hp_lehninger2", "hp_thomas_dill2", "sdm12", "gbmr4"]
        {
            assert_eq!(canonical_moltype(moltype), moltype, "{moltype}");
        }
    }

    /// Spellings that were never sourmash's get no compatibility path; they are rejected
    /// further along rather than quietly mapped.
    #[test]
    fn earlier_kmerseek_spellings_are_not_translated() {
        for moltype in ["raw", "hp_lehninger", "hp_thomas_dill", "reduced_sdm12"] {
            assert_eq!(canonical_moltype(moltype), moltype, "{moltype}");
            assert!(Alphabet::from_moltype(moltype).is_none(), "{moltype}");
            assert!(Alphabet::from_moltype(moltype).is_none(), "{moltype}");
        }
    }
}
