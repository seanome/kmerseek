//! Every reduced amino acid alphabet kmerseek can index with.
//!
//! Two families live here. The HP tables split the 20 residues into hydrophobic and polar,
//! differing only on the borderline residues C, G, P, W and Y, plus one three-class variant
//! that gives cysteine its own symbol. The multi-letter alphabets keep 4 to 18 classes and
//! come from the fold-recognition literature.
//!
//! An alphabet is a partition of the 20 canonical residues, stored as a residue-to-class
//! table. [`alphabet_table`] looks one up by moltype. It returns `None` for the three
//! alphabets sourmash encodes for us, which [`crate::hash_functions`] handles instead.

use std::collections::HashMap;
use std::sync::LazyLock;

// ----------------------------------------------------------------------------
// HP alphabets: two classes, or three where cysteine is split out.
// HP alphabet mappings for the Kmerseek robustness sweep.
//
// Each alphabet partitions the 20 canonical amino acids into
// hydrophobic (h) and polar (p) classes. The alphabets differ on
// borderline residues — principally C, G, P, W, and Y.
// ----------------------------------------------------------------------------

/// HP alphabet variants evaluated in the robustness sweep.
///
/// See module docs for motivation. The default production alphabet is
/// determined by the sweep results (see Supplementary Figure N of the
/// Kmerseek paper).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HpAlphabet {
    Lehninger,
    ThomasDill,
    KyteDoolittle,
    ThomasDillNoC,
    LehningerCNonpolar,
    LehningerHpc,
    PBotC1stEd,
    RandomControl,
    /// Seeded random control; seed must be 1-10.
    Random(u64),
}

/// Translate a sourmash moltype into the kmerseek name for the same alphabet.
///
/// sourmash writes `protein`, `dayhoff` and `hp`. Each names a partition kmerseek also has,
/// and kmerseek hashes all three through sourmash's own hash function rather than
/// pre-encoding them, so a sketch or index carrying one of these names holds hashes
/// kmerseek can read directly. Reading the names keeps that data usable.
///
/// Every other string passes through unchanged, including kmerseek's own earlier spellings
/// (`raw`, `hp_<name>` without a class count). Those were never sourmash's and have no
/// interop argument, so they are rejected further along.
pub fn canonical_moltype(moltype: &str) -> &str {
    match moltype {
        "protein" => "protein20",
        "dayhoff" => "dayhoff6",
        "hp" => "hp_lehninger2",
        other => other,
    }
}

impl HpAlphabet {
    pub fn table(&self) -> &'static HashMap<u8, u8> {
        match self {
            Self::Lehninger => &LEHNINGER_HP,
            Self::ThomasDill => &THOMAS_DILL_HP,
            Self::KyteDoolittle => &KYTE_DOOLITTLE_HP,
            Self::ThomasDillNoC => &THOMAS_DILL_NO_C_HP,
            Self::LehningerCNonpolar => &LEHNINGER_C_NONPOLAR_HP,
            Self::LehningerHpc => &LEHNINGER_HPC,
            Self::PBotC1stEd => &PBOTC_1ST_ED_HP,
            Self::RandomControl => &RANDOM_CONTROL_HP,
            Self::Random(1) => &RANDOM_HP_1,
            Self::Random(2) => &RANDOM_HP_2,
            Self::Random(3) => &RANDOM_HP_3,
            Self::Random(4) => &RANDOM_HP_4,
            Self::Random(5) => &RANDOM_HP_5,
            Self::Random(6) => &RANDOM_HP_6,
            Self::Random(7) => &RANDOM_HP_7,
            Self::Random(8) => &RANDOM_HP_8,
            Self::Random(9) => &RANDOM_HP_9,
            Self::Random(10) => &RANDOM_HP_10,
            Self::Random(n) => {
                panic!("random-control seed {n} not pre-computed (only 1-10 supported)")
            }
        }
    }

    /// Short identifier used in output filenames and Nextflow traces.
    pub fn name(&self) -> &'static str {
        match self {
            Self::Lehninger => "lehninger",
            Self::ThomasDill => "thomas_dill",
            Self::KyteDoolittle => "kyte_doolittle",
            Self::ThomasDillNoC => "thomas_dill_no_c",
            Self::LehningerCNonpolar => "lehninger_c_nonpolar",
            Self::LehningerHpc => "lehninger_hpc",
            Self::PBotC1stEd => "pbotc_1st_ed",
            Self::RandomControl => "random_control",
            Self::Random(1) => "random_control_1",
            Self::Random(2) => "random_control_2",
            Self::Random(3) => "random_control_3",
            Self::Random(4) => "random_control_4",
            Self::Random(5) => "random_control_5",
            Self::Random(6) => "random_control_6",
            Self::Random(7) => "random_control_7",
            Self::Random(8) => "random_control_8",
            Self::Random(9) => "random_control_9",
            Self::Random(10) => "random_control_10",
            Self::Random(n) => {
                panic!("random-control seed {n} not pre-computed (only 1-10 supported)")
            }
        }
    }

    pub fn all_named() -> &'static [HpAlphabet] {
        &[
            HpAlphabet::Lehninger,
            HpAlphabet::ThomasDill,
            HpAlphabet::KyteDoolittle,
            HpAlphabet::ThomasDillNoC,
            HpAlphabet::LehningerCNonpolar,
            HpAlphabet::LehningerHpc,
            HpAlphabet::PBotC1stEd,
            HpAlphabet::RandomControl,
        ]
    }

    /// Number of classes this alphabet partitions the 20 residues into.
    pub fn size(&self) -> usize {
        match self {
            Self::LehningerHpc => 3,
            _ => 2,
        }
    }

    /// Moltype string stored in the index (e.g. `"hp_thomas_dill2"`).
    ///
    /// The `hp_` prefix keeps the family greppable; the trailing digit is the class count,
    /// so these read the same way as the multi-letter alphabets in
    /// the multi-letter alphabets below (`sdm12`, `gbmr4`). Every alphabet states its size.
    pub fn to_moltype(&self) -> String {
        match self {
            // The seed goes after the class count, so "…control2_3" is seed 3 of a
            // 2-class control and never reads as class count 23.
            Self::Random(n) => format!("hp_random_control2_{n}"),
            _ => format!("hp_{}{}", self.name(), self.size()),
        }
    }

    /// Parse from a moltype string. Returns `None` for unrecognized strings.
    ///
    /// Takes kmerseek names only. sourmash's `hp` reaches [`Self::Lehninger`] through
    /// [`canonical_moltype`], which every entry point applies first.
    pub fn from_moltype(s: &str) -> Option<HpAlphabet> {
        if let Some(alphabet) = Self::all_named().iter().find(|a| a.to_moltype() == s) {
            return Some(*alphabet);
        }
        let seed = s.strip_prefix("hp_random_control2_")?.parse().ok()?;
        Some(HpAlphabet::Random(seed))
    }

    /// Whether sourmash encodes this alphabet itself.
    ///
    /// Lehninger is sourmash's own `aa_to_hp` partition, so `hp_lehninger2` is sketched
    /// through sourmash's Murmur64Hp rather than pre-encoded here. That keeps its hashes
    /// identical to sourmash's for the same sequence. The other HP tables have no sourmash
    /// equivalent and are pre-encoded.
    pub fn uses_sourmash_encoder(&self) -> bool {
        matches!(self, Self::Lehninger)
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

#[cfg(test)]
mod hp_tests {
    use super::*;

    const ALL_RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";

    fn all_named_alphabets() -> &'static [HpAlphabet] {
        HpAlphabet::all_named()
    }

    #[test]
    fn all_alphabets_cover_20_residues() {
        for alpha in all_named_alphabets() {
            let t = alpha.table();
            // 20 amino acids + 1 stop codon
            assert_eq!(t.len(), 21, "alphabet {:?} has wrong table size", alpha);
            for &r in ALL_RESIDUES {
                assert!(t.contains_key(&r), "alphabet {:?} missing residue {}", alpha, r as char);
            }
        }
    }

    #[test]
    fn all_mappings_are_h_p_or_c() {
        for alpha in all_named_alphabets() {
            for &r in ALL_RESIDUES {
                let encoded = alpha.table()[&r];
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
                alpha.table()[&b'*'],
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
            HpAlphabet::Lehninger.table()[&b'C'],
            b'p',
            "Lehninger: C polar (Nelson & Cox 2021, Ch. 3)"
        );
        assert_eq!(
            HpAlphabet::ThomasDillNoC.table()[&b'C'],
            b'p',
            "ThomasDillNoC: C polar by construction"
        );

        // C = hydrophobic in ThomasDill, KyteDoolittle, LehningerCNonpolar, PBotC1stEd.
        assert_eq!(
            HpAlphabet::ThomasDill.table()[&b'C'],
            b'h',
            "ThomasDill: C hydrophobic (Thomas & Dill 1996, Fig. 3)"
        );
        assert_eq!(
            HpAlphabet::KyteDoolittle.table()[&b'C'],
            b'h',
            "KyteDoolittle: C hydrophobic (hydropathy +2.5, Kyte & Doolittle 1982)"
        );
        assert_eq!(
            HpAlphabet::LehningerCNonpolar.table()[&b'C'],
            b'h',
            "LehningerCNonpolar: C hydrophobic by construction"
        );
        assert_eq!(
            HpAlphabet::PBotC1stEd.table()[&b'C'],
            b'h',
            "PBotC1stEd: C borderline-h (Phillips et al. 2008, Fig. 8.30)"
        );

        // C = its own third class (cystine) in LehningerHpc, distinct from h/p.
        assert_eq!(
            HpAlphabet::LehningerHpc.table()[&b'C'],
            b'c',
            "LehningerHpc: C is cystine, its own class distinct from h/p"
        );
    }

    #[test]
    fn glycine_placement() {
        // G = hydrophobic in Lehninger (nonpolar aliphatic) and LehningerCNonpolar.
        assert_eq!(
            HpAlphabet::Lehninger.table()[&b'G'],
            b'h',
            "Lehninger: G nonpolar aliphatic (Nelson & Cox 2021, Ch. 3)"
        );
        assert_eq!(
            HpAlphabet::LehningerCNonpolar.table()[&b'G'],
            b'h',
            "LehningerCNonpolar: G inherits Lehninger placement"
        );
        assert_eq!(
            HpAlphabet::LehningerHpc.table()[&b'G'],
            b'h',
            "LehningerHpc: G inherits Lehninger placement"
        );

        // G = polar in ThomasDill, KyteDoolittle, ThomasDillNoC, PBotC1stEd.
        assert_eq!(
            HpAlphabet::ThomasDill.table()[&b'G'],
            b'p',
            "ThomasDill: G polar (Thomas & Dill 1996, Fig. 3)"
        );
        assert_eq!(
            HpAlphabet::KyteDoolittle.table()[&b'G'],
            b'p',
            "KyteDoolittle: G polar (hydropathy -0.4, Kyte & Doolittle 1982)"
        );
        assert_eq!(
            HpAlphabet::ThomasDillNoC.table()[&b'G'],
            b'p',
            "ThomasDillNoC: G polar (inherits ThomasDill)"
        );
        assert_eq!(
            HpAlphabet::PBotC1stEd.table()[&b'G'],
            b'p',
            "PBotC1stEd: G polar (Phillips et al. 2008, Fig. 8.30)"
        );
    }

    #[test]
    fn proline_placement() {
        // P = hydrophobic in Lehninger (nonpolar aliphatic), LehningerCNonpolar, PBotC1stEd.
        assert_eq!(
            HpAlphabet::Lehninger.table()[&b'P'],
            b'h',
            "Lehninger: P nonpolar aliphatic (Nelson & Cox 2021, Ch. 3)"
        );
        assert_eq!(
            HpAlphabet::LehningerCNonpolar.table()[&b'P'],
            b'h',
            "LehningerCNonpolar: P inherits Lehninger placement"
        );
        assert_eq!(
            HpAlphabet::PBotC1stEd.table()[&b'P'],
            b'h',
            "PBotC1stEd: P hydrophobic (Phillips et al. 2008, Fig. 8.30; moved to p in 2nd ed)"
        );
        assert_eq!(
            HpAlphabet::LehningerHpc.table()[&b'P'],
            b'h',
            "LehningerHpc: P inherits Lehninger placement"
        );

        // P = polar in ThomasDill, KyteDoolittle, ThomasDillNoC.
        assert_eq!(
            HpAlphabet::ThomasDill.table()[&b'P'],
            b'p',
            "ThomasDill: P polar (Thomas & Dill 1996, Fig. 3)"
        );
        assert_eq!(
            HpAlphabet::KyteDoolittle.table()[&b'P'],
            b'p',
            "KyteDoolittle: P polar (hydropathy -1.6, Kyte & Doolittle 1982)"
        );
        assert_eq!(
            HpAlphabet::ThomasDillNoC.table()[&b'P'],
            b'p',
            "ThomasDillNoC: P polar (inherits ThomasDill)"
        );
    }

    #[test]
    fn tryptophan_placement() {
        // W = polar in KyteDoolittle only (hydropathy -0.9 due to indole NH).
        assert_eq!(
            HpAlphabet::KyteDoolittle.table()[&b'W'],
            b'p',
            "KyteDoolittle: W polar (hydropathy -0.9, Kyte & Doolittle 1982)"
        );

        // W = hydrophobic in all other named alphabets.
        for alpha in [
            HpAlphabet::Lehninger,
            HpAlphabet::ThomasDill,
            HpAlphabet::ThomasDillNoC,
            HpAlphabet::LehningerCNonpolar,
            HpAlphabet::LehningerHpc,
            HpAlphabet::PBotC1stEd,
        ] {
            assert_eq!(alpha.table()[&b'W'], b'h', "{:?}: W hydrophobic (aromatic)", alpha);
        }
    }

    #[test]
    fn tyrosine_placement() {
        // Y = polar in KyteDoolittle only (hydropathy -1.3 due to hydroxyl).
        assert_eq!(
            HpAlphabet::KyteDoolittle.table()[&b'Y'],
            b'p',
            "KyteDoolittle: Y polar (hydropathy -1.3, Kyte & Doolittle 1982)"
        );

        // Y = hydrophobic in all other named alphabets.
        for alpha in [
            HpAlphabet::Lehninger,
            HpAlphabet::ThomasDill,
            HpAlphabet::ThomasDillNoC,
            HpAlphabet::LehningerCNonpolar,
            HpAlphabet::LehningerHpc,
            HpAlphabet::PBotC1stEd,
        ] {
            assert_eq!(alpha.table()[&b'Y'], b'h', "{:?}: Y hydrophobic (aromatic)", alpha);
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
        for &r in ALL_RESIDUES {
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
        let t = HpAlphabet::RandomControl.table();
        for &r in b"ADGKLMQRWY" {
            assert_eq!(t[&r], b'h', "RandomControl: {} should be h", r as char);
        }
        for &r in b"CEFHINPSTV" {
            assert_eq!(t[&r], b'p', "RandomControl: {} should be p", r as char);
        }
    }

    #[test]
    fn all_alphabet_names_are_unique() {
        let names: Vec<_> = HpAlphabet::all_named().iter().map(|a| a.name()).collect();
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
            (HpAlphabet::Lehninger, "hp_lehninger2"),
            (HpAlphabet::ThomasDill, "hp_thomas_dill2"),
            (HpAlphabet::KyteDoolittle, "hp_kyte_doolittle2"),
            (HpAlphabet::ThomasDillNoC, "hp_thomas_dill_no_c2"),
            (HpAlphabet::LehningerCNonpolar, "hp_lehninger_c_nonpolar2"),
            (HpAlphabet::LehningerHpc, "hp_lehninger_hpc3"),
            (HpAlphabet::PBotC1stEd, "hp_pbotc_1st_ed2"),
            (HpAlphabet::RandomControl, "hp_random_control2"),
        ];

        for (alphabet, moltype) in expected {
            assert_eq!(alphabet.to_moltype(), moltype);
            assert_eq!(HpAlphabet::from_moltype(moltype), Some(alphabet));
        }
    }

    /// The count in the name has to be the number of distinct symbols the table emits,
    /// otherwise the name lies about what the index stores.
    #[test]
    fn class_count_matches_distinct_table_symbols() {
        for alphabet in HpAlphabet::all_named() {
            let symbols: std::collections::HashSet<u8> = alphabet
                .table()
                .iter()
                .filter(|(residue, _)| **residue != b'*')
                .map(|(_, symbol)| *symbol)
                .collect();
            assert_eq!(symbols.len(), alphabet.size(), "{}", alphabet.name());
            assert!(
                alphabet.to_moltype().ends_with(&alphabet.size().to_string()),
                "{}",
                alphabet.name()
            );
        }
    }

    /// Seeded controls put the seed after the class count, so seed 3 of a 2-class control
    /// cannot be misread as a 23-class alphabet.
    #[test]
    fn seeded_random_control_moltype_separates_count_from_seed() {
        assert_eq!(HpAlphabet::Random(3).to_moltype(), "hp_random_control2_3");
        assert_eq!(HpAlphabet::from_moltype("hp_random_control2_3"), Some(HpAlphabet::Random(3)));
        assert_eq!(HpAlphabet::from_moltype("hp_random_control2_10"), Some(HpAlphabet::Random(10)));
    }

    #[test]
    fn from_moltype_rejects_unrelated_strings() {
        for moltype in ["protein", "dayhoff6", "sdm12", "hp_", "hp_nope2"] {
            assert_eq!(HpAlphabet::from_moltype(moltype), None, "{moltype}");
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

        for alphabet in HpAlphabet::all_named() {
            // The randomized control is a negative control, so it has no reason to agree.
            if matches!(alphabet, HpAlphabet::RandomControl) {
                continue;
            }
            let table = alphabet.table();
            for &residue in ALWAYS_HYDROPHOBIC {
                assert_eq!(
                    table[&residue],
                    b'h',
                    "{}: {} should be hydrophobic in every scheme",
                    alphabet.name(),
                    residue as char
                );
            }
            for &residue in ALWAYS_POLAR {
                assert_eq!(
                    table[&residue],
                    b'p',
                    "{}: {} should be polar in every scheme",
                    alphabet.name(),
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
mod multi_letter_tests {
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

/// The residue-to-class table for alphabets kmerseek encodes itself.
///
/// Returns `None` for the three alphabets sourmash encodes: `protein20` (no reduction),
/// `dayhoff6`, and `hp_lehninger2` (sourmash's own `aa_to_hp` partition). Those are hashed
/// by sourmash directly; everything else is pre-encoded through the table returned here.
pub fn alphabet_table(moltype: &str) -> Option<&'static HashMap<u8, u8>> {
    let moltype = canonical_moltype(moltype);
    if let Some(alphabet) = HpAlphabet::from_moltype(moltype) {
        if alphabet.uses_sourmash_encoder() {
            return None;
        }
        return Some(alphabet.table());
    }
    ReducedAlphabet::from_moltype(moltype).map(|alphabet| alphabet.table())
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
            assert!(HpAlphabet::from_moltype(moltype).is_none(), "{moltype}");
            assert!(ReducedAlphabet::from_moltype(moltype).is_none(), "{moltype}");
        }
    }
}
