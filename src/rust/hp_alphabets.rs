//! HP alphabet mappings for the Kmerseek robustness sweep.
//!
//! Each alphabet partitions the 20 canonical amino acids into
//! hydrophobic (h) and polar (p) classes. The alphabets differ on
//! borderline residues — principally C, G, P, W, and Y.

use std::collections::HashMap;
use std::sync::LazyLock;

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
    /// [`crate::reduced_alphabets`] (`sdm12`, `gbmr4`). Every alphabet states its size.
    pub fn to_moltype(&self) -> String {
        match self {
            // The seed goes after the class count, so "…control2_3" is seed 3 of a
            // 2-class control and never reads as class count 23.
            Self::Random(n) => format!("hp_random_control2_{n}"),
            _ => format!("hp_{}{}", self.name(), self.size()),
        }
    }

    /// Parse from a moltype string. Returns `None` for unrecognized strings.
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
mod tests {
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
