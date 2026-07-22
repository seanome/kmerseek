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
    LehningerPlusC,
    PBotC1stEd,
    ShuffledControl,
    /// Seeded shuffled control; seed must be 1-10.
    Shuffled(u64),
}

impl HpAlphabet {
    pub fn table(&self) -> &'static HashMap<u8, u8> {
        match self {
            Self::Lehninger => &LEHNINGER_HP,
            Self::ThomasDill => &THOMAS_DILL_HP,
            Self::KyteDoolittle => &KYTE_DOOLITTLE_HP,
            Self::ThomasDillNoC => &THOMAS_DILL_NO_C_HP,
            Self::LehningerPlusC => &LEHNINGER_PLUS_C_HP,
            Self::PBotC1stEd => &PBOTC_1ST_ED_HP,
            Self::ShuffledControl => &SHUFFLED_CONTROL_HP,
            Self::Shuffled(1) => &SHUFFLED_HP_1,
            Self::Shuffled(2) => &SHUFFLED_HP_2,
            Self::Shuffled(3) => &SHUFFLED_HP_3,
            Self::Shuffled(4) => &SHUFFLED_HP_4,
            Self::Shuffled(5) => &SHUFFLED_HP_5,
            Self::Shuffled(6) => &SHUFFLED_HP_6,
            Self::Shuffled(7) => &SHUFFLED_HP_7,
            Self::Shuffled(8) => &SHUFFLED_HP_8,
            Self::Shuffled(9) => &SHUFFLED_HP_9,
            Self::Shuffled(10) => &SHUFFLED_HP_10,
            Self::Shuffled(n) => panic!("shuffled seed {n} not pre-computed (only 1-10 supported)"),
        }
    }

    /// Short identifier used in output filenames and Nextflow traces.
    pub fn name(&self) -> &'static str {
        match self {
            Self::Lehninger => "lehninger",
            Self::ThomasDill => "thomas_dill",
            Self::KyteDoolittle => "kyte_doolittle",
            Self::ThomasDillNoC => "thomas_dill_no_c",
            Self::LehningerPlusC => "lehninger_plus_c",
            Self::PBotC1stEd => "pbotc_1st_ed",
            Self::ShuffledControl => "shuffled_control",
            Self::Shuffled(1) => "shuffled_control_1",
            Self::Shuffled(2) => "shuffled_control_2",
            Self::Shuffled(3) => "shuffled_control_3",
            Self::Shuffled(4) => "shuffled_control_4",
            Self::Shuffled(5) => "shuffled_control_5",
            Self::Shuffled(6) => "shuffled_control_6",
            Self::Shuffled(7) => "shuffled_control_7",
            Self::Shuffled(8) => "shuffled_control_8",
            Self::Shuffled(9) => "shuffled_control_9",
            Self::Shuffled(10) => "shuffled_control_10",
            Self::Shuffled(n) => panic!("shuffled seed {n} not pre-computed (only 1-10 supported)"),
        }
    }

    pub fn all_named() -> &'static [HpAlphabet] {
        &[
            HpAlphabet::Lehninger,
            HpAlphabet::ThomasDill,
            HpAlphabet::KyteDoolittle,
            HpAlphabet::ThomasDillNoC,
            HpAlphabet::LehningerPlusC,
            HpAlphabet::PBotC1stEd,
            HpAlphabet::ShuffledControl,
        ]
    }

    /// Moltype string stored in the index (e.g. `"hp_thomas_dill"`).
    pub fn to_moltype(&self) -> String {
        format!("hp_{}", self.name())
    }

    /// Parse from a moltype string. Returns `None` for `"hp"` (sourmash built-in)
    /// and for unrecognized strings.
    pub fn from_moltype(s: &str) -> Option<HpAlphabet> {
        let name = s.strip_prefix("hp_")?;
        // Try named alphabets first.
        if let Some(a) = Self::all_named().iter().find(|a| a.name() == name) {
            return Some(*a);
        }
        // Try seeded shuffled controls: "shuffled_control_N" -> Shuffled(N).
        if let Some(n_str) = name.strip_prefix("shuffled_control_") {
            if let Ok(n) = n_str.parse::<u64>() {
                return Some(HpAlphabet::Shuffled(n));
            }
        }
        None
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
static LEHNINGER_PLUS_C_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"ACFGILMPVWY", b"DEHKNQRST"));

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
// Shuffled negative control.
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
// statistically robust. See `shuffled_hp()` below.
// -----------------------------------------------------------------------------
static SHUFFLED_CONTROL_HP: LazyLock<HashMap<u8, u8>> =
    LazyLock::new(|| build_hp(b"ADGKLMQRWY", b"CEFHINPSTV"));

// Seeded shuffled controls — 10 independent random partitions for null-distribution estimation.
static SHUFFLED_HP_1: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(1));
static SHUFFLED_HP_2: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(2));
static SHUFFLED_HP_3: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(3));
static SHUFFLED_HP_4: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(4));
static SHUFFLED_HP_5: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(5));
static SHUFFLED_HP_6: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(6));
static SHUFFLED_HP_7: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(7));
static SHUFFLED_HP_8: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(8));
static SHUFFLED_HP_9: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(9));
static SHUFFLED_HP_10: LazyLock<HashMap<u8, u8>> = LazyLock::new(|| shuffled_hp(10));

/// Generate a seeded random HP partition for negative-control runs.
/// Pass multiple seeds to characterize the null distribution.
pub fn shuffled_hp(seed: u64) -> HashMap<u8, u8> {
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
    fn all_mappings_are_h_or_p() {
        for alpha in all_named_alphabets() {
            for &r in ALL_RESIDUES {
                let encoded = alpha.table()[&r];
                assert!(
                    encoded == b'h' || encoded == b'p',
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

        // C = hydrophobic in ThomasDill, KyteDoolittle, LehningerPlusC, PBotC1stEd.
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
            HpAlphabet::LehningerPlusC.table()[&b'C'],
            b'h',
            "LehningerPlusC: C hydrophobic by construction"
        );
        assert_eq!(
            HpAlphabet::PBotC1stEd.table()[&b'C'],
            b'h',
            "PBotC1stEd: C borderline-h (Phillips et al. 2008, Fig. 8.30)"
        );
    }

    #[test]
    fn glycine_placement() {
        // G = hydrophobic in Lehninger (nonpolar aliphatic) and LehningerPlusC.
        assert_eq!(
            HpAlphabet::Lehninger.table()[&b'G'],
            b'h',
            "Lehninger: G nonpolar aliphatic (Nelson & Cox 2021, Ch. 3)"
        );
        assert_eq!(
            HpAlphabet::LehningerPlusC.table()[&b'G'],
            b'h',
            "LehningerPlusC: G inherits Lehninger placement"
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
        // P = hydrophobic in Lehninger (nonpolar aliphatic), LehningerPlusC, PBotC1stEd.
        assert_eq!(
            HpAlphabet::Lehninger.table()[&b'P'],
            b'h',
            "Lehninger: P nonpolar aliphatic (Nelson & Cox 2021, Ch. 3)"
        );
        assert_eq!(
            HpAlphabet::LehningerPlusC.table()[&b'P'],
            b'h',
            "LehningerPlusC: P inherits Lehninger placement"
        );
        assert_eq!(
            HpAlphabet::PBotC1stEd.table()[&b'P'],
            b'h',
            "PBotC1stEd: P hydrophobic (Phillips et al. 2008, Fig. 8.30; moved to p in 2nd ed)"
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
            HpAlphabet::LehningerPlusC,
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
            HpAlphabet::LehningerPlusC,
            HpAlphabet::PBotC1stEd,
        ] {
            assert_eq!(alpha.table()[&b'Y'], b'h', "{:?}: Y hydrophobic (aromatic)", alpha);
        }
    }

    #[test]
    fn shuffled_hp_deterministic() {
        let a = shuffled_hp(42);
        let b = shuffled_hp(42);
        assert_eq!(a, b, "shuffled_hp must be deterministic for the same seed");
    }

    #[test]
    fn shuffled_hp_covers_all_residues() {
        let m = shuffled_hp(0);
        assert_eq!(m.len(), 21);
        for &r in ALL_RESIDUES {
            assert!(m.contains_key(&r), "shuffled_hp missing residue {}", r as char);
            let v = m[&r];
            assert!(v == b'h' || v == b'p');
        }
    }

    #[test]
    fn shuffled_hp_differs_by_seed() {
        let a = shuffled_hp(1);
        let b = shuffled_hp(2);
        // Two different seeds should produce different partitions
        // (astronomically unlikely to collide with 20 elements)
        assert_ne!(a, b, "shuffled_hp with different seeds should differ");
    }

    #[test]
    fn shuffled_control_static_matches_documented_partition() {
        // h: ADGKLMQRWY  p: CEFHINPSTV — frozen in static for reproducibility.
        let t = HpAlphabet::ShuffledControl.table();
        for &r in b"ADGKLMQRWY" {
            assert_eq!(t[&r], b'h', "ShuffledControl: {} should be h", r as char);
        }
        for &r in b"CEFHINPSTV" {
            assert_eq!(t[&r], b'p', "ShuffledControl: {} should be p", r as char);
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
}
