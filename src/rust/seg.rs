//! Low-complexity masking in the manner of SEG (Wootton & Federhen 1993, Computers &
//! Chemistry 17:149-163), for `--seg-mask` in `kmerseek search`.
//!
//! SEG slides a window of 12 residues along a protein and computes each window's Shannon
//! entropy over the 20 amino acids, in bits. A window at or below 2.2 bits starts a
//! low-complexity segment; the segment grows over every neighbouring window at or below
//! 2.5 bits; and every residue those windows cover is masked. These are the defaults of NCBI
//! `segmasker` and of BLASTP's `-seg yes` (window 12, locut 2.2, hicut 2.5).
//!
//! SEG then trims each segment to its least probable subsequence. That step is not done
//! here, so a masked stretch can run a few residues past where `segmasker` would stop it.

/// Residues per window.
pub const SEG_WINDOW: usize = 12;
/// A window at or below this entropy, in bits, starts a segment.
pub const SEG_LOCUT: f64 = 2.2;
/// A window at or below this entropy, in bits, extends a segment it touches.
pub const SEG_HICUT: f64 = 2.5;

/// One flag per residue of `residues`: true where it lies in a low-complexity segment.
/// Letters are compared case-insensitively; any byte counts as its own residue type. A
/// protein shorter than one window is never masked.
pub fn seg_mask(residues: &[u8]) -> Vec<bool> {
    let n = residues.len();
    let mut mask = vec![false; n];
    if n < SEG_WINDOW {
        return mask;
    }
    let entropies = window_entropies(residues);
    let mut w = 0;
    while w < entropies.len() {
        if entropies[w] > SEG_HICUT {
            w += 1;
            continue;
        }
        // A run of windows at or below hicut; masked only if one of them is at or below locut.
        let first = w;
        while w < entropies.len() && entropies[w] <= SEG_HICUT {
            w += 1;
        }
        if entropies[first..w].iter().any(|&h| h <= SEG_LOCUT) {
            mask[first..w - 1 + SEG_WINDOW].fill(true);
        }
    }
    mask
}

/// Shannon entropy in bits of every `SEG_WINDOW`-residue window, in order of start.
///
/// H = log2(W) - (1/W) sum_c n_c log2(n_c), with n_c the window's count of residue c. Sliding
/// the window by one changes two counts, so the sum is updated from those two terms rather
/// than recounted over every residue type.
fn window_entropies(residues: &[u8]) -> Vec<f64> {
    let mut counts = [0u32; 256];
    let key = |b: u8| b.to_ascii_uppercase() as usize;
    let term = |n: u32| if n == 0 { 0.0 } else { n as f64 * (n as f64).log2() };
    let w = SEG_WINDOW as f64;
    let mut sum = 0.0;
    let mut add = |counts: &mut [u32; 256], c: usize, delta: i32| {
        sum -= term(counts[c]);
        counts[c] = counts[c].wrapping_add_signed(delta);
        sum += term(counts[c]);
        sum
    };
    let mut total = 0.0;
    for &b in &residues[..SEG_WINDOW] {
        total = add(&mut counts, key(b), 1);
    }
    let mut out = Vec::with_capacity(residues.len() + 1 - SEG_WINDOW);
    out.push((w.log2() - total / w).max(0.0));
    for i in SEG_WINDOW..residues.len() {
        add(&mut counts, key(residues[i - SEG_WINDOW]), -1);
        let total = add(&mut counts, key(residues[i]), 1);
        out.push((w.log2() - total / w).max(0.0));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Twelve different residues have log2(12) = 3.58 bits, above both cutoffs: nothing is
    /// masked. A run of one residue has 0 bits and is masked end to end.
    #[test]
    fn masks_a_homopolymer_and_not_varied_sequence() {
        let varied = b"ACDEFGHIKLMNPQRSTVWY";
        assert!(seg_mask(varied).iter().all(|&m| !m));
        let mut seq = varied.to_vec();
        seq.extend_from_slice(&[b'Q'; 20]);
        seq.extend_from_slice(varied);
        let mask = seg_mask(&seq);
        // Every Q is masked; the varied ends are masked only where a window that reaches
        // into the Q run is still at or below 2.5 bits.
        assert!(mask[20..40].iter().all(|&m| m));
        assert!(!mask[0] && !mask[59]);
    }

    /// The window entropies, by hand: "SDSSDSSSDSSD" holds 8 S and 4 D, so
    /// -(2/3 log2 2/3 + 1/3 log2 1/3) = 0.918 bits.
    #[test]
    fn window_entropy_by_hand() {
        let h = window_entropies(b"SDSSDSSSDSSD");
        assert_eq!(h.len(), 1);
        assert!((h[0] - 0.918_295_834_054_489_6).abs() < 1e-12);
    }

    /// The sliding update gives the entropy a fresh count gives, on every window of a real
    /// low-complexity stretch (DSPP's decoy, notebook 344) and of varied sequence.
    #[test]
    fn sliding_entropy_matches_a_recount() {
        let seq =
            b"ITDSDSQRESDSSDSNSSDSDSSESNSSDDSSDSDSDSDSSSDSDSSSDSSDSSSDSTIQACDEFGHIKLMNPQRSTVWY";
        let fresh = |window: &[u8]| -> f64 {
            let mut counts = [0u32; 256];
            window.iter().for_each(|&b| counts[b as usize] += 1);
            let w = window.len() as f64;
            counts
                .iter()
                .filter(|&&c| c > 0)
                .map(|&c| -(c as f64 / w) * (c as f64 / w).log2())
                .sum()
        };
        let sliding = window_entropies(seq);
        assert_eq!(sliding.len(), seq.len() - SEG_WINDOW + 1);
        for (i, h) in sliding.iter().enumerate() {
            assert!((h - fresh(&seq[i..i + SEG_WINDOW])).abs() < 1e-9, "window {i}");
        }
    }

    /// Lower case is the same residue, and a protein shorter than one window is left alone.
    #[test]
    fn short_and_lower_case() {
        assert_eq!(seg_mask(b"QQQQQ"), vec![false; 5]);
        assert_eq!(seg_mask(b"qqqqqqqqqqqqQQQQ"), seg_mask(b"QQQQQQQQQQQQQQQQ"));
    }
}
