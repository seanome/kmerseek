use std::collections::{HashMap, HashSet};

/// Calculate ANI (Average Nucleotide Identity) from containment
pub fn ani(containment: f64, _size: usize) -> f64 {
    // Simplified ANI calculation based on containment
    // This is a rough approximation - in practice, ANI calculation is more complex
    if containment <= 0.0 {
        0.0
    } else {
        // Use a logarithmic relationship for ANI
        let ani = 1.0 - (-containment.ln()).exp();
        ani.clamp(0.0, 1.0)
    }
}

/// Calculate abundance statistics for intersecting k-mers
///
/// WHY: This function takes the mins arrays because abundances are stored in the same order
/// as the mins. We need to find the position of each intersecting hash in both the query
/// and target mins arrays to get the correct corresponding abundances. Using enumerate()
/// on a HashSet would give arbitrary indices that don't correspond to the actual positions
/// in the minhash arrays.
///
/// # Arguments
/// * `intersection` - Set of hash values that appear in both query and target
/// * `query_mins` - Array of hash values from the query minhash (aligned with query_abunds)
/// * `query_abunds` - Array of abundances from the query minhash
/// * `target_mins` - Array of hash values from the target minhash (aligned with target_abunds)
/// * `target_abunds` - Array of abundances from the target minhash
pub fn abundance_stats(
    intersection: &HashSet<u64>,
    query_mins: &[u64],
    query_abunds: &[u64],
    target_mins: &[u64],
    target_abunds: &[u64],
) -> (f64, f64, f64) {
    let mut abunds = Vec::new();

    // Build index maps for O(1) lookup: hash value -> position in mins array
    // WHY: This is more efficient than calling position() for each hash in the intersection.
    // Building the maps is O(n) and lookups are O(1), making the overall complexity O(n)
    // instead of O(n*m) where n is intersection size and m is mins array size.
    let query_index: HashMap<u64, usize> =
        query_mins.iter().enumerate().map(|(i, &hash)| (hash, i)).collect();
    let target_index: HashMap<u64, usize> =
        target_mins.iter().enumerate().map(|(i, &hash)| (hash, i)).collect();

    // Get abundances for intersecting k-mers
    // WHY: We iterate over the intersection HashSet and for each hash, we look up its position
    // in both the query and target index maps. This gives us the correct indices to use
    // for accessing the abundance arrays. This is necessary because HashSet iteration order
    // is undefined, so we can't use enumerate() indices.
    for &hashval in intersection {
        // Look up the position of this hash in both index maps
        if let (Some(&q_idx), Some(&t_idx)) =
            (query_index.get(&hashval), target_index.get(&hashval))
        {
            // Get the corresponding abundances using the correct indices
            if let (Some(&query_abund), Some(&target_abund)) =
                (query_abunds.get(q_idx), target_abunds.get(t_idx))
            {
                // Average the abundances from query and target
                abunds.push((query_abund + target_abund) as f64 / 2.0);
            }
        }
    }

    if abunds.is_empty() {
        return (1.0, 1.0, 0.0);
    }

    // Calculate statistics
    let sum: f64 = abunds.iter().sum();
    let average = sum / abunds.len() as f64;

    abunds.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = if abunds.len() % 2 == 0 {
        (abunds[abunds.len() / 2 - 1] + abunds[abunds.len() / 2]) / 2.0
    } else {
        abunds[abunds.len() / 2]
    };

    let variance: f64 =
        abunds.iter().map(|&x| (x - average).powi(2)).sum::<f64>() / abunds.len() as f64;
    let std_dev = variance.sqrt();

    (average, median, std_dev)
}

/// Calculate weighted metrics, weighted by abundances:
/// - Weighted fraction of target in query
pub fn weighted_fraction_target_in_query(
    query_abunds: Option<&[u64]>,
    target_abunds: Option<&[u64]>,
) -> f64 {
    let f_weighted_target_in_query =
        if let (Some(query_abunds), Some(target_abunds)) = (query_abunds, target_abunds) {
            let query_weight: f64 = query_abunds.iter().sum::<u64>() as f64;
            let target_weight: f64 = target_abunds.iter().sum::<u64>() as f64;
            if query_weight > 0.0 {
                target_weight / query_weight
            } else {
                0.0
            }
        } else {
            1.0
        };

    f_weighted_target_in_query
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn test_ani_non_positive_containment() {
        // Both the zero and negative branches return 0.0.
        assert_eq!(ani(0.0, 100), 0.0);
        assert_eq!(ani(-0.25, 100), 0.0);
    }

    #[test]
    fn test_ani_clamps_and_computes_exact() {
        // For containment in (0, 1] the raw value 1 - exp(-ln c) is <= 0, so it clamps to 0.0.
        assert_eq!(ani(0.5, 100), 0.0); // raw = 1 - exp(ln 2) = -1.0, clamped
        assert_eq!(ani(1.0, 100), 0.0); // raw = 1 - exp(0) = 0.0
                                        // For containment > 1 the value lands in (0, 1): 1 - exp(-ln 2) = 0.5.
        assert!((ani(2.0, 100) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_abundance_stats_empty_intersection_returns_defaults() {
        let empty: HashSet<u64> = HashSet::new();
        let out = abundance_stats(&empty, &[1, 2], &[3, 4], &[1, 2], &[5, 6]);
        assert_eq!(out, (1.0, 1.0, 0.0));
    }

    #[test]
    fn test_abundance_stats_odd_count() {
        // Per-hash averaged abundances: 10->(2+4)/2=3, 20->(4+4)/2=4, 30->(6+4)/2=5.
        let inter: HashSet<u64> = [10, 20, 30].into_iter().collect();
        let (avg, med, sd) =
            abundance_stats(&inter, &[10, 20, 30], &[2, 4, 6], &[30, 20, 10], &[4, 4, 4]);
        assert!((avg - 4.0).abs() < EPS);
        assert!((med - 4.0).abs() < EPS);
        // variance of {3,4,5} = 2/3
        assert!((sd - (2.0f64 / 3.0).sqrt()).abs() < EPS);
    }

    #[test]
    fn test_abundance_stats_even_count_median_is_midpoint() {
        // 10->(2+4)/2=3, 20->(4+8)/2=6; median = (3+6)/2 = 4.5.
        let inter: HashSet<u64> = [10, 20].into_iter().collect();
        let (avg, med, _sd) = abundance_stats(&inter, &[10, 20], &[2, 4], &[10, 20], &[4, 8]);
        assert!((avg - 4.5).abs() < EPS);
        assert!((med - 4.5).abs() < EPS);
    }

    #[test]
    fn test_weighted_fraction_target_in_query() {
        // Both present, positive query weight: target/query = 4/6.
        let f = weighted_fraction_target_in_query(Some(&[1, 2, 3]), Some(&[2, 2]));
        assert!((f - 4.0 / 6.0).abs() < EPS);
        // Zero query weight -> 0.0.
        assert_eq!(weighted_fraction_target_in_query(Some(&[0, 0]), Some(&[2, 2])), 0.0);
        // Missing abundances -> 1.0 (either side).
        assert_eq!(weighted_fraction_target_in_query(None, Some(&[2, 2])), 1.0);
        assert_eq!(weighted_fraction_target_in_query(Some(&[1, 2]), None), 1.0);
    }
}
