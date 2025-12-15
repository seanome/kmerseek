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
