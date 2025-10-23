use std::collections::HashSet;

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
pub fn abundance_stats(
    intersection: &HashSet<u64>,
    query_abunds: &[u64],
    target_abunds: &[u64],
) -> (f64, f64, f64) {
    let mut abunds = Vec::new();

    // Get abundances for intersecting k-mers
    for (i, &_min) in intersection.iter().enumerate() {
        if let (Some(&query_abund), Some(&target_abund)) =
            (query_abunds.get(i), target_abunds.get(i))
        {
            abunds.push((query_abund + target_abund) as f64 / 2.0);
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
