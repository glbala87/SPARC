//! Accuracy metrics for truthset validation
//!
//! Classification metrics (barcode detection), correlation metrics (expression),
//! and clustering metrics (ARI, NMI).

use ahash::AHashMap;

// ─── Classification metrics ──────────────────────────────────────────

/// Sensitivity (true positive rate) = TP / (TP + FN)
pub fn sensitivity(tp: u64, fn_count: u64) -> f64 {
    let denom = tp + fn_count;
    if denom == 0 {
        return 0.0;
    }
    tp as f64 / denom as f64
}

/// Specificity (true negative rate) = TN / (TN + FP)
pub fn specificity(tn: u64, fp: u64) -> f64 {
    let denom = tn + fp;
    if denom == 0 {
        return 0.0;
    }
    tn as f64 / denom as f64
}

/// Precision = TP / (TP + FP)
pub fn precision(tp: u64, fp: u64) -> f64 {
    let denom = tp + fp;
    if denom == 0 {
        return 0.0;
    }
    tp as f64 / denom as f64
}

/// Recall = TP / (TP + FN)  (same as sensitivity)
pub fn recall(tp: u64, fn_count: u64) -> f64 {
    sensitivity(tp, fn_count)
}

/// F1 score = 2 * precision * recall / (precision + recall)
pub fn f1_score(prec: f64, rec: f64) -> f64 {
    let denom = prec + rec;
    if denom == 0.0 {
        return 0.0;
    }
    2.0 * prec * rec / denom
}

// ─── Expression quantification metrics ───────────────────────────────

/// Pearson correlation coefficient
pub fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "Vectors must have equal length");
    let n = x.len() as f64;
    if n < 2.0 {
        return 0.0;
    }

    let mean_x = x.iter().sum::<f64>() / n;
    let mean_y = y.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_x = 0.0;
    let mut var_y = 0.0;

    for (&xi, &yi) in x.iter().zip(y.iter()) {
        let dx = xi - mean_x;
        let dy = yi - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    let denom = (var_x * var_y).sqrt();
    if denom < 1e-15 {
        return 0.0;
    }
    cov / denom
}

/// Spearman rank correlation coefficient
pub fn spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "Vectors must have equal length");
    let ranks_x = rank(x);
    let ranks_y = rank(y);
    pearson_correlation(&ranks_x, &ranks_y)
}

/// Compute ranks (average rank for ties)
fn rank(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && (indexed[j].1 - indexed[i].1).abs() < 1e-15 {
            j += 1;
        }
        // Average rank for ties (1-indexed)
        let avg_rank = (i + j) as f64 / 2.0 + 0.5;
        for item in indexed.iter().take(j).skip(i) {
            ranks[item.0] = avg_rank;
        }
        i = j;
    }
    ranks
}

/// Mean absolute error
pub fn mean_absolute_error(observed: &[f64], expected: &[f64]) -> f64 {
    assert_eq!(observed.len(), expected.len());
    if observed.is_empty() {
        return 0.0;
    }
    let sum: f64 = observed
        .iter()
        .zip(expected.iter())
        .map(|(o, e)| (o - e).abs())
        .sum();
    sum / observed.len() as f64
}

/// Root mean squared error
pub fn rmse(observed: &[f64], expected: &[f64]) -> f64 {
    assert_eq!(observed.len(), expected.len());
    if observed.is_empty() {
        return 0.0;
    }
    let mse: f64 = observed
        .iter()
        .zip(expected.iter())
        .map(|(o, e)| (o - e).powi(2))
        .sum::<f64>()
        / observed.len() as f64;
    mse.sqrt()
}

// ─── Clustering metrics ──────────────────────────────────────────────

/// Adjusted Rand Index (ARI)
///
/// Measures similarity between two clusterings, adjusted for chance.
/// Range: [-1, 1], where 1 = perfect agreement, 0 = random.
pub fn adjusted_rand_index(labels_true: &[usize], labels_pred: &[usize]) -> f64 {
    assert_eq!(labels_true.len(), labels_pred.len());
    let n = labels_true.len();
    if n == 0 {
        return 0.0;
    }

    // Build contingency table
    let mut contingency: AHashMap<(usize, usize), u64> = AHashMap::new();
    for (&t, &p) in labels_true.iter().zip(labels_pred.iter()) {
        *contingency.entry((t, p)).or_insert(0) += 1;
    }

    // Row and column sums
    let mut row_sums: AHashMap<usize, u64> = AHashMap::new();
    let mut col_sums: AHashMap<usize, u64> = AHashMap::new();
    for (&(t, p), &count) in &contingency {
        *row_sums.entry(t).or_insert(0) += count;
        *col_sums.entry(p).or_insert(0) += count;
    }

    // Sum of C(n_ij, 2) over all cells
    let sum_comb: f64 = contingency.values().map(|&v| comb2(v)).sum();

    // Sum of C(a_i, 2) for row sums
    let sum_comb_rows: f64 = row_sums.values().map(|&v| comb2(v)).sum();

    // Sum of C(b_j, 2) for column sums
    let sum_comb_cols: f64 = col_sums.values().map(|&v| comb2(v)).sum();

    let n_comb = comb2(n as u64);
    if n_comb == 0.0 {
        return 0.0;
    }

    let expected = sum_comb_rows * sum_comb_cols / n_comb;
    let max_index = (sum_comb_rows + sum_comb_cols) / 2.0;
    let denom = max_index - expected;

    if denom.abs() < 1e-15 {
        return if (sum_comb - expected).abs() < 1e-15 {
            1.0
        } else {
            0.0
        };
    }

    (sum_comb - expected) / denom
}

/// Normalized Mutual Information (NMI)
///
/// Measures mutual information between two clusterings, normalized to [0, 1].
pub fn normalized_mutual_info(labels_true: &[usize], labels_pred: &[usize]) -> f64 {
    assert_eq!(labels_true.len(), labels_pred.len());
    let n = labels_true.len() as f64;
    if n < 1.0 {
        return 0.0;
    }

    // Build contingency table
    let mut contingency: AHashMap<(usize, usize), u64> = AHashMap::new();
    for (&t, &p) in labels_true.iter().zip(labels_pred.iter()) {
        *contingency.entry((t, p)).or_insert(0) += 1;
    }

    let mut row_sums: AHashMap<usize, u64> = AHashMap::new();
    let mut col_sums: AHashMap<usize, u64> = AHashMap::new();
    for (&(t, p), &count) in &contingency {
        *row_sums.entry(t).or_insert(0) += count;
        *col_sums.entry(p).or_insert(0) += count;
    }

    // Entropy of true labels
    let h_true: f64 = row_sums
        .values()
        .map(|&a| {
            let p = a as f64 / n;
            if p > 0.0 { -p * p.ln() } else { 0.0 }
        })
        .sum();

    // Entropy of predicted labels
    let h_pred: f64 = col_sums
        .values()
        .map(|&b| {
            let p = b as f64 / n;
            if p > 0.0 { -p * p.ln() } else { 0.0 }
        })
        .sum();

    // Mutual information
    let mut mi = 0.0;
    for (&(t, p), &nij) in &contingency {
        if nij == 0 {
            continue;
        }
        let pij = nij as f64 / n;
        let pi = row_sums[&t] as f64 / n;
        let pj = col_sums[&p] as f64 / n;
        mi += pij * (pij / (pi * pj)).ln();
    }

    let denom = ((h_true + h_pred) / 2.0).max(1e-15);
    (mi / denom).clamp(0.0, 1.0)
}

/// Binomial coefficient C(n, 2) = n*(n-1)/2
fn comb2(n: u64) -> f64 {
    if n < 2 {
        return 0.0;
    }
    (n as f64) * (n as f64 - 1.0) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classification_metrics() {
        // TP=80, FP=10, TN=90, FN=20
        let tp = 80;
        let fp = 10;
        let tn = 90;
        let fn_count = 20;

        assert!((sensitivity(tp, fn_count) - 0.8).abs() < 1e-10);
        assert!((specificity(tn, fp) - 0.9).abs() < 1e-10);
        assert!((precision(tp, fp) - 80.0 / 90.0).abs() < 1e-10);
        assert!((recall(tp, fn_count) - 0.8).abs() < 1e-10);

        let p = precision(tp, fp);
        let r = recall(tp, fn_count);
        let f1 = f1_score(p, r);
        assert!(f1 > 0.0 && f1 < 1.0);
    }

    #[test]
    fn test_pearson_perfect() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        assert!((pearson_correlation(&x, &y) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pearson_negative() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![10.0, 8.0, 6.0, 4.0, 2.0];
        assert!((pearson_correlation(&x, &y) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_spearman_monotonic() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![10.0, 20.0, 30.0, 40.0, 50.0];
        assert!((spearman_correlation(&x, &y) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mae_and_rmse() {
        let obs = vec![1.0, 2.0, 3.0];
        let exp = vec![1.0, 2.0, 3.0];
        assert!(mean_absolute_error(&obs, &exp).abs() < 1e-10);
        assert!(rmse(&obs, &exp).abs() < 1e-10);

        let obs2 = vec![1.0, 2.0, 3.0];
        let exp2 = vec![2.0, 3.0, 4.0];
        assert!((mean_absolute_error(&obs2, &exp2) - 1.0).abs() < 1e-10);
        assert!((rmse(&obs2, &exp2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ari_perfect() {
        let t = vec![0, 0, 0, 1, 1, 1];
        let p = vec![0, 0, 0, 1, 1, 1];
        assert!((adjusted_rand_index(&t, &p) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_ari_permuted() {
        // ARI should be 1.0 even with relabeled clusters
        let t = vec![0, 0, 0, 1, 1, 1];
        let p = vec![1, 1, 1, 0, 0, 0];
        assert!((adjusted_rand_index(&t, &p) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_nmi_perfect() {
        let t = vec![0, 0, 0, 1, 1, 1];
        let p = vec![0, 0, 0, 1, 1, 1];
        assert!((normalized_mutual_info(&t, &p) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_nmi_permuted() {
        let t = vec![0, 0, 0, 1, 1, 1];
        let p = vec![1, 1, 1, 0, 0, 0];
        assert!((normalized_mutual_info(&t, &p) - 1.0).abs() < 1e-10);
    }
}
