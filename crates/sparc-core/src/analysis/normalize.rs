//! Normalization and log transformation

/// Normalize each cell to target_sum total counts
/// data is cells x genes
pub fn normalize_total(data: &mut [Vec<f64>], target_sum: f64) {
    for cell in data.iter_mut() {
        let total: f64 = cell.iter().sum();
        if total > 0.0 {
            let scale = target_sum / total;
            for val in cell.iter_mut() {
                *val *= scale;
            }
        }
    }
}

/// Apply log1p transformation
pub fn log1p_transform(data: &mut [Vec<f64>]) {
    for cell in data.iter_mut() {
        for val in cell.iter_mut() {
            *val = (*val + 1.0).ln();
        }
    }
}

/// Scale data to zero mean and unit variance per gene
pub fn scale(data: &mut [Vec<f64>], max_value: Option<f64>) {
    if data.is_empty() {
        return;
    }
    let n_cells = data.len();
    let n_genes = data[0].len();

    for g in 0..n_genes {
        let mean: f64 = data.iter().map(|c| c[g]).sum::<f64>() / n_cells as f64;
        let var: f64 =
            data.iter().map(|c| (c[g] - mean).powi(2)).sum::<f64>() / n_cells as f64;
        let std = var.sqrt().max(1e-12);

        for cell in data.iter_mut() {
            cell[g] = (cell[g] - mean) / std;
            if let Some(max) = max_value {
                cell[g] = cell[g].min(max).max(-max);
            }
        }
    }
}

/// Find highly variable genes by mean-variance relationship
/// Returns indices of highly variable genes
pub fn highly_variable_genes(
    data: &[Vec<f64>],
    min_mean: f64,
    max_mean: f64,
    min_disp: f64,
) -> Vec<usize> {
    if data.is_empty() {
        return Vec::new();
    }
    let n_cells = data.len();
    let n_genes = data[0].len();
    let mut hvg_indices = Vec::new();

    for g in 0..n_genes {
        let mean: f64 = data.iter().map(|c| c[g]).sum::<f64>() / n_cells as f64;
        let var: f64 =
            data.iter().map(|c| (c[g] - mean).powi(2)).sum::<f64>() / n_cells as f64;

        let disp = if mean > 0.0 { var / mean } else { 0.0 };

        if mean >= min_mean && mean <= max_mean && disp >= min_disp {
            hvg_indices.push(g);
        }
    }

    hvg_indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize_total() {
        let mut data = vec![vec![10.0, 20.0, 30.0], vec![5.0, 5.0, 10.0]];

        normalize_total(&mut data, 100.0);

        let sum0: f64 = data[0].iter().sum();
        let sum1: f64 = data[1].iter().sum();
        assert!((sum0 - 100.0).abs() < 0.001);
        assert!((sum1 - 100.0).abs() < 0.001);
    }

    #[test]
    fn test_log1p() {
        let mut data = vec![vec![0.0, 1.0, 9.0]];
        log1p_transform(&mut data);
        assert!((data[0][0] - 0.0_f64.ln_1p()).abs() < 0.001);
        assert!((data[0][1] - 1.0_f64.ln_1p()).abs() < 0.001);
    }

    #[test]
    fn test_scale() {
        let mut data = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
        scale(&mut data, None);

        // Mean should be ~0 for each gene
        let mean_g0: f64 = data.iter().map(|c| c[0]).sum::<f64>() / 3.0;
        assert!(mean_g0.abs() < 0.001);
    }
}
