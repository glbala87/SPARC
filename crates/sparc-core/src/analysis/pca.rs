//! Principal Component Analysis via power iteration

/// Simple PCA implementation using power iteration
/// Returns (components, explained_variance) where components is n_cells x n_pcs
pub fn pca(data: &[Vec<f64>], n_components: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    let n_cells = data.len();
    let n_genes = if n_cells > 0 { data[0].len() } else { 0 };
    let n_pcs = n_components.min(n_cells).min(n_genes);

    let mut components = Vec::with_capacity(n_pcs);
    let mut variances = Vec::with_capacity(n_pcs);

    // Work with a mutable copy for deflation
    let mut residual: Vec<Vec<f64>> = data.to_vec();

    for _ in 0..n_pcs {
        // Initialize vector
        let mut v = vec![1.0 / (n_genes as f64).sqrt(); n_genes];

        // Power iteration
        for _ in 0..30 {
            // u = X * v
            let mut u: Vec<f64> = residual
                .iter()
                .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
                .collect();

            let u_norm: f64 = u.iter().map(|x| x * x).sum::<f64>().sqrt();
            if u_norm > 1e-12 {
                for val in u.iter_mut() {
                    *val /= u_norm;
                }
            }

            // v = X^T * u
            v = vec![0.0; n_genes];
            for (i, row) in residual.iter().enumerate() {
                for (j, val) in row.iter().enumerate() {
                    v[j] += val * u[i];
                }
            }

            let v_norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            if v_norm > 1e-12 {
                for val in v.iter_mut() {
                    *val /= v_norm;
                }
            }
        }

        // Compute projection scores
        let scores: Vec<f64> = residual
            .iter()
            .map(|row| row.iter().zip(v.iter()).map(|(a, b)| a * b).sum())
            .collect();

        let variance: f64 = scores.iter().map(|x| x * x).sum::<f64>() / n_cells as f64;
        variances.push(variance);
        components.push(scores);

        // Deflate
        for (i, row) in residual.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                *val -= components.last().unwrap()[i] * v[j];
            }
        }
    }

    // Transpose to n_cells x n_pcs
    let mut result = vec![vec![0.0; n_pcs]; n_cells];
    for (pc, scores) in components.iter().enumerate() {
        for (cell, score) in scores.iter().enumerate() {
            result[cell][pc] = *score;
        }
    }

    (result, variances)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pca_basic() {
        let data = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![1.0, 1.0, 0.0],
        ];

        let (components, variances) = pca(&data, 2);

        assert_eq!(components.len(), 4);
        assert_eq!(components[0].len(), 2);
        assert_eq!(variances.len(), 2);
        // First PC should explain more variance
        assert!(variances[0] >= variances[1]);
    }

    #[test]
    fn test_pca_empty() {
        let data: Vec<Vec<f64>> = Vec::new();
        let (components, variances) = pca(&data, 2);
        assert!(components.is_empty());
        assert!(variances.is_empty());
    }
}
