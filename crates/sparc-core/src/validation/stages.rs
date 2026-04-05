//! Per-stage validators for truthset validation
//!
//! Each validator runs a pipeline stage against synthetic data and compares
//! results to ground truth, producing accuracy metrics.

use serde::{Deserialize, Serialize};

use crate::barcode::{BarcodeCorrector, BarcodeMatch};
use crate::count::CountMatrix;
use crate::validation::metrics;
use crate::validation::synthetic::SyntheticDataset;

// ─── Extract validation ──────────────────────────────────────────────

/// Results from barcode extraction validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtractValidationResult {
    /// True positives (valid barcodes correctly detected)
    pub tp: u64,
    /// False positives (invalid barcodes incorrectly accepted)
    pub fp: u64,
    /// True negatives (invalid barcodes correctly rejected)
    pub tn: u64,
    /// False negatives (valid barcodes missed)
    pub fn_count: u64,
    /// Sensitivity = TP / (TP + FN)
    pub sensitivity: f64,
    /// Specificity = TN / (TN + FP)
    pub specificity: f64,
    /// Precision = TP / (TP + FP)
    pub precision: f64,
    /// Recall = TP / (TP + FN)
    pub recall: f64,
    /// F1 score
    pub f1: f64,
    /// Barcode correction accuracy (correctly corrected / total mutated)
    pub correction_accuracy: f64,
    /// Total reads processed
    pub total_reads: usize,
}

/// Validate barcode extraction against ground truth
pub fn validate_extract(dataset: &SyntheticDataset) -> ExtractValidationResult {
    let corrector = BarcodeCorrector::new(dataset.whitelist.clone(), 1);

    let mut tp: u64 = 0;
    let mut fp: u64 = 0;
    let mut tn: u64 = 0;
    let mut fn_count: u64 = 0;
    let mut correct_corrections: u64 = 0;
    let mut total_mutations: u64 = 0;

    // Build lookup for mutated -> original barcode
    let mut mutation_map = ahash::AHashMap::new();
    for (original, mutated, _) in &dataset.truth.mutated_barcodes {
        mutation_map.insert(mutated.clone(), original.clone());
    }

    let invalid_set: ahash::AHashSet<&str> = dataset
        .truth
        .invalid_barcodes
        .iter()
        .map(|s| s.as_str())
        .collect();

    for r1 in &dataset.r1_records {
        // Extract barcode from R1
        let barcode_bytes = match r1.subsequence(
            dataset.read_structure.barcode_start,
            dataset.read_structure.barcode_len,
        ) {
            Some(b) => b,
            None => continue,
        };
        let barcode_str = String::from_utf8_lossy(barcode_bytes).to_string();

        let result = corrector.match_barcode(&barcode_str);

        let is_in_whitelist = dataset.whitelist.contains(&barcode_str);
        let is_mutated = mutation_map.contains_key(&barcode_str);
        let is_invalid = invalid_set.contains(barcode_str.as_str());

        match &result {
            BarcodeMatch::Exact(_) if is_in_whitelist => tp += 1,
            BarcodeMatch::Exact(_) => fp += 1,
            BarcodeMatch::Corrected(_, corrected, _) => {
                if is_mutated {
                    let expected_original = &mutation_map[&barcode_str];
                    if corrected == expected_original {
                        tp += 1;
                        correct_corrections += 1;
                    } else {
                        fp += 1;
                    }
                    total_mutations += 1;
                } else if is_invalid {
                    fp += 1;
                } else {
                    tp += 1;
                }
            }
            BarcodeMatch::NoMatch(_) => {
                if is_invalid {
                    tn += 1;
                } else if is_mutated {
                    fn_count += 1;
                    total_mutations += 1;
                } else {
                    fn_count += 1;
                }
            }
        }
    }

    let sens = metrics::sensitivity(tp, fn_count);
    let spec = metrics::specificity(tn, fp);
    let prec = metrics::precision(tp, fp);
    let rec = metrics::recall(tp, fn_count);
    let f1 = metrics::f1_score(prec, rec);
    let correction_accuracy = if total_mutations > 0 {
        correct_corrections as f64 / total_mutations as f64
    } else {
        1.0
    };

    ExtractValidationResult {
        tp,
        fp,
        tn,
        fn_count,
        sensitivity: sens,
        specificity: spec,
        precision: prec,
        recall: rec,
        f1,
        correction_accuracy,
        total_reads: dataset.r1_records.len(),
    }
}

// ─── Count validation ────────────────────────────────────────────────

/// Results from count matrix validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CountValidationResult {
    /// Pearson correlation between observed and truth
    pub pearson_r: f64,
    /// Spearman rank correlation
    pub spearman_rho: f64,
    /// Mean absolute error
    pub mae: f64,
    /// Root mean squared error
    pub rmse: f64,
    /// Fraction of cells with concordant counts
    pub cells_concordant: f64,
    /// Fraction of genes with concordant counts
    pub genes_concordant: f64,
    /// Number of cells in truth
    pub truth_cells: usize,
    /// Number of cells in observed
    pub observed_cells: usize,
}

/// Validate a count matrix against ground truth
pub fn validate_count(
    truth: &CountMatrix,
    observed: &CountMatrix,
) -> CountValidationResult {
    // Build dense vectors from both matrices for comparison
    // Align by barcode and gene names
    let (truth_vec, obs_vec) = align_matrices(truth, observed);

    let pearson_r = if truth_vec.len() > 1 {
        metrics::pearson_correlation(&truth_vec, &obs_vec)
    } else {
        0.0
    };

    let spearman_rho = if truth_vec.len() > 1 {
        metrics::spearman_correlation(&truth_vec, &obs_vec)
    } else {
        0.0
    };

    let mae = metrics::mean_absolute_error(&obs_vec, &truth_vec);
    let rmse_val = metrics::rmse(&obs_vec, &truth_vec);

    // Per-cell concordance: fraction of cells where total counts match within 10%
    let cells_concordant = compute_cell_concordance(truth, observed, 0.1);
    let genes_concordant = compute_gene_concordance(truth, observed, 0.1);

    CountValidationResult {
        pearson_r,
        spearman_rho,
        mae,
        rmse: rmse_val,
        cells_concordant,
        genes_concordant,
        truth_cells: truth.n_cols,
        observed_cells: observed.n_cols,
    }
}

/// Align two count matrices into paired vectors for correlation
fn align_matrices(truth: &CountMatrix, observed: &CountMatrix) -> (Vec<f64>, Vec<f64>) {
    // Build lookup for observed values
    let mut obs_lookup = ahash::AHashMap::new();
    let obs_barcode_idx: ahash::AHashMap<&str, usize> = observed
        .barcodes
        .iter()
        .enumerate()
        .map(|(i, b)| (b.as_str(), i))
        .collect();
    let obs_gene_idx: ahash::AHashMap<&str, usize> = observed
        .genes
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();

    for (idx, (&r, &c)) in observed.rows.iter().zip(observed.cols.iter()).enumerate() {
        obs_lookup.insert((r, c), observed.values[idx]);
    }

    let mut truth_vec = Vec::new();
    let mut obs_vec = Vec::new();

    for (t_gene_idx, gene) in truth.genes.iter().enumerate() {
        for (t_cell_idx, barcode) in truth.barcodes.iter().enumerate() {
            let t_val = truth.get(t_gene_idx, t_cell_idx) as f64;

            let o_val = if let (Some(&oi_gene), Some(&oi_cell)) =
                (obs_gene_idx.get(gene.as_str()), obs_barcode_idx.get(barcode.as_str()))
            {
                obs_lookup.get(&(oi_gene, oi_cell)).copied().unwrap_or(0) as f64
            } else {
                0.0
            };

            truth_vec.push(t_val);
            obs_vec.push(o_val);
        }
    }

    (truth_vec, obs_vec)
}

/// Compute fraction of cells where total counts are concordant within tolerance
fn compute_cell_concordance(truth: &CountMatrix, observed: &CountMatrix, tolerance: f64) -> f64 {
    let truth_counts = truth.counts_per_cell();
    let obs_barcode_map: ahash::AHashMap<&str, usize> = observed
        .barcodes
        .iter()
        .enumerate()
        .map(|(i, b)| (b.as_str(), i))
        .collect();
    let obs_counts = observed.counts_per_cell();

    let mut concordant = 0;
    let mut total = 0;

    for (i, bc) in truth.barcodes.iter().enumerate() {
        if let Some(&obs_idx) = obs_barcode_map.get(bc.as_str()) {
            total += 1;
            let t = truth_counts[i] as f64;
            let o = obs_counts[obs_idx] as f64;
            if t == 0.0 && o == 0.0 {
                concordant += 1;
            } else if t > 0.0 && ((o - t).abs() / t) <= tolerance {
                concordant += 1;
            }
        }
    }

    if total == 0 {
        return 0.0;
    }
    concordant as f64 / total as f64
}

/// Compute fraction of genes where total counts are concordant within tolerance
fn compute_gene_concordance(truth: &CountMatrix, observed: &CountMatrix, tolerance: f64) -> f64 {
    let truth_counts = truth.counts_per_gene();
    let obs_gene_map: ahash::AHashMap<&str, usize> = observed
        .genes
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();
    let obs_counts = observed.counts_per_gene();

    let mut concordant = 0;
    let mut total = 0;

    for (i, gene) in truth.genes.iter().enumerate() {
        if let Some(&obs_idx) = obs_gene_map.get(gene.as_str()) {
            total += 1;
            let t = truth_counts[i] as f64;
            let o = obs_counts[obs_idx] as f64;
            if t == 0.0 && o == 0.0 {
                concordant += 1;
            } else if t > 0.0 && ((o - t).abs() / t) <= tolerance {
                concordant += 1;
            }
        }
    }

    if total == 0 {
        return 0.0;
    }
    concordant as f64 / total as f64
}

// ─── Analysis validation ─────────────────────────────────────────────

/// Results from clustering/analysis validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisValidationResult {
    /// Adjusted Rand Index
    pub ari: f64,
    /// Normalized Mutual Information
    pub nmi: f64,
    /// Expected number of clusters
    pub n_clusters_expected: usize,
    /// Found number of clusters
    pub n_clusters_found: usize,
}

/// Validate clustering results against ground truth cell types
pub fn validate_analysis(
    truth_labels: &[usize],
    predicted_labels: &[usize],
    n_expected: usize,
) -> AnalysisValidationResult {
    let ari = metrics::adjusted_rand_index(truth_labels, predicted_labels);
    let nmi = metrics::normalized_mutual_info(truth_labels, predicted_labels);

    let n_found = predicted_labels.iter().copied().collect::<ahash::AHashSet<_>>().len();

    AnalysisValidationResult {
        ari,
        nmi,
        n_clusters_expected: n_expected,
        n_clusters_found: n_found,
    }
}

/// Run Rust-side clustering validation using label propagation
pub fn validate_clustering_with_knn(
    dataset: &SyntheticDataset,
) -> AnalysisValidationResult {
    use crate::analysis::cluster::label_propagation;

    let matrix = &dataset.truth.expression_matrix;
    let n_cells = matrix.n_cols;
    let n_genes = matrix.n_rows;

    // Build a dense representation for distance computation
    let mut dense = vec![vec![0.0f64; n_genes]; n_cells];
    for (idx, (&r, &c)) in matrix.rows.iter().zip(matrix.cols.iter()).enumerate() {
        dense[c][r] = matrix.values[idx] as f64;
    }

    // Simple normalization: log1p
    for cell in &mut dense {
        let total: f64 = cell.iter().sum();
        if total > 0.0 {
            for val in cell.iter_mut() {
                *val = (*val / total * 10000.0 + 1.0).ln();
            }
        }
    }

    // Build KNN graph (k=15)
    let k = 15.min(n_cells - 1);
    let mut graph: Vec<Vec<(usize, f64)>> = Vec::with_capacity(n_cells);

    for i in 0..n_cells {
        let mut dists: Vec<(usize, f64)> = (0..n_cells)
            .filter(|&j| j != i)
            .map(|j| {
                let dist: f64 = dense[i]
                    .iter()
                    .zip(dense[j].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                (j, dist)
            })
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        dists.truncate(k);
        graph.push(dists);
    }

    let predicted = label_propagation(&graph, 1.0, 100);

    // Get truth labels in barcode order
    let truth_labels: Vec<usize> = matrix
        .barcodes
        .iter()
        .map(|bc| dataset.truth.cell_types[bc])
        .collect();

    validate_analysis(
        &truth_labels,
        &predicted,
        dataset.config.n_cell_types,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::validation::synthetic::SyntheticConfig;

    #[test]
    fn test_extract_validation() {
        let config = SyntheticConfig {
            n_cells: 50,
            n_genes: 10,
            n_cell_types: 3,
            seed: 42,
            ..Default::default()
        };
        let dataset = SyntheticDataset::generate(config);
        let result = validate_extract(&dataset);

        // Should have high sensitivity for exact barcode matches
        assert!(
            result.sensitivity > 0.5,
            "Sensitivity too low: {}",
            result.sensitivity
        );
        assert!(result.total_reads > 0);
    }

    #[test]
    fn test_count_validation_perfect() {
        // Validating a matrix against itself should give perfect scores
        let barcodes = vec!["CELL1".to_string(), "CELL2".to_string()];
        let genes = vec!["GENE1".to_string(), "GENE2".to_string()];
        let data = vec![vec![10, 5], vec![3, 8]];

        let matrix = CountMatrix::from_dense(barcodes, genes, data);
        let result = validate_count(&matrix, &matrix);

        assert!(
            (result.pearson_r - 1.0).abs() < 1e-10,
            "Pearson should be 1.0 for identical matrices"
        );
        assert!(
            (result.cells_concordant - 1.0).abs() < 1e-10,
            "All cells should be concordant"
        );
    }

    #[test]
    fn test_analysis_validation_perfect() {
        let truth = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
        let pred = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
        let result = validate_analysis(&truth, &pred, 3);

        assert!((result.ari - 1.0).abs() < 1e-10);
        assert!((result.nmi - 1.0).abs() < 1e-10);
        assert_eq!(result.n_clusters_expected, 3);
        assert_eq!(result.n_clusters_found, 3);
    }
}
