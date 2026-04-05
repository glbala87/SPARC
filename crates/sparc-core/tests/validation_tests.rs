//! Integration tests for the truthset validation framework

use sparc_core::validation::metrics;
use sparc_core::validation::report::{ValidationReport, ValidationThresholds};
use sparc_core::validation::stages;
use sparc_core::validation::synthetic::{SyntheticConfig, SyntheticDataset};

// ─── Synthetic data generation tests ─────────────────────────────────

#[test]
fn test_synthetic_dataset_dimensions() {
    let config = SyntheticConfig {
        n_cells: 100,
        n_genes: 50,
        n_cell_types: 4,
        n_markers_per_type: 5,
        seed: 42,
        ..Default::default()
    };
    let dataset = SyntheticDataset::generate(config);

    assert_eq!(dataset.truth.barcodes.len(), 100);
    assert_eq!(dataset.truth.genes.len(), 50);
    assert_eq!(dataset.truth.cell_type_names.len(), 4);
    assert_eq!(dataset.truth.expression_matrix.n_rows, 50); // genes
    assert_eq!(dataset.truth.expression_matrix.n_cols, 100); // cells
    assert_eq!(dataset.truth.expression_profiles.len(), 4);
    assert!(dataset.truth.expression_profiles[0].len() == 50);
}

#[test]
fn test_synthetic_barcodes_unique() {
    let config = SyntheticConfig {
        n_cells: 200,
        seed: 99,
        ..Default::default()
    };
    let dataset = SyntheticDataset::generate(config);

    let set: ahash::AHashSet<&str> = dataset.truth.barcodes.iter().map(|s| s.as_str()).collect();
    assert_eq!(set.len(), 200, "All barcodes should be unique");
}

#[test]
fn test_synthetic_deterministic_output() {
    let config = SyntheticConfig {
        n_cells: 30,
        n_genes: 15,
        seed: 77,
        ..Default::default()
    };
    let d1 = SyntheticDataset::generate(config.clone());
    let d2 = SyntheticDataset::generate(config);

    assert_eq!(d1.truth.barcodes, d2.truth.barcodes);
    assert_eq!(d1.truth.genes, d2.truth.genes);
    assert_eq!(d1.truth.expression_matrix.values.len(), d2.truth.expression_matrix.values.len());
    assert_eq!(d1.r1_records.len(), d2.r1_records.len());
}

#[test]
fn test_synthetic_mutation_count() {
    let config = SyntheticConfig {
        n_cells: 100,
        mutation_rate: 0.1,
        seed: 42,
        ..Default::default()
    };
    let dataset = SyntheticDataset::generate(config);

    assert_eq!(dataset.truth.mutated_barcodes.len(), 10);
    for (_orig, _mutated, dist) in &dataset.truth.mutated_barcodes {
        assert_eq!(*dist, 1, "All mutations should be Hamming distance 1");
    }
}

#[test]
fn test_synthetic_invalid_barcodes() {
    let config = SyntheticConfig {
        n_cells: 100,
        invalid_barcode_rate: 0.05,
        seed: 42,
        ..Default::default()
    };
    let dataset = SyntheticDataset::generate(config);

    assert_eq!(dataset.truth.invalid_barcodes.len(), 5);
    for inv in &dataset.truth.invalid_barcodes {
        assert!(
            !dataset.whitelist.contains(inv),
            "Invalid barcode should not be in whitelist"
        );
    }
}

// ─── Metric unit tests ──────────────────────────────────────────────

#[test]
fn test_classification_edge_cases() {
    assert_eq!(metrics::sensitivity(0, 0), 0.0);
    assert_eq!(metrics::specificity(0, 0), 0.0);
    assert_eq!(metrics::precision(0, 0), 0.0);
    assert_eq!(metrics::f1_score(0.0, 0.0), 0.0);

    assert!((metrics::sensitivity(100, 0) - 1.0).abs() < 1e-10);
    assert!((metrics::precision(100, 0) - 1.0).abs() < 1e-10);
}

#[test]
fn test_pearson_constant_vectors() {
    let x = vec![5.0, 5.0, 5.0, 5.0];
    let y = vec![3.0, 3.0, 3.0, 3.0];
    // Constant vectors -> correlation undefined, should return 0
    assert_eq!(metrics::pearson_correlation(&x, &y), 0.0);
}

#[test]
fn test_spearman_with_ties() {
    let x = vec![1.0, 1.0, 2.0, 3.0];
    let y = vec![1.0, 2.0, 3.0, 4.0];
    let rho = metrics::spearman_correlation(&x, &y);
    assert!(rho > 0.8, "Should have strong positive correlation: {}", rho);
}

#[test]
fn test_ari_random_labels() {
    // Random-ish assignments should give ARI near 0
    let truth = vec![0, 0, 1, 1, 2, 2, 0, 1, 2, 0];
    let pred = vec![1, 2, 0, 2, 1, 0, 2, 0, 1, 2];
    let ari = metrics::adjusted_rand_index(&truth, &pred);
    assert!(ari.abs() < 0.5, "ARI for random labels should be near 0: {}", ari);
}

#[test]
fn test_nmi_range() {
    let truth = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
    let pred = vec![0, 0, 1, 1, 1, 2, 2, 2, 0];
    let nmi = metrics::normalized_mutual_info(&truth, &pred);
    assert!(nmi >= 0.0 && nmi <= 1.0, "NMI should be in [0, 1]: {}", nmi);
}

// ─── Stage validation tests ─────────────────────────────────────────

#[test]
fn test_extract_validation_runs() {
    let config = SyntheticConfig {
        n_cells: 30,
        n_genes: 10,
        n_cell_types: 2,
        seed: 42,
        ..Default::default()
    };
    let dataset = SyntheticDataset::generate(config);
    let result = stages::validate_extract(&dataset);

    assert!(result.total_reads > 0, "Should process some reads");
    assert!(result.tp > 0, "Should have true positives");
    assert!(result.sensitivity >= 0.0 && result.sensitivity <= 1.0);
    assert!(result.f1 >= 0.0 && result.f1 <= 1.0);
}

#[test]
fn test_count_validation_self_comparison() {
    let config = SyntheticConfig {
        n_cells: 20,
        n_genes: 10,
        seed: 42,
        ..Default::default()
    };
    let dataset = SyntheticDataset::generate(config);

    let result = stages::validate_count(
        &dataset.truth.expression_matrix,
        &dataset.truth.expression_matrix,
    );

    assert!(
        (result.pearson_r - 1.0).abs() < 1e-10,
        "Self-comparison Pearson should be 1.0: {}",
        result.pearson_r
    );
    assert!(
        result.mae.abs() < 1e-10,
        "Self-comparison MAE should be 0: {}",
        result.mae
    );
}

#[test]
fn test_analysis_validation_runs() {
    let config = SyntheticConfig {
        n_cells: 60,
        n_genes: 30,
        n_cell_types: 3,
        n_markers_per_type: 8,
        marker_fold_change: 10.0,
        seed: 42,
        ..Default::default()
    };
    let dataset = SyntheticDataset::generate(config);
    let result = stages::validate_clustering_with_knn(&dataset);

    assert!(result.ari >= -1.0 && result.ari <= 1.0);
    assert!(result.nmi >= 0.0 && result.nmi <= 1.0);
    assert_eq!(result.n_clusters_expected, 3);
    assert!(result.n_clusters_found > 0);
}

// ─── Report tests ────────────────────────────────────────────────────

#[test]
fn test_full_validation_pipeline() {
    let config = SyntheticConfig {
        n_cells: 40,
        n_genes: 15,
        n_cell_types: 2,
        n_markers_per_type: 5,
        marker_fold_change: 8.0,
        seed: 42,
        ..Default::default()
    };
    let thresholds = ValidationThresholds {
        min_barcode_f1: 0.5,  // Lenient for testing
        min_expression_pearson: 0.5,
        min_clustering_ari: 0.0,
    };

    let dataset = SyntheticDataset::generate(config.clone());
    let mut report = ValidationReport::new(config, thresholds);

    // Run all stages
    let extract_result = stages::validate_extract(&dataset);
    report.set_extract_results(extract_result);

    let count_result = stages::validate_count(
        &dataset.truth.expression_matrix,
        &dataset.truth.expression_matrix,
    );
    report.set_count_results(count_result);

    let analysis_result = stages::validate_clustering_with_knn(&dataset);
    report.set_analysis_results(analysis_result);

    // Verify report structure
    assert!(report.extract_results.is_some());
    assert!(report.count_results.is_some());
    assert!(report.analysis_results.is_some());
    assert!(report.stage_results.extract.is_some());
    assert!(report.stage_results.count.is_some());
    assert!(report.stage_results.analysis.is_some());

    // Should be valid JSON
    let json = report.to_json().unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
    assert!(parsed.is_object());
    assert!(parsed["extract_results"].is_object());
    assert!(parsed["count_results"].is_object());
    assert!(parsed["analysis_results"].is_object());

    // Summary should be non-empty
    let summary = report.summary();
    assert!(summary.contains("SPARC Truthset Validation Report"));
    assert!(summary.contains("Extract"));
    assert!(summary.contains("Count"));
    assert!(summary.contains("Analysis"));
}
