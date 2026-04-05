"""Tests for the SPARC truthset validation framework."""

import numpy as np
import pytest


class TestSyntheticGeneration:
    """Tests for synthetic data generation."""

    def test_generate_default_config(self):
        from sparc.validation.synthetic import SyntheticConfig, generate_synthetic_dataset

        config = SyntheticConfig(n_cells=50, n_genes=20, n_cell_types=3, seed=42)
        dataset = generate_synthetic_dataset(config)

        assert dataset.count_matrix.shape == (50, 20)
        assert len(dataset.barcodes) == 50
        assert len(dataset.genes) == 20
        assert len(dataset.cell_type_names) == 3
        assert len(dataset.cell_type_labels) == 50

    def test_barcodes_unique(self):
        from sparc.validation.synthetic import SyntheticConfig, generate_synthetic_dataset

        config = SyntheticConfig(n_cells=100, seed=42)
        dataset = generate_synthetic_dataset(config)
        assert len(set(dataset.barcodes)) == 100

    def test_deterministic_output(self):
        from sparc.validation.synthetic import SyntheticConfig, generate_synthetic_dataset

        config = SyntheticConfig(n_cells=30, n_genes=10, seed=77)
        d1 = generate_synthetic_dataset(config)
        d2 = generate_synthetic_dataset(config)

        assert d1.barcodes == d2.barcodes
        assert d1.genes == d2.genes
        np.testing.assert_array_equal(d1.count_matrix, d2.count_matrix)

    def test_cell_type_distribution(self):
        from sparc.validation.synthetic import SyntheticConfig, generate_synthetic_dataset

        config = SyntheticConfig(n_cells=100, n_cell_types=5, seed=42)
        dataset = generate_synthetic_dataset(config)

        # Round-robin: each type should have 20 cells
        for ct in range(5):
            count = np.sum(dataset.cell_type_labels == ct)
            assert count == 20, f"CellType_{ct} should have 20 cells, got {count}"

    def test_marker_genes_enriched(self):
        from sparc.validation.synthetic import SyntheticConfig, generate_synthetic_dataset

        config = SyntheticConfig(
            n_cells=200,
            n_genes=50,
            n_cell_types=3,
            n_markers_per_type=5,
            marker_fold_change=10.0,
            base_expression=2.0,
            seed=42,
        )
        dataset = generate_synthetic_dataset(config)

        for ct in range(3):
            profile = dataset.expression_profiles[ct]
            marker_start = ct * 5
            marker_mean = np.mean(profile[marker_start:marker_start + 5])
            non_markers = np.concatenate([profile[:marker_start], profile[marker_start + 5:]])
            nonmarker_mean = np.mean(non_markers)
            assert marker_mean > nonmarker_mean * 3, (
                f"Marker genes for CellType_{ct} should be enriched: "
                f"marker_mean={marker_mean:.2f}, nonmarker_mean={nonmarker_mean:.2f}"
            )

    def test_mutated_barcodes(self):
        from sparc.validation.synthetic import SyntheticConfig, generate_synthetic_dataset

        config = SyntheticConfig(n_cells=100, mutation_rate=0.1, seed=42)
        dataset = generate_synthetic_dataset(config)

        assert len(dataset.mutated_barcodes) == 10
        for orig, mutated, dist in dataset.mutated_barcodes:
            assert dist == 1
            assert len(orig) == len(mutated)
            diffs = sum(a != b for a, b in zip(orig, mutated))
            assert diffs == 1, f"Expected 1 difference, got {diffs}"

    def test_invalid_barcodes(self):
        from sparc.validation.synthetic import SyntheticConfig, generate_synthetic_dataset

        config = SyntheticConfig(n_cells=100, invalid_barcode_rate=0.05, seed=42)
        dataset = generate_synthetic_dataset(config)

        assert len(dataset.invalid_barcodes) == 5
        barcode_set = set(dataset.barcodes)
        for inv in dataset.invalid_barcodes:
            assert inv not in barcode_set

    def test_to_anndata(self):
        from sparc.validation.synthetic import SyntheticConfig, generate_synthetic_dataset

        config = SyntheticConfig(n_cells=50, n_genes=20, n_cell_types=3, seed=42)
        dataset = generate_synthetic_dataset(config)

        adata = dataset.to_anndata()
        assert adata.n_obs == 50
        assert adata.n_vars == 20
        assert "cell_type" in adata.obs.columns
        assert "cell_type_id" in adata.obs.columns


class TestMetrics:
    """Tests for accuracy metrics."""

    def test_barcode_detection_perfect(self):
        from sparc.validation.metrics import barcode_detection_metrics

        truth = {"A", "B", "C"}
        detected = {"A", "B", "C"}
        all_obs = {"A", "B", "C", "D", "E"}

        result = barcode_detection_metrics(truth, detected, all_obs)
        assert result["tp"] == 3
        assert result["fp"] == 0
        assert result["fn"] == 0
        assert result["f1"] == 1.0
        assert result["sensitivity"] == 1.0
        assert result["precision"] == 1.0

    def test_barcode_detection_partial(self):
        from sparc.validation.metrics import barcode_detection_metrics

        truth = {"A", "B", "C"}
        detected = {"A", "B", "D"}  # missed C, falsely detected D
        all_obs = {"A", "B", "C", "D", "E"}

        result = barcode_detection_metrics(truth, detected, all_obs)
        assert result["tp"] == 2
        assert result["fp"] == 1
        assert result["fn"] == 1
        assert result["precision"] == 2 / 3
        assert result["recall"] == 2 / 3

    def test_expression_correlation_perfect(self):
        from sparc.validation.metrics import expression_correlation_metrics

        matrix = np.array([[1, 2, 3], [4, 5, 6]])
        result = expression_correlation_metrics(matrix, matrix)

        assert abs(result["pearson_r"] - 1.0) < 1e-10
        assert abs(result["mae"]) < 1e-10
        assert abs(result["rmse"]) < 1e-10

    def test_expression_correlation_different(self):
        from sparc.validation.metrics import expression_correlation_metrics

        truth = np.array([[1, 2, 3], [4, 5, 6]], dtype=float)
        obs = np.array([[2, 3, 4], [5, 6, 7]], dtype=float)
        result = expression_correlation_metrics(truth, obs)

        assert abs(result["pearson_r"] - 1.0) < 1e-10  # Same trend
        assert abs(result["mae"] - 1.0) < 1e-10  # Offset by 1

    def test_clustering_metrics_perfect(self):
        from sparc.validation.metrics import clustering_metrics

        truth = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        pred = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        result = clustering_metrics(truth, pred)

        assert abs(result["ari"] - 1.0) < 1e-10
        assert abs(result["nmi"] - 1.0) < 1e-10
        assert result["n_clusters_expected"] == 3
        assert result["n_clusters_found"] == 3

    def test_clustering_metrics_permuted(self):
        """ARI/NMI should be 1.0 even with relabeled clusters."""
        from sparc.validation.metrics import clustering_metrics

        truth = np.array([0, 0, 0, 1, 1, 1])
        pred = np.array([1, 1, 1, 0, 0, 0])
        result = clustering_metrics(truth, pred)

        assert abs(result["ari"] - 1.0) < 1e-10

    def test_clustering_metrics_random(self):
        from sparc.validation.metrics import clustering_metrics

        rng = np.random.default_rng(42)
        truth = rng.integers(0, 5, size=100)
        pred = rng.integers(0, 5, size=100)
        result = clustering_metrics(truth, pred)

        assert abs(result["ari"]) < 0.3, f"ARI for random labels should be near 0: {result['ari']}"


class TestValidationRunner:
    """Tests for the validation runner."""

    def test_full_validation(self):
        from sparc.validation.runner import run_validation, ValidationThresholds
        from sparc.validation.synthetic import SyntheticConfig

        config = SyntheticConfig(n_cells=50, n_genes=20, n_cell_types=3, seed=42)
        thresholds = ValidationThresholds(
            min_barcode_f1=0.5,
            min_expression_pearson=0.5,
            min_clustering_ari=0.0,
        )

        report = run_validation(config=config, thresholds=thresholds, stages="extract,count")

        assert report.extract_results is not None
        assert report.count_results is not None
        assert report.extract_pass is not None
        assert report.count_pass is not None

    def test_count_self_validation(self):
        """Count validation against itself should give perfect scores."""
        from sparc.validation.runner import run_validation
        from sparc.validation.synthetic import SyntheticConfig

        config = SyntheticConfig(n_cells=30, n_genes=10, seed=42)
        report = run_validation(config=config, stages="count")

        assert report.count_results is not None
        assert abs(report.count_results["pearson_r"] - 1.0) < 1e-10
        assert report.count_pass is True

    def test_report_json_serialization(self):
        from sparc.validation.runner import run_validation
        from sparc.validation.synthetic import SyntheticConfig

        import json

        config = SyntheticConfig(n_cells=20, n_genes=10, seed=42)
        report = run_validation(config=config, stages="extract,count")

        json_str = report.to_json()
        parsed = json.loads(json_str)
        assert "extract_results" in parsed
        assert "count_results" in parsed
        assert "overall_pass" in parsed

    def test_report_summary(self):
        from sparc.validation.runner import run_validation
        from sparc.validation.synthetic import SyntheticConfig

        config = SyntheticConfig(n_cells=20, n_genes=10, seed=42)
        report = run_validation(config=config, stages="extract,count")

        summary = report.summary()
        assert "SPARC Truthset Validation Report" in summary
        assert "Extract" in summary
        assert "Count" in summary
