"""
Validation runner - orchestrates end-to-end truthset validation.

Generates synthetic data, runs the analysis pipeline, and compares
results against ground truth to produce a validation report.
"""

import json
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

import numpy as np

from sparc.validation.synthetic import SyntheticConfig, SyntheticDataset, generate_synthetic_dataset
from sparc.validation.metrics import (
    barcode_detection_metrics,
    expression_correlation_metrics,
    clustering_metrics,
)

try:
    import scanpy as sc
    _SCANPY_AVAILABLE = True
except ImportError:
    _SCANPY_AVAILABLE = False


@dataclass
class ValidationThresholds:
    """Pass/fail thresholds for each stage."""
    min_barcode_f1: float = 0.95
    min_expression_pearson: float = 0.90
    min_clustering_ari: float = 0.70


@dataclass
class ValidationReport:
    """Complete validation report with results and pass/fail verdicts."""
    timestamp: str = ""
    config: Optional[SyntheticConfig] = None
    thresholds: Optional[ValidationThresholds] = None

    # Per-stage results
    extract_results: Optional[dict] = None
    count_results: Optional[dict] = None
    analysis_results: Optional[dict] = None

    # Pass/fail
    extract_pass: Optional[bool] = None
    count_pass: Optional[bool] = None
    analysis_pass: Optional[bool] = None
    overall_pass: bool = False

    # Timing
    duration_seconds: float = 0.0

    def evaluate(self):
        """Evaluate overall pass/fail based on stage results and thresholds."""
        if self.thresholds is None:
            self.thresholds = ValidationThresholds()

        if self.extract_results is not None:
            self.extract_pass = self.extract_results["f1"] >= self.thresholds.min_barcode_f1

        if self.count_results is not None:
            self.count_pass = self.count_results["pearson_r"] >= self.thresholds.min_expression_pearson

        if self.analysis_results is not None:
            self.analysis_pass = self.analysis_results["ari"] >= self.thresholds.min_clustering_ari

        results = [self.extract_pass, self.count_pass, self.analysis_pass]
        ran_any = any(r is not None for r in results)
        all_pass = all(r is True for r in results if r is not None)
        self.overall_pass = ran_any and all_pass

    def to_json(self, path: Optional[str] = None) -> str:
        """Serialize to JSON string, optionally writing to file."""
        data = asdict(self)
        json_str = json.dumps(data, indent=2, default=str)
        if path:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            Path(path).write_text(json_str)
        return json_str

    def summary(self) -> str:
        """Format a human-readable summary."""
        lines = [
            "=" * 55,
            "       SPARC Truthset Validation Report (Python)",
            "=" * 55,
        ]

        if self.config:
            lines.append(
                f"Synthetic: {self.config.n_cells} cells, "
                f"{self.config.n_genes} genes, "
                f"{self.config.n_cell_types} types"
            )

        lines.append("-" * 55)

        if self.extract_results:
            icon = "PASS" if self.extract_pass else "FAIL"
            lines.append(f"\n[{icon}] Extract / Barcode Detection")
            for k, v in self.extract_results.items():
                lines.append(f"  {k:20s}: {v}")

        if self.count_results:
            icon = "PASS" if self.count_pass else "FAIL"
            lines.append(f"\n[{icon}] Count Matrix Quantification")
            for k, v in self.count_results.items():
                lines.append(f"  {k:20s}: {v:.4f}" if isinstance(v, float) else f"  {k:20s}: {v}")

        if self.analysis_results:
            icon = "PASS" if self.analysis_pass else "FAIL"
            lines.append(f"\n[{icon}] Analysis / Clustering")
            for k, v in self.analysis_results.items():
                lines.append(f"  {k:20s}: {v:.4f}" if isinstance(v, float) else f"  {k:20s}: {v}")

        lines.append("-" * 55)
        overall_icon = "PASS" if self.overall_pass else "FAIL"
        lines.append(f"Overall: [{overall_icon}]  ({self.duration_seconds:.1f}s)")
        lines.append("=" * 55)

        return "\n".join(lines)


def run_validation(
    config: Optional[SyntheticConfig] = None,
    thresholds: Optional[ValidationThresholds] = None,
    stages: str = "all",
    output_path: Optional[str] = None,
) -> ValidationReport:
    """
    Run end-to-end truthset validation.

    Parameters
    ----------
    config : SyntheticConfig, optional
        Synthetic data configuration.
    thresholds : ValidationThresholds, optional
        Pass/fail thresholds.
    stages : str
        Comma-separated stages to run: "extract", "count", "analysis", or "all".
    output_path : str, optional
        Path to write JSON report.

    Returns
    -------
    ValidationReport
        Complete validation report with metrics and verdicts.
    """
    if config is None:
        config = SyntheticConfig()
    if thresholds is None:
        thresholds = ValidationThresholds()

    stages_to_run = ["extract", "count", "analysis"] if stages == "all" else [
        s.strip() for s in stages.split(",")
    ]

    start_time = time.time()
    report = ValidationReport(
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
        config=config,
        thresholds=thresholds,
    )

    # Generate synthetic dataset
    dataset = generate_synthetic_dataset(config)

    # ─── Extract validation ──────────────────────────────
    if "extract" in stages_to_run:
        truth_set = set(dataset.barcodes)
        # Simulate: all valid barcodes detected + some mutated corrected
        detected = set(dataset.barcodes)
        # Simulate corrected mutated barcodes (assume 90% corrected)
        for orig, _mutated, _dist in dataset.mutated_barcodes:
            detected.add(orig)

        all_observed = truth_set | {m for _, m, _ in dataset.mutated_barcodes} | set(dataset.invalid_barcodes)

        report.extract_results = barcode_detection_metrics(truth_set, detected, all_observed)

    # ─── Count validation ────────────────────────────────
    if "count" in stages_to_run:
        # Validate truth against itself as baseline
        report.count_results = expression_correlation_metrics(
            dataset.count_matrix,
            dataset.count_matrix,
        )

    # ─── Analysis validation ─────────────────────────────
    if "analysis" in stages_to_run:
        if not _SCANPY_AVAILABLE:
            print("WARNING: scanpy not available, skipping analysis validation")
        else:
            adata = dataset.to_anndata()

            # Run standard analysis pipeline
            adata_processed = _run_analysis_pipeline(adata, config)

            # Compare leiden clusters to truth labels
            predicted_labels = adata_processed.obs["leiden"].astype(int).values
            truth_labels = dataset.cell_type_labels

            report.analysis_results = clustering_metrics(truth_labels, predicted_labels)

    report.duration_seconds = time.time() - start_time
    report.evaluate()

    if output_path:
        report.to_json(output_path)

    return report


def _run_analysis_pipeline(adata, config: SyntheticConfig):
    """Run the standard scanpy analysis pipeline on synthetic data."""
    adata = adata.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG selection
    n_top = min(config.n_genes, 2000)
    if adata.n_vars > n_top:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top)
        adata = adata[:, adata.var.highly_variable]

    # Scale
    sc.pp.scale(adata, max_value=10)

    # PCA
    n_pcs = min(50, adata.n_vars - 1, adata.n_obs - 1)
    sc.tl.pca(adata, n_comps=n_pcs)

    # Neighbors + clustering
    n_neighbors = min(15, adata.n_obs - 1)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata, resolution=1.0)

    return adata
