"""
SPARC Truthset Validation Framework

Generates synthetic ground-truth datasets and validates pipeline stages
(extract, count, analysis) against expected results with accuracy metrics.
"""

from sparc.validation.synthetic import SyntheticDataset, SyntheticConfig, generate_synthetic_dataset
from sparc.validation.metrics import (
    barcode_detection_metrics,
    expression_correlation_metrics,
    clustering_metrics,
)
from sparc.validation.runner import run_validation, ValidationReport

__all__ = [
    "SyntheticDataset",
    "SyntheticConfig",
    "generate_synthetic_dataset",
    "barcode_detection_metrics",
    "expression_correlation_metrics",
    "clustering_metrics",
    "run_validation",
    "ValidationReport",
]
