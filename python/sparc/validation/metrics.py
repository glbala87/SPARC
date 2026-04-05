"""
Accuracy metrics for truthset validation.

Classification metrics (barcode detection), correlation metrics (expression),
and clustering metrics (ARI, NMI).
"""

from typing import Optional

import numpy as np
from scipy import stats

try:
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    _SKLEARN_AVAILABLE = True
except ImportError:
    _SKLEARN_AVAILABLE = False


def barcode_detection_metrics(
    truth_barcodes: set[str],
    detected_barcodes: set[str],
    all_observed: set[str],
) -> dict[str, float]:
    """
    Compute classification metrics for barcode detection.

    Parameters
    ----------
    truth_barcodes : set of str
        Known valid barcodes (ground truth).
    detected_barcodes : set of str
        Barcodes that were detected/accepted by the pipeline.
    all_observed : set of str
        All barcodes observed in the data (valid + invalid + mutated).

    Returns
    -------
    dict with keys: tp, fp, tn, fn, sensitivity, specificity, precision, recall, f1
    """
    tp = len(truth_barcodes & detected_barcodes)
    fp = len(detected_barcodes - truth_barcodes)
    fn_ = len(truth_barcodes - detected_barcodes)
    tn = len(all_observed - truth_barcodes - detected_barcodes)

    sensitivity = tp / (tp + fn_) if (tp + fn_) > 0 else 0.0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = sensitivity
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0

    return {
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn_,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def expression_correlation_metrics(
    truth_matrix: np.ndarray,
    observed_matrix: np.ndarray,
) -> dict[str, float]:
    """
    Compute correlation metrics between truth and observed expression matrices.

    Parameters
    ----------
    truth_matrix : np.ndarray
        Ground truth count matrix (cells x genes).
    observed_matrix : np.ndarray
        Observed count matrix (cells x genes), same shape as truth.

    Returns
    -------
    dict with keys: pearson_r, spearman_rho, mae, rmse
    """
    truth_flat = truth_matrix.ravel().astype(np.float64)
    obs_flat = observed_matrix.ravel().astype(np.float64)

    if len(truth_flat) < 2 or np.std(truth_flat) == 0 or np.std(obs_flat) == 0:
        pearson_r = 0.0
        spearman_rho = 0.0
    else:
        pearson_r, _ = stats.pearsonr(truth_flat, obs_flat)
        spearman_rho, _ = stats.spearmanr(truth_flat, obs_flat)

    mae = float(np.mean(np.abs(truth_flat - obs_flat)))
    rmse = float(np.sqrt(np.mean((truth_flat - obs_flat) ** 2)))

    return {
        "pearson_r": float(pearson_r),
        "spearman_rho": float(spearman_rho),
        "mae": mae,
        "rmse": rmse,
    }


def clustering_metrics(
    truth_labels: np.ndarray,
    predicted_labels: np.ndarray,
) -> dict[str, float]:
    """
    Compute clustering accuracy metrics.

    Parameters
    ----------
    truth_labels : array-like
        Ground truth cluster/cell-type labels.
    predicted_labels : array-like
        Predicted cluster labels.

    Returns
    -------
    dict with keys: ari, nmi, n_clusters_expected, n_clusters_found
    """
    if not _SKLEARN_AVAILABLE:
        raise ImportError("scikit-learn is required for clustering metrics: pip install scikit-learn")

    truth_labels = np.asarray(truth_labels)
    predicted_labels = np.asarray(predicted_labels)

    ari = float(adjusted_rand_score(truth_labels, predicted_labels))
    nmi = float(normalized_mutual_info_score(truth_labels, predicted_labels))

    n_expected = len(np.unique(truth_labels))
    n_found = len(np.unique(predicted_labels))

    return {
        "ari": ari,
        "nmi": nmi,
        "n_clusters_expected": n_expected,
        "n_clusters_found": n_found,
    }
