//! Python bindings for validation module

use pyo3::prelude::*;
use sparc_core::validation::metrics;

/// Compute Adjusted Rand Index between two label vectors.
#[pyfunction]
pub fn adjusted_rand_index(labels_true: Vec<usize>, labels_pred: Vec<usize>) -> f64 {
    metrics::adjusted_rand_index(&labels_true, &labels_pred)
}

/// Compute Normalized Mutual Information between two label vectors.
#[pyfunction]
pub fn normalized_mutual_info(labels_true: Vec<usize>, labels_pred: Vec<usize>) -> f64 {
    metrics::normalized_mutual_info(&labels_true, &labels_pred)
}

/// Compute Pearson correlation between two vectors.
#[pyfunction]
pub fn pearson_correlation(x: Vec<f64>, y: Vec<f64>) -> f64 {
    metrics::pearson_correlation(&x, &y)
}

/// Compute Spearman rank correlation between two vectors.
#[pyfunction]
pub fn spearman_correlation(x: Vec<f64>, y: Vec<f64>) -> f64 {
    metrics::spearman_correlation(&x, &y)
}

/// Compute F1 score from precision and recall.
#[pyfunction]
pub fn f1_score(precision: f64, recall: f64) -> f64 {
    metrics::f1_score(precision, recall)
}
