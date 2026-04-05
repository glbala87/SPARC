//! Python bindings for analysis modules (normalize, PCA, clustering)

use numpy::{PyArray1, PyArray2, ToPyArray};
use pyo3::prelude::*;

use sparc_core::analysis::{
    cluster::label_propagation,
    neighbors::build_knn_graph,
    normalize::{highly_variable_genes, log1p_transform, normalize_total, scale},
    pca::pca,
};

/// Normalize cells to a target sum, apply log1p, and return processed matrix.
#[pyfunction]
#[pyo3(signature = (data, target_sum = 10000.0))]
pub fn py_normalize_total<'py>(py: Python<'py>, data: Vec<Vec<f64>>, target_sum: f64) -> &'py PyArray2<f64> {
    let mut d = data;
    normalize_total(&mut d, target_sum);
    log1p_transform(&mut d);
    vec_to_array2(py, &d)
}

/// Scale data to zero mean and unit variance per gene.
#[pyfunction]
#[pyo3(signature = (data, max_value = None))]
pub fn py_scale<'py>(py: Python<'py>, data: Vec<Vec<f64>>, max_value: Option<f64>) -> &'py PyArray2<f64> {
    let mut d = data;
    scale(&mut d, max_value);
    vec_to_array2(py, &d)
}

/// Find highly variable genes. Returns list of gene indices.
#[pyfunction]
#[pyo3(signature = (data, min_mean = 0.0125, max_mean = 3.0, min_disp = 0.5))]
pub fn py_highly_variable_genes(data: Vec<Vec<f64>>, min_mean: f64, max_mean: f64, min_disp: f64) -> Vec<usize> {
    highly_variable_genes(&data, min_mean, max_mean, min_disp)
}

/// Run PCA. Returns (scores: n_cells x n_pcs, variances: n_pcs).
#[pyfunction]
#[pyo3(signature = (data, n_components = 50))]
pub fn py_pca<'py>(py: Python<'py>, data: Vec<Vec<f64>>, n_components: usize) -> (&'py PyArray2<f64>, &'py PyArray1<f64>) {
    let (scores, variances) = pca(&data, n_components);
    (vec_to_array2(py, &scores), variances.to_pyarray(py))
}

/// Build KNN graph. Returns list of (neighbor_index, distance) per cell.
#[pyfunction]
#[pyo3(signature = (coords, k = 15))]
pub fn py_build_knn_graph(coords: Vec<Vec<f64>>, k: usize) -> Vec<Vec<(usize, f64)>> {
    build_knn_graph(&coords, k)
}

/// Run label propagation clustering. Returns cluster labels.
#[pyfunction]
#[pyo3(signature = (graph, resolution = 1.0, max_iterations = 100))]
pub fn py_label_propagation(graph: Vec<Vec<(usize, f64)>>, resolution: f64, max_iterations: usize) -> Vec<usize> {
    label_propagation(&graph, resolution, max_iterations)
}

/// Run full analysis pipeline: normalize -> HVG -> scale -> PCA -> KNN -> cluster.
/// Returns (pca_scores, cluster_labels, hvg_indices).
#[pyfunction]
#[pyo3(signature = (data, n_pcs = 50, n_neighbors = 15, resolution = 1.0))]
pub fn py_run_analysis<'py>(
    py: Python<'py>,
    data: Vec<Vec<f64>>,
    n_pcs: usize,
    n_neighbors: usize,
    resolution: f64,
) -> (&'py PyArray2<f64>, Vec<usize>, Vec<usize>) {
    let mut d = data;
    let n_cells = d.len();

    normalize_total(&mut d, 10000.0);
    log1p_transform(&mut d);

    let hvg = highly_variable_genes(&d, 0.0125, 3.0, 0.5);
    let mut analysis_data = if hvg.is_empty() {
        d
    } else {
        d.iter().map(|cell| hvg.iter().map(|&g| cell[g]).collect()).collect()
    };

    scale(&mut analysis_data, Some(10.0));

    let actual_pcs = n_pcs.min(analysis_data[0].len().saturating_sub(1)).min(n_cells.saturating_sub(1));
    let (pca_coords, _) = pca(&analysis_data, actual_pcs);

    let k = n_neighbors.min(n_cells.saturating_sub(1));
    let graph = build_knn_graph(&pca_coords, k);
    let labels = label_propagation(&graph, resolution, 100);

    (vec_to_array2(py, &pca_coords), labels, hvg)
}

fn vec_to_array2<'py>(py: Python<'py>, data: &[Vec<f64>]) -> &'py PyArray2<f64> {
    if data.is_empty() {
        return PyArray2::zeros(py, (0, 0), false);
    }
    let n_rows = data.len();
    let n_cols = data[0].len();
    let flat: Vec<f64> = data.iter().flat_map(|r| r.iter().copied()).collect();
    PyArray1::from_vec(py, flat)
        .reshape((n_rows, n_cols))
        .expect("reshape dimensions match")
}
