//! Count matrix Python bindings

use numpy::{PyArray1, PyArray2, ToPyArray};
use pyo3::prelude::*;
use sparc_core::count::{CountMatrix, GeneCounter};

/// Python wrapper for CountMatrix
#[pyclass(name = "CountMatrix")]
pub struct PyCountMatrix {
    inner: CountMatrix,
}

#[pymethods]
impl PyCountMatrix {
    /// Create empty count matrix
    #[new]
    fn new() -> Self {
        Self {
            inner: CountMatrix::new(),
        }
    }

    /// Get barcodes (column names)
    #[getter]
    fn barcodes(&self) -> Vec<String> {
        self.inner.barcodes.clone()
    }

    /// Get genes (row names)
    #[getter]
    fn genes(&self) -> Vec<String> {
        self.inner.genes.clone()
    }

    /// Get number of rows (genes)
    #[getter]
    fn n_rows(&self) -> usize {
        self.inner.n_rows
    }

    /// Get number of columns (cells)
    #[getter]
    fn n_cols(&self) -> usize {
        self.inner.n_cols
    }

    /// Get number of non-zero entries
    #[getter]
    fn nnz(&self) -> usize {
        self.inner.values.len()
    }

    /// Get shape as tuple
    #[getter]
    fn shape(&self) -> (usize, usize) {
        (self.inner.n_rows, self.inner.n_cols)
    }

    /// Get row indices as numpy array
    fn row_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<usize>> {
        self.inner.rows.to_pyarray(py)
    }

    /// Get column indices as numpy array
    fn col_indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<usize>> {
        self.inner.cols.to_pyarray(py)
    }

    /// Get values as numpy array
    fn values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        self.inner.values.to_pyarray(py)
    }

    /// Get total counts per cell as numpy array
    fn counts_per_cell<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        self.inner.counts_per_cell().to_pyarray(py)
    }

    /// Get total counts per gene as numpy array
    fn counts_per_gene<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        self.inner.counts_per_gene().to_pyarray(py)
    }

    /// Get number of genes detected per cell
    fn genes_per_cell<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        self.inner.genes_per_cell().to_pyarray(py)
    }

    /// Get number of cells expressing each gene
    fn cells_per_gene<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        self.inner.cells_per_gene().to_pyarray(py)
    }

    /// Convert to dense numpy array (for small matrices)
    fn to_dense<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        let mut dense = vec![vec![0u32; self.inner.n_cols]; self.inner.n_rows];
        for (i, ((&r, &c), &v)) in self
            .inner
            .rows
            .iter()
            .zip(self.inner.cols.iter())
            .zip(self.inner.values.iter())
            .enumerate()
        {
            let _ = i; // suppress warning
            dense[r][c] = v;
        }

        let flat: Vec<u32> = dense.into_iter().flatten().collect();
        PyArray1::from_vec(py, flat)
            .reshape((self.inner.n_rows, self.inner.n_cols))
            .unwrap()
    }

    /// Write to Matrix Market format
    fn write_mtx(&self, path: &str) -> PyResult<()> {
        self.inner
            .write_mtx(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
    }

    /// Write barcodes to file
    fn write_barcodes(&self, path: &str) -> PyResult<()> {
        self.inner
            .write_barcodes(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
    }

    /// Write genes to file
    fn write_genes(&self, path: &str) -> PyResult<()> {
        self.inner
            .write_genes(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "CountMatrix(genes={}, cells={}, nnz={})",
            self.inner.n_rows,
            self.inner.n_cols,
            self.inner.values.len()
        )
    }
}

/// Python wrapper for GeneCounter
#[pyclass(name = "GeneCounter")]
pub struct PyGeneCounter {
    inner: GeneCounter,
}

#[pymethods]
impl PyGeneCounter {
    #[new]
    fn new() -> Self {
        Self {
            inner: GeneCounter::new(),
        }
    }

    /// Add a count for a barcode-gene pair
    fn add_count(&mut self, barcode: &str, gene: &str, count: u32) {
        self.inner.add_count(barcode, gene, count);
    }

    /// Increment count by 1
    fn increment(&mut self, barcode: &str, gene: &str) {
        self.inner.increment(barcode, gene);
    }

    /// Get number of cells
    fn num_cells(&self) -> usize {
        self.inner.num_cells()
    }

    /// Get number of genes
    fn num_genes(&self) -> usize {
        self.inner.num_genes()
    }

    /// Build the count matrix
    fn build(&mut self) -> PyCountMatrix {
        // We need to take ownership, so create a new counter
        let counter = std::mem::take(&mut self.inner);
        PyCountMatrix {
            inner: counter.build(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "GeneCounter(genes={}, cells={})",
            self.inner.num_genes(),
            self.inner.num_cells()
        )
    }
}
