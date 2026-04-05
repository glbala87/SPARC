//! Python bindings for QC metrics

use pyo3::prelude::*;
use sparc_core::qc::{QcMetrics, QcReport};

/// Python wrapper for QcMetrics
#[pyclass(name = "QcMetrics")]
pub struct PyQcMetrics {
    inner: QcMetrics,
}

#[pymethods]
impl PyQcMetrics {
    #[new]
    fn new() -> Self {
        Self { inner: QcMetrics::new() }
    }

    #[getter]
    fn total_reads(&self) -> u64 { self.inner.total_reads }

    #[getter]
    fn valid_barcode_reads(&self) -> u64 { self.inner.valid_barcode_reads }

    #[getter]
    fn mapped_reads(&self) -> u64 { self.inner.mapped_reads }

    #[getter]
    fn assigned_reads(&self) -> u64 { self.inner.assigned_reads }

    #[getter]
    fn num_cells(&self) -> u64 { self.inner.num_cells }

    #[getter]
    fn total_genes(&self) -> u64 { self.inner.total_genes }

    #[getter]
    fn mean_reads_per_cell(&self) -> f64 { self.inner.mean_reads_per_cell }

    #[getter]
    fn median_reads_per_cell(&self) -> f64 { self.inner.median_reads_per_cell }

    #[getter]
    fn mean_genes_per_cell(&self) -> f64 { self.inner.mean_genes_per_cell }

    #[getter]
    fn median_genes_per_cell(&self) -> f64 { self.inner.median_genes_per_cell }

    #[getter]
    fn median_umi_per_cell(&self) -> f64 { self.inner.median_umi_per_cell }

    #[getter]
    fn sequencing_saturation(&self) -> f64 { self.inner.sequencing_saturation }

    fn barcode_validity_rate(&self) -> f64 { self.inner.barcode_validity_rate() }

    fn mapping_rate(&self) -> f64 { self.inner.mapping_rate() }

    fn assignment_rate(&self) -> f64 { self.inner.assignment_rate() }

    fn __repr__(&self) -> String {
        format!(
            "QcMetrics(cells={}, genes={}, reads={})",
            self.inner.num_cells, self.inner.total_genes, self.inner.total_reads
        )
    }
}

/// Python wrapper for QcReport
#[pyclass(name = "QcReport")]
pub struct PyQcReport {
    inner: QcReport,
}

#[pymethods]
impl PyQcReport {
    #[new]
    fn new(sample_name: String) -> Self {
        Self { inner: QcReport::new(sample_name) }
    }

    #[getter]
    fn sample_name(&self) -> &str { &self.inner.sample_name }

    #[getter]
    fn warnings(&self) -> Vec<String> { self.inner.warnings.clone() }

    fn generate_warnings(&mut self) {
        self.inner.generate_warnings();
    }

    fn to_json(&self) -> PyResult<String> {
        self.inner.to_json()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!(
            "QcReport(sample='{}', warnings={})",
            self.inner.sample_name, self.inner.warnings.len()
        )
    }
}
