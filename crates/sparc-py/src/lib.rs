//! Python bindings for SPARC

mod analysis;
mod bam;
mod barcode;
mod fastq;
mod matrix;
mod qc;
mod validation_py;

use pyo3::prelude::*;

/// SPARC Python module
#[pymodule]
fn sparc_py(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // Core I/O classes
    m.add_class::<fastq::PyFastqParser>()?;
    m.add_class::<fastq::PyFastqRecord>()?;
    m.add_class::<fastq::PyFastqWriter>()?;
    m.add_class::<bam::PyBamParser>()?;
    m.add_class::<bam::PyBamRecord>()?;

    // Barcode classes
    m.add_class::<barcode::PyWhitelist>()?;
    m.add_class::<barcode::PyBarcodeCorrector>()?;

    // Count matrix classes
    m.add_class::<matrix::PyCountMatrix>()?;
    m.add_class::<matrix::PyGeneCounter>()?;

    // QC classes
    m.add_class::<qc::PyQcMetrics>()?;
    m.add_class::<qc::PyQcReport>()?;

    // Analysis functions
    m.add_function(wrap_pyfunction!(analysis::py_normalize_total, m)?)?;
    m.add_function(wrap_pyfunction!(analysis::py_scale, m)?)?;
    m.add_function(wrap_pyfunction!(analysis::py_highly_variable_genes, m)?)?;
    m.add_function(wrap_pyfunction!(analysis::py_pca, m)?)?;
    m.add_function(wrap_pyfunction!(analysis::py_build_knn_graph, m)?)?;
    m.add_function(wrap_pyfunction!(analysis::py_label_propagation, m)?)?;
    m.add_function(wrap_pyfunction!(analysis::py_run_analysis, m)?)?;

    // Validation functions
    m.add_function(wrap_pyfunction!(validation_py::adjusted_rand_index, m)?)?;
    m.add_function(wrap_pyfunction!(validation_py::normalized_mutual_info, m)?)?;
    m.add_function(wrap_pyfunction!(validation_py::pearson_correlation, m)?)?;
    m.add_function(wrap_pyfunction!(validation_py::spearman_correlation, m)?)?;
    m.add_function(wrap_pyfunction!(validation_py::f1_score, m)?)?;

    // Module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
