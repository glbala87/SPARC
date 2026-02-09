//! Python bindings for SPARC

mod bam;
mod barcode;
mod fastq;
mod matrix;

use pyo3::prelude::*;

/// SPARC Python module
#[pymodule]
fn sparc_py(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<fastq::PyFastqParser>()?;
    m.add_class::<fastq::PyFastqRecord>()?;
    m.add_class::<fastq::PyFastqWriter>()?;
    m.add_class::<bam::PyBamParser>()?;
    m.add_class::<bam::PyBamRecord>()?;
    m.add_class::<barcode::PyWhitelist>()?;
    m.add_class::<barcode::PyBarcodeCorrector>()?;
    m.add_class::<matrix::PyCountMatrix>()?;
    m.add_class::<matrix::PyGeneCounter>()?;

    // Module metadata
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
