//! BAM Python bindings

use pyo3::prelude::*;
use sparc_core::bam::{BamParser, BamRecord};

/// Python wrapper for BamRecord
#[pyclass(name = "BamRecord")]
pub struct PyBamRecord {
    inner: BamRecord,
}

#[pymethods]
impl PyBamRecord {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn seq(&self) -> &[u8] {
        &self.inner.seq
    }

    #[getter]
    fn qual(&self) -> &[u8] {
        &self.inner.qual
    }

    #[getter]
    fn mapq(&self) -> u8 {
        self.inner.mapq
    }

    #[getter]
    fn tid(&self) -> i32 {
        self.inner.tid
    }

    #[getter]
    fn pos(&self) -> i64 {
        self.inner.pos
    }

    #[getter]
    fn cigar(&self) -> &str {
        &self.inner.cigar
    }

    #[getter]
    fn cell_barcode(&self) -> Option<&str> {
        self.inner.cell_barcode.as_deref()
    }

    #[getter]
    fn umi(&self) -> Option<&str> {
        self.inner.umi.as_deref()
    }

    #[getter]
    fn gene_name(&self) -> Option<&str> {
        self.inner.gene_name.as_deref()
    }

    #[getter]
    fn gene_id(&self) -> Option<&str> {
        self.inner.gene_id.as_deref()
    }

    #[getter]
    fn is_mapped(&self) -> bool {
        self.inner.is_mapped
    }

    #[getter]
    fn is_reverse(&self) -> bool {
        self.inner.is_reverse
    }

    /// Check if record has valid cell barcode and UMI
    fn has_valid_tags(&self) -> bool {
        self.inner.has_valid_tags()
    }

    /// Check if record is assigned to a gene
    fn is_assigned(&self) -> bool {
        self.inner.is_assigned()
    }

    fn __repr__(&self) -> String {
        format!(
            "BamRecord(name='{}', pos={}, mapq={})",
            self.inner.name, self.inner.pos, self.inner.mapq
        )
    }
}

/// Python wrapper for BamParser
#[pyclass(name = "BamParser")]
pub struct PyBamParser {
    inner: BamParser,
}

#[pymethods]
impl PyBamParser {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let inner = BamParser::open(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Get reference names from header
    fn reference_names(&self) -> Vec<String> {
        self.inner.reference_names()
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<PyBamRecord>> {
        match self.inner.next() {
            Some(Ok(record)) => Ok(Some(PyBamRecord { inner: record })),
            Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string())),
            None => Ok(None),
        }
    }

    /// Read all records into a list
    fn read_all(&mut self) -> PyResult<Vec<PyBamRecord>> {
        let records = self
            .inner
            .read_all()
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(records.into_iter().map(|r| PyBamRecord { inner: r }).collect())
    }

    /// Filter records by mapping quality
    fn filter_by_mapq(&mut self, min_mapq: u8) -> PyResult<Vec<PyBamRecord>> {
        let records = self
            .inner
            .filter_by_mapq(min_mapq)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(records.into_iter().map(|r| PyBamRecord { inner: r }).collect())
    }
}
