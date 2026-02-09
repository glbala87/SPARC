//! FASTQ Python bindings

use pyo3::prelude::*;
use sparc_core::fastq::{FastqParser, FastqRecord, FastqWriter};

/// Python wrapper for FastqRecord
#[pyclass(name = "FastqRecord")]
pub struct PyFastqRecord {
    inner: FastqRecord,
}

#[pymethods]
impl PyFastqRecord {
    #[new]
    fn new(id: String, seq: Vec<u8>, qual: Vec<u8>) -> Self {
        Self {
            inner: FastqRecord::new(id, seq, qual),
        }
    }

    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }

    #[getter]
    fn seq(&self) -> &[u8] {
        &self.inner.seq
    }

    #[getter]
    fn qual(&self) -> &[u8] {
        &self.inner.qual
    }

    /// Get sequence as string
    fn seq_str(&self) -> String {
        String::from_utf8_lossy(&self.inner.seq).to_string()
    }

    /// Get quality as string
    fn qual_str(&self) -> String {
        String::from_utf8_lossy(&self.inner.qual).to_string()
    }

    /// Get mean quality score
    fn mean_quality(&self) -> f64 {
        self.inner.mean_quality()
    }

    /// Get subsequence
    fn subsequence(&self, start: usize, len: usize) -> Option<Vec<u8>> {
        self.inner.subsequence(start, len).map(|s| s.to_vec())
    }

    fn __repr__(&self) -> String {
        format!(
            "FastqRecord(id='{}', len={})",
            self.inner.id,
            self.inner.seq.len()
        )
    }
}

/// Python wrapper for FastqParser
#[pyclass(name = "FastqParser")]
pub struct PyFastqParser {
    inner: FastqParser,
}

#[pymethods]
impl PyFastqParser {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let inner = FastqParser::open(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(Self { inner })
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> PyResult<Option<PyFastqRecord>> {
        match self.inner.next() {
            Some(Ok(record)) => Ok(Some(PyFastqRecord { inner: record })),
            Some(Err(e)) => Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string())),
            None => Ok(None),
        }
    }

    /// Read all records into a list
    fn read_all(&mut self) -> PyResult<Vec<PyFastqRecord>> {
        let mut records = Vec::new();
        while let Some(result) = self.inner.next() {
            let record = result
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
            records.push(PyFastqRecord { inner: record });
        }
        Ok(records)
    }
}

/// Python wrapper for FastqWriter
#[pyclass(name = "FastqWriter")]
pub struct PyFastqWriter {
    inner: Option<FastqWriter>,
}

#[pymethods]
impl PyFastqWriter {
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let inner = FastqWriter::new(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(Self { inner: Some(inner) })
    }

    /// Write a record
    fn write(&mut self, record: &PyFastqRecord) -> PyResult<()> {
        if let Some(ref mut writer) = self.inner {
            writer
                .write_record(&record.inner)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                "Writer is closed",
            ))
        }
    }

    /// Flush and close the writer
    fn close(&mut self) -> PyResult<()> {
        if let Some(mut writer) = self.inner.take() {
            writer
                .flush()
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))
        } else {
            Ok(())
        }
    }

    fn __enter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __exit__(
        &mut self,
        _exc_type: Option<PyObject>,
        _exc_value: Option<PyObject>,
        _traceback: Option<PyObject>,
    ) -> PyResult<bool> {
        self.close()?;
        Ok(false)
    }
}
