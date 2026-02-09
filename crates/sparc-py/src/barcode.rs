//! Barcode Python bindings

use pyo3::prelude::*;
use sparc_core::barcode::{BarcodeCorrector, BarcodeMatch, Whitelist};

/// Python wrapper for Whitelist
#[pyclass(name = "Whitelist")]
pub struct PyWhitelist {
    inner: Whitelist,
}

#[pymethods]
impl PyWhitelist {
    /// Create whitelist from file
    #[new]
    fn new(path: &str) -> PyResult<Self> {
        let inner = Whitelist::from_file(path)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Create whitelist from list of barcodes
    #[staticmethod]
    fn from_list(barcodes: Vec<String>) -> PyResult<Self> {
        let inner = Whitelist::from_vec(barcodes)
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Check if barcode is in whitelist
    fn contains(&self, barcode: &str) -> bool {
        self.inner.contains(barcode)
    }

    /// Get number of barcodes
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Get expected barcode length
    fn barcode_len(&self) -> usize {
        self.inner.barcode_len()
    }

    /// Get all barcodes as list
    fn to_list(&self) -> Vec<String> {
        self.inner.to_vec()
    }

    fn __repr__(&self) -> String {
        format!("Whitelist(n={}, len={})", self.inner.len(), self.inner.barcode_len())
    }
}

/// Python wrapper for BarcodeCorrector
#[pyclass(name = "BarcodeCorrector")]
pub struct PyBarcodeCorrector {
    inner: BarcodeCorrector,
}

#[pymethods]
impl PyBarcodeCorrector {
    /// Create a barcode corrector
    #[new]
    fn new(whitelist: &PyWhitelist, max_distance: u32) -> Self {
        let inner = BarcodeCorrector::new(whitelist.inner.clone(), max_distance);
        Self { inner }
    }

    /// Match a barcode, returning (status, corrected_barcode, distance)
    /// status: "exact", "corrected", or "no_match"
    fn match_barcode(&self, barcode: &str) -> (String, Option<String>, u32) {
        match self.inner.match_barcode(barcode) {
            BarcodeMatch::Exact(bc) => ("exact".to_string(), Some(bc), 0),
            BarcodeMatch::Corrected(_, bc, dist) => ("corrected".to_string(), Some(bc), dist),
            BarcodeMatch::NoMatch(_) => ("no_match".to_string(), None, 0),
        }
    }

    /// Check if barcode is valid (exact or correctable)
    fn is_valid(&self, barcode: &str) -> bool {
        self.inner.match_barcode(barcode).is_valid()
    }

    /// Get corrected barcode or None
    fn correct(&self, barcode: &str) -> Option<String> {
        self.inner.match_barcode(barcode).barcode().map(|s| s.to_string())
    }

    /// Batch correct barcodes
    fn correct_batch(&self, barcodes: Vec<String>) -> Vec<Option<String>> {
        barcodes
            .iter()
            .map(|bc| self.inner.match_barcode(bc).barcode().map(|s| s.to_string()))
            .collect()
    }
}
