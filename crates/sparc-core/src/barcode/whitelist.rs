//! Barcode whitelist handling

use crate::{Error, Result};
use ahash::AHashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Barcode whitelist for exact matching
#[derive(Debug, Clone)]
pub struct Whitelist {
    barcodes: AHashSet<String>,
    barcode_len: usize,
}

impl Whitelist {
    /// Create an empty whitelist
    pub fn new() -> Self {
        Self {
            barcodes: AHashSet::new(),
            barcode_len: 0,
        }
    }

    /// Load whitelist from file (one barcode per line)
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        let reader = BufReader::new(file);

        let mut barcodes = AHashSet::new();
        let mut barcode_len = 0;

        for line in reader.lines() {
            let barcode = line?;
            let barcode = barcode.trim();
            if barcode.is_empty() || barcode.starts_with('#') {
                continue;
            }

            if barcode_len == 0 {
                barcode_len = barcode.len();
            } else if barcode.len() != barcode_len {
                return Err(Error::Barcode(format!(
                    "Inconsistent barcode length: expected {}, got {}",
                    barcode_len,
                    barcode.len()
                )));
            }

            barcodes.insert(barcode.to_string());
        }

        Ok(Self {
            barcodes,
            barcode_len,
        })
    }

    /// Create whitelist from a vector of barcodes
    pub fn from_vec(barcodes: Vec<String>) -> Result<Self> {
        if barcodes.is_empty() {
            return Ok(Self::new());
        }

        let barcode_len = barcodes[0].len();
        for bc in &barcodes {
            if bc.len() != barcode_len {
                return Err(Error::Barcode(format!(
                    "Inconsistent barcode length: expected {}, got {}",
                    barcode_len,
                    bc.len()
                )));
            }
        }

        Ok(Self {
            barcodes: barcodes.into_iter().collect(),
            barcode_len,
        })
    }

    /// Check if a barcode is in the whitelist
    pub fn contains(&self, barcode: &str) -> bool {
        self.barcodes.contains(barcode)
    }

    /// Get the number of barcodes
    pub fn len(&self) -> usize {
        self.barcodes.len()
    }

    /// Check if whitelist is empty
    pub fn is_empty(&self) -> bool {
        self.barcodes.is_empty()
    }

    /// Get expected barcode length
    pub fn barcode_len(&self) -> usize {
        self.barcode_len
    }

    /// Get all barcodes as a vector
    pub fn to_vec(&self) -> Vec<String> {
        self.barcodes.iter().cloned().collect()
    }

    /// Iterate over barcodes
    pub fn iter(&self) -> impl Iterator<Item = &String> {
        self.barcodes.iter()
    }
}

impl Default for Whitelist {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_whitelist_from_vec() {
        let barcodes = vec![
            "AAACCCAAGAAACACT".to_string(),
            "AAACCCAAGAAACCAT".to_string(),
            "AAACCCAAGAAACCCA".to_string(),
        ];

        let whitelist = Whitelist::from_vec(barcodes).unwrap();
        assert_eq!(whitelist.len(), 3);
        assert_eq!(whitelist.barcode_len(), 16);
        assert!(whitelist.contains("AAACCCAAGAAACACT"));
        assert!(!whitelist.contains("AAACCCAAGAAACXXX"));
    }
}
