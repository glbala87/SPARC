//! 10x Genomics protocol implementations

mod tenx_3prime;
mod tenx_5prime;

pub use tenx_3prime::TenX3Prime;
pub use tenx_5prime::TenX5Prime;

use crate::{ReadStructure, Result};

/// Extracted read components
#[derive(Debug, Clone)]
pub struct ReadComponents {
    /// Cell barcode sequence
    pub barcode: Vec<u8>,
    /// UMI sequence
    pub umi: Vec<u8>,
    /// cDNA sequence
    pub cdna: Vec<u8>,
    /// Barcode quality scores
    pub barcode_qual: Vec<u8>,
    /// UMI quality scores
    pub umi_qual: Vec<u8>,
    /// cDNA quality scores
    pub cdna_qual: Vec<u8>,
}

impl ReadComponents {
    /// Get barcode as string
    pub fn barcode_str(&self) -> String {
        String::from_utf8_lossy(&self.barcode).to_string()
    }

    /// Get UMI as string
    pub fn umi_str(&self) -> String {
        String::from_utf8_lossy(&self.umi).to_string()
    }

    /// Check if barcode has good quality (mean Q >= threshold)
    pub fn barcode_quality_ok(&self, min_qual: u8) -> bool {
        if self.barcode_qual.is_empty() {
            return false;
        }
        let mean_qual: f64 = self.barcode_qual.iter().map(|&q| (q - 33) as f64).sum::<f64>()
            / self.barcode_qual.len() as f64;
        mean_qual >= min_qual as f64
    }

    /// Check if UMI has good quality (mean Q >= threshold)
    pub fn umi_quality_ok(&self, min_qual: u8) -> bool {
        if self.umi_qual.is_empty() {
            return false;
        }
        let mean_qual: f64 =
            self.umi_qual.iter().map(|&q| (q - 33) as f64).sum::<f64>() / self.umi_qual.len() as f64;
        mean_qual >= min_qual as f64
    }
}

/// Protocol trait for different single-cell sequencing kits
pub trait Protocol: Send + Sync {
    /// Get the read structure for this protocol
    fn read_structure(&self) -> &ReadStructure;

    /// Extract components from R1 read
    fn extract_r1(&self, seq: &[u8], qual: &[u8]) -> Result<ReadComponents>;

    /// Protocol name
    fn name(&self) -> &str;

    /// Protocol version
    fn version(&self) -> &str;
}
