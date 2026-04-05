//! 10x Genomics 3' Gene Expression kit implementation

use super::{Protocol, ReadComponents};
use crate::{Error, ReadStructure, Result};

/// 10x Genomics 3' v3 protocol
///
/// Read structure:
/// - R1: Barcode (16bp) + UMI (12bp)
/// - R2: cDNA
pub struct TenX3Prime {
    read_structure: ReadStructure,
    version: String,
}

impl TenX3Prime {
    /// Create a new 10x 3' v3 protocol
    pub fn v3() -> Self {
        Self {
            read_structure: ReadStructure::new(0, 16, 16, 12, 0),
            version: "v3".to_string(),
        }
    }

    /// Create a new 10x 3' v2 protocol (shorter UMI)
    pub fn v2() -> Self {
        Self {
            read_structure: ReadStructure::new(0, 16, 16, 10, 0),
            version: "v2".to_string(),
        }
    }

    /// Create with custom read structure
    pub fn custom(read_structure: ReadStructure) -> Self {
        Self {
            read_structure,
            version: "custom".to_string(),
        }
    }
}

impl Protocol for TenX3Prime {
    fn read_structure(&self) -> &ReadStructure {
        &self.read_structure
    }

    fn extract_r1(&self, seq: &[u8], qual: &[u8]) -> Result<ReadComponents> {
        let rs = &self.read_structure;
        let min_len = rs.barcode_start + rs.barcode_len + rs.umi_len;

        if seq.len() < min_len {
            return Err(Error::Protocol(format!(
                "R1 too short: {} < {} required",
                seq.len(),
                min_len
            )));
        }

        let barcode_end = rs.barcode_start + rs.barcode_len;
        let umi_end = rs.umi_start + rs.umi_len;

        Ok(ReadComponents {
            barcode: seq[rs.barcode_start..barcode_end].to_vec(),
            umi: seq[rs.umi_start..umi_end].to_vec(),
            cdna: Vec::new(), // cDNA is on R2
            barcode_qual: qual[rs.barcode_start..barcode_end].to_vec(),
            umi_qual: qual[rs.umi_start..umi_end].to_vec(),
            cdna_qual: Vec::new(),
        })
    }

    fn name(&self) -> &str {
        "10x Genomics 3' Gene Expression"
    }

    fn version(&self) -> &str {
        &self.version
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_v3() {
        let protocol = TenX3Prime::v3();

        // 16bp barcode + 12bp UMI
        let seq = b"AAACCCAAGAAACACTGGGGTTTTAAAA";
        let qual = b"IIIIIIIIIIIIIIIIIIIIIIIIIIII";

        let components = protocol.extract_r1(seq, qual).unwrap();

        assert_eq!(components.barcode_str(), "AAACCCAAGAAACACT");
        assert_eq!(components.umi_str(), "GGGGTTTTAAAA");
        assert_eq!(components.barcode.len(), 16);
        assert_eq!(components.umi.len(), 12);
    }

    #[test]
    fn test_extract_too_short() {
        let protocol = TenX3Prime::v3();

        let seq = b"AAACCCAAGAAACACT"; // Only 16bp, need 28
        let qual = b"IIIIIIIIIIIIIIII";

        let result = protocol.extract_r1(seq, qual);
        assert!(result.is_err());
    }
}
