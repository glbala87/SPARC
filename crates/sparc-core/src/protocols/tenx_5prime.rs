//! 10x Genomics 5' Gene Expression kit implementation

use super::{Protocol, ReadComponents};
use crate::{Error, ReadStructure, Result};

/// 10x Genomics 5' v2 protocol
///
/// Read structure:
/// - R1: Barcode (16bp) + UMI (10bp)
/// - R2: cDNA (5' end)
pub struct TenX5Prime {
    read_structure: ReadStructure,
    version: String,
}

impl TenX5Prime {
    /// Create a new 10x 5' v2 protocol
    pub fn v2() -> Self {
        Self {
            read_structure: ReadStructure::new(0, 16, 16, 10, 0),
            version: "v2".to_string(),
        }
    }

    /// Create a new 10x 5' v1.1 protocol
    pub fn v1_1() -> Self {
        Self {
            read_structure: ReadStructure::new(0, 16, 16, 10, 0),
            version: "v1.1".to_string(),
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

impl Protocol for TenX5Prime {
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
        "10x Genomics 5' Gene Expression"
    }

    fn version(&self) -> &str {
        &self.version
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_v2() {
        let protocol = TenX5Prime::v2();

        // 16bp barcode + 10bp UMI
        let seq = b"AAACCCAAGAAACACTGGGGTTTTAA";
        let qual = b"IIIIIIIIIIIIIIIIIIIIIIIIII";

        let components = protocol.extract_r1(seq, qual).unwrap();

        assert_eq!(components.barcode_str(), "AAACCCAAGAAACACT");
        assert_eq!(components.umi_str(), "GGGGTTTTAA");
        assert_eq!(components.barcode.len(), 16);
        assert_eq!(components.umi.len(), 10);
    }
}
