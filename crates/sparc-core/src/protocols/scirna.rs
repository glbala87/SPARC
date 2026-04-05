//! sci-RNA-seq protocol implementation

use super::{Protocol, ReadComponents};
use crate::{Error, ReadStructure, Result};

/// sci-RNA-seq protocol
///
/// Read structure:
/// - R1: Ligation barcode (10bp) + UMI (8bp)
/// - R2: cDNA
pub struct SciRNA {
    read_structure: ReadStructure,
}

impl SciRNA {
    pub fn new() -> Self {
        Self {
            read_structure: ReadStructure::new(0, 10, 10, 8, 0),
        }
    }

    pub fn custom(read_structure: ReadStructure) -> Self {
        Self { read_structure }
    }
}

impl Default for SciRNA {
    fn default() -> Self {
        Self::new()
    }
}

impl Protocol for SciRNA {
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
            cdna: Vec::new(),
            barcode_qual: qual[rs.barcode_start..barcode_end].to_vec(),
            umi_qual: qual[rs.umi_start..umi_end].to_vec(),
            cdna_qual: Vec::new(),
        })
    }

    fn name(&self) -> &str {
        "sci-RNA-seq"
    }

    fn version(&self) -> &str {
        "v1"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scirna_extraction() {
        let protocol = SciRNA::new();
        // 10bp barcode + 8bp UMI = 18bp minimum
        let seq = b"AAAAGGGGCCTTTTAAAA";
        let qual = b"IIIIIIIIIIIIIIIIII";

        let components = protocol.extract_r1(seq, qual).unwrap();
        assert_eq!(components.barcode.len(), 10);
        assert_eq!(components.umi.len(), 8);
        assert_eq!(components.barcode_str(), "AAAAGGGGCC");
        assert_eq!(components.umi_str(), "TTTTAAAA");
    }
}
