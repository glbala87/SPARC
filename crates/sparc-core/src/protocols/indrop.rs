//! inDrop protocol implementation

use super::{Protocol, ReadComponents};
use crate::{Error, ReadStructure, Result};

/// inDrop protocol
///
/// Read structure:
/// - R1: Barcode Part 1 (8bp) + Barcode Part 2 (8bp) + UMI (6bp)
/// - R2: cDNA
/// - Combined barcode = part1 + part2 (16bp total)
pub struct InDrop {
    read_structure: ReadStructure,
}

impl InDrop {
    pub fn new() -> Self {
        Self {
            // Total barcode = 16bp (two 8bp halves), UMI = 6bp
            read_structure: ReadStructure::new(0, 16, 16, 6, 0),
        }
    }

    pub fn custom(read_structure: ReadStructure) -> Self {
        Self { read_structure }
    }
}

impl Default for InDrop {
    fn default() -> Self {
        Self::new()
    }
}

impl Protocol for InDrop {
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
        "inDrop"
    }

    fn version(&self) -> &str {
        "v3"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_indrop_extraction() {
        let protocol = InDrop::new();
        // 16bp barcode (8+8) + 6bp UMI = 22bp minimum
        let seq = b"AAAAGGGGCCCCTTTTAAAAGG";
        let qual = b"IIIIIIIIIIIIIIIIIIIIII";

        let components = protocol.extract_r1(seq, qual).unwrap();
        assert_eq!(components.barcode.len(), 16);
        assert_eq!(components.umi.len(), 6);
        assert_eq!(components.barcode_str(), "AAAAGGGGCCCCTTTT");
        assert_eq!(components.umi_str(), "AAAAGG");
    }

    #[test]
    fn test_indrop_too_short() {
        let protocol = InDrop::new();
        let seq = b"AAAAGGGGCCCC"; // Only 12bp
        let qual = b"IIIIIIIIIIII";
        assert!(protocol.extract_r1(seq, qual).is_err());
    }
}
