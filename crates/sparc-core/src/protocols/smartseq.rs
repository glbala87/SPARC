//! SMART-seq2 protocol implementation

use super::{Protocol, ReadComponents};
use crate::{ReadStructure, Result};

/// SMART-seq2 protocol (plate-based, no barcode/UMI per read)
///
/// Each file represents one cell. The barcode is the sample/well name.
/// Reads are full-length cDNA with no embedded barcode or UMI.
pub struct SmartSeq2 {
    read_structure: ReadStructure,
    sample_name: String,
}

impl SmartSeq2 {
    pub fn new(sample_name: String) -> Self {
        Self {
            read_structure: ReadStructure::new(0, 0, 0, 0, 0),
            sample_name,
        }
    }

    pub fn with_name(name: &str) -> Self {
        Self::new(name.to_string())
    }
}

impl Protocol for SmartSeq2 {
    fn read_structure(&self) -> &ReadStructure {
        &self.read_structure
    }

    fn extract_r1(&self, seq: &[u8], qual: &[u8]) -> Result<ReadComponents> {
        // No barcode or UMI - entire read is cDNA
        // Use sample name as barcode
        Ok(ReadComponents {
            barcode: self.sample_name.as_bytes().to_vec(),
            umi: Vec::new(),
            cdna: seq.to_vec(),
            barcode_qual: vec![b'I'; self.sample_name.len()],
            umi_qual: Vec::new(),
            cdna_qual: qual.to_vec(),
        })
    }

    fn name(&self) -> &str {
        "SMART-seq2"
    }

    fn version(&self) -> &str {
        "v2"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smartseq2_extraction() {
        let protocol = SmartSeq2::new("WellA01".to_string());
        let seq = b"ACGTACGTACGTACGTACGTACGT";
        let qual = b"IIIIIIIIIIIIIIIIIIIIIIII";

        let components = protocol.extract_r1(seq, qual).unwrap();
        assert_eq!(components.barcode_str(), "WellA01");
        assert!(components.umi.is_empty());
        assert_eq!(components.cdna.len(), 24);
    }

    #[test]
    fn test_smartseq2_empty_read() {
        let protocol = SmartSeq2::new("Well".to_string());
        let seq = b"";
        let qual = b"";

        let components = protocol.extract_r1(seq, qual).unwrap();
        assert_eq!(components.barcode_str(), "Well");
        assert!(components.cdna.is_empty());
    }
}
