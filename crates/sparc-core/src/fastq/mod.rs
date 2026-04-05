//! FASTQ parsing and writing module

mod parser;
mod writer;

pub use parser::FastqParser;
pub use writer::FastqWriter;

/// A FASTQ record
#[derive(Debug, Clone)]
pub struct FastqRecord {
    /// Read identifier
    pub id: String,
    /// Sequence data
    pub seq: Vec<u8>,
    /// Quality scores (Phred+33 encoded)
    pub qual: Vec<u8>,
}

impl FastqRecord {
    pub fn new(id: String, seq: Vec<u8>, qual: Vec<u8>) -> Self {
        Self { id, seq, qual }
    }

    /// Extract a subsequence from the record
    pub fn subsequence(&self, start: usize, len: usize) -> Option<&[u8]> {
        if start + len <= self.seq.len() {
            Some(&self.seq[start..start + len])
        } else {
            None
        }
    }

    /// Extract quality scores for a subsequence
    pub fn subqual(&self, start: usize, len: usize) -> Option<&[u8]> {
        if start + len <= self.qual.len() {
            Some(&self.qual[start..start + len])
        } else {
            None
        }
    }

    /// Calculate mean quality score for the entire read
    pub fn mean_quality(&self) -> f64 {
        if self.qual.is_empty() {
            return 0.0;
        }
        let sum: u64 = self.qual.iter().map(|&q| (q - 33) as u64).sum();
        sum as f64 / self.qual.len() as f64
    }

    /// Calculate mean quality score for a region
    pub fn mean_quality_region(&self, start: usize, len: usize) -> Option<f64> {
        let region = self.subqual(start, len)?;
        if region.is_empty() {
            return Some(0.0);
        }
        let sum: u64 = region.iter().map(|&q| (q - 33) as u64).sum();
        Some(sum as f64 / region.len() as f64)
    }
}
