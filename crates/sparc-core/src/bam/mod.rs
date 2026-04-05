//! BAM parsing and writing module

mod parser;
mod writer;

pub use parser::BamParser;
pub use writer::BamWriter;

/// A BAM record with extracted single-cell tags
#[derive(Debug, Clone)]
pub struct BamRecord {
    /// Read name
    pub name: String,
    /// Sequence
    pub seq: Vec<u8>,
    /// Quality scores
    pub qual: Vec<u8>,
    /// Mapping quality
    pub mapq: u8,
    /// Reference ID (-1 for unmapped)
    pub tid: i32,
    /// Position (0-based)
    pub pos: i64,
    /// CIGAR string
    pub cigar: String,
    /// Cell barcode (CB tag)
    pub cell_barcode: Option<String>,
    /// UMI (UB tag)
    pub umi: Option<String>,
    /// Gene name (GN tag)
    pub gene_name: Option<String>,
    /// Gene ID (GX tag)
    pub gene_id: Option<String>,
    /// Is mapped
    pub is_mapped: bool,
    /// Is reverse strand
    pub is_reverse: bool,
}

impl BamRecord {
    pub fn new(name: String, seq: Vec<u8>, qual: Vec<u8>) -> Self {
        Self {
            name,
            seq,
            qual,
            mapq: 0,
            tid: -1,
            pos: -1,
            cigar: String::new(),
            cell_barcode: None,
            umi: None,
            gene_name: None,
            gene_id: None,
            is_mapped: false,
            is_reverse: false,
        }
    }

    /// Check if this record has valid cell barcode and UMI
    pub fn has_valid_tags(&self) -> bool {
        self.cell_barcode.is_some() && self.umi.is_some()
    }

    /// Check if this record is assigned to a gene
    pub fn is_assigned(&self) -> bool {
        self.gene_name.is_some() || self.gene_id.is_some()
    }
}
