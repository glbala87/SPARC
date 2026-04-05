//! # SPARC Core
//!
//! SPARC: Single-cell Pipeline Accelerated in Rust Core
//!
//! This crate provides the core functionality for processing single-cell sequencing data,
//! including FASTQ/BAM parsing, barcode detection, UMI deduplication, count matrix generation,
//! alignment integration, streaming processing, and downstream analysis.

pub mod aligner;
pub mod analysis;
pub mod bam;
pub mod barcode;
pub mod count;
pub mod fastq;
pub mod protocols;
pub mod qc;
pub mod streaming;
pub mod umi;
pub mod validation;

pub use aligner::{Aligner, AlignerConfig, AlignerType};
pub use bam::{BamParser, BamRecord, BamWriter};
pub use barcode::{BarcodeCorrector, BarcodeMatcher, Whitelist};
pub use count::{CountMatrix, CsrMatrix, GeneCounter};
pub use fastq::{FastqParser, FastqRecord, FastqWriter};
pub use protocols::{DropSeq, InDrop, Protocol, SciRNA, SmartSeq2, TenX3Prime, TenX5Prime};
pub use qc::{QcMetrics, QcReport};
pub use streaming::{StreamConfig, StreamStats, StreamingProcessor};
pub use umi::{UmiDeduplicator, UmiGraph};
pub use validation::{ValidationReport, SyntheticConfig, SyntheticDataset, TruthSet};

/// Error types for SPARC core
#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("FASTQ parsing error: {0}")]
    FastqParse(String),

    #[error("BAM parsing error: {0}")]
    BamParse(String),

    #[error("Barcode error: {0}")]
    Barcode(String),

    #[error("UMI error: {0}")]
    Umi(String),

    #[error("Protocol error: {0}")]
    Protocol(String),

    #[error("Invalid read structure: {0}")]
    ReadStructure(String),
}

pub type Result<T> = std::result::Result<T, Error>;

/// Read structure definition for parsing sequencing reads
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ReadStructure {
    /// Barcode start position (0-indexed)
    pub barcode_start: usize,
    /// Barcode length
    pub barcode_len: usize,
    /// UMI start position (0-indexed)
    pub umi_start: usize,
    /// UMI length
    pub umi_len: usize,
    /// cDNA start position (0-indexed)
    pub cdna_start: usize,
}

impl ReadStructure {
    pub fn new(
        barcode_start: usize,
        barcode_len: usize,
        umi_start: usize,
        umi_len: usize,
        cdna_start: usize,
    ) -> Self {
        Self {
            barcode_start,
            barcode_len,
            umi_start,
            umi_len,
            cdna_start,
        }
    }

    /// 10x Genomics 3' v3 read structure
    pub fn tenx_3prime_v3() -> Self {
        Self::new(0, 16, 16, 12, 0)
    }

    /// 10x Genomics 5' v2 read structure
    pub fn tenx_5prime_v2() -> Self {
        Self::new(0, 16, 16, 10, 0)
    }

    /// Drop-seq read structure
    pub fn dropseq() -> Self {
        Self::new(0, 12, 12, 8, 0)
    }

    /// inDrop read structure
    pub fn indrop() -> Self {
        Self::new(0, 16, 16, 6, 0)
    }

    /// sci-RNA-seq read structure
    pub fn scirna() -> Self {
        Self::new(0, 10, 10, 8, 0)
    }
}
