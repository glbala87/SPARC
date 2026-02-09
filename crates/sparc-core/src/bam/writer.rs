//! BAM file writer using rust-htslib

use crate::{Error, Result};
use rust_htslib::bam::{self, header::HeaderRecord, Header, Write};
use std::path::Path;

/// BAM file writer
pub struct BamWriter {
    writer: bam::Writer,
}

impl BamWriter {
    /// Create a new BAM writer with the given header
    pub fn new<P: AsRef<Path>>(path: P, header: &Header) -> Result<Self> {
        let writer = bam::Writer::from_path(path.as_ref(), header, bam::Format::Bam)
            .map_err(|e| Error::BamParse(format!("Failed to create BAM writer: {}", e)))?;
        Ok(Self { writer })
    }

    /// Create a default header for single-cell data
    pub fn create_default_header() -> Header {
        let mut header = Header::new();

        let mut hd = HeaderRecord::new(b"HD");
        hd.push_tag(b"VN", "1.6");
        hd.push_tag(b"SO", "unsorted");
        header.push_record(&hd);

        let mut pg = HeaderRecord::new(b"PG");
        pg.push_tag(b"ID", "sparc");
        pg.push_tag(b"PN", "sparc");
        pg.push_tag(b"VN", env!("CARGO_PKG_VERSION"));
        header.push_record(&pg);

        header
    }

    /// Write a record
    pub fn write(&mut self, record: &bam::Record) -> Result<()> {
        self.writer
            .write(record)
            .map_err(|e| Error::BamParse(format!("Failed to write record: {}", e)))
    }
}
