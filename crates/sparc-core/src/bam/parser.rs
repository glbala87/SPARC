//! BAM file parser using rust-htslib

use super::BamRecord;
use crate::{Error, Result};
use rust_htslib::bam::{self, Read};
use std::path::Path;

/// BAM file parser
pub struct BamParser {
    reader: bam::Reader,
    header: bam::Header,
}

impl BamParser {
    /// Open a BAM file
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = bam::Reader::from_path(path.as_ref())
            .map_err(|e| Error::BamParse(format!("Failed to open BAM: {}", e)))?;
        let header = bam::Header::from_template(reader.header());
        Ok(Self { reader, header })
    }

    /// Get the header
    pub fn header(&self) -> &bam::Header {
        &self.header
    }

    /// Get reference names
    pub fn reference_names(&self) -> Vec<String> {
        self.reader
            .header()
            .target_names()
            .iter()
            .map(|n| String::from_utf8_lossy(n).to_string())
            .collect()
    }

    /// Convert rust-htslib record to our BamRecord
    fn convert_record(&self, record: &bam::Record) -> BamRecord {
        let name = String::from_utf8_lossy(record.qname()).to_string();
        let seq = record.seq().as_bytes();
        let qual = record.qual().to_vec();

        let mut bam_record = BamRecord::new(name, seq, qual);
        bam_record.mapq = record.mapq();
        bam_record.tid = record.tid();
        bam_record.pos = record.pos();
        bam_record.is_mapped = !record.is_unmapped();
        bam_record.is_reverse = record.is_reverse();

        // Extract CIGAR
        bam_record.cigar = record
            .cigar()
            .iter()
            .map(|c| format!("{}", c))
            .collect::<Vec<_>>()
            .join("");

        // Extract single-cell tags
        if let Ok(aux) = record.aux(b"CB") {
            if let rust_htslib::bam::record::Aux::String(s) = aux {
                bam_record.cell_barcode = Some(s.to_string());
            }
        }
        if let Ok(aux) = record.aux(b"UB") {
            if let rust_htslib::bam::record::Aux::String(s) = aux {
                bam_record.umi = Some(s.to_string());
            }
        }
        if let Ok(aux) = record.aux(b"GN") {
            if let rust_htslib::bam::record::Aux::String(s) = aux {
                bam_record.gene_name = Some(s.to_string());
            }
        }
        if let Ok(aux) = record.aux(b"GX") {
            if let rust_htslib::bam::record::Aux::String(s) = aux {
                bam_record.gene_id = Some(s.to_string());
            }
        }

        bam_record
    }

    /// Read all records
    pub fn read_all(&mut self) -> Result<Vec<BamRecord>> {
        let mut records = Vec::new();
        let mut record = bam::Record::new();

        while self
            .reader
            .read(&mut record)
            .map_err(|e| Error::BamParse(e.to_string()))?
        {
            records.push(self.convert_record(&record));
        }

        Ok(records)
    }

    /// Filter records by mapping quality
    pub fn filter_by_mapq(&mut self, min_mapq: u8) -> Result<Vec<BamRecord>> {
        let mut records = Vec::new();
        let mut record = bam::Record::new();

        while self
            .reader
            .read(&mut record)
            .map_err(|e| Error::BamParse(e.to_string()))?
        {
            if record.mapq() >= min_mapq {
                records.push(self.convert_record(&record));
            }
        }

        Ok(records)
    }
}

impl Iterator for BamParser {
    type Item = Result<BamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut record = bam::Record::new();
        match self.reader.read(&mut record) {
            Ok(true) => Some(Ok(self.convert_record(&record))),
            Ok(false) => None,
            Err(e) => Some(Err(Error::BamParse(e.to_string()))),
        }
    }
}
