//! FASTQ file writer with compression support

use super::FastqRecord;
use crate::{Error, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// FASTQ writer supporting plain text and gzip compression
pub struct FastqWriter {
    writer: Box<dyn Write>,
}

impl FastqWriter {
    /// Create a new FASTQ writer
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let file = File::create(path)?;

        let writer: Box<dyn Write> = if path
            .extension()
            .map_or(false, |ext| ext == "gz" || ext == "gzip")
        {
            Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
        } else {
            Box::new(BufWriter::new(file))
        };

        Ok(Self { writer })
    }

    /// Write a FASTQ record
    pub fn write_record(&mut self, record: &FastqRecord) -> Result<()> {
        writeln!(self.writer, "@{}", record.id)?;
        self.writer.write_all(&record.seq)?;
        writeln!(self.writer)?;
        writeln!(self.writer, "+")?;
        self.writer.write_all(&record.qual)?;
        writeln!(self.writer)?;
        Ok(())
    }

    /// Write multiple records
    pub fn write_records(&mut self, records: &[FastqRecord]) -> Result<()> {
        for record in records {
            self.write_record(record)?;
        }
        Ok(())
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().map_err(Error::from)
    }
}

impl Drop for FastqWriter {
    fn drop(&mut self) {
        let _ = self.flush();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_write_fastq() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.fastq");

        let record = FastqRecord::new(
            "read1".to_string(),
            b"ACGTACGT".to_vec(),
            b"IIIIIIII".to_vec(),
        );

        let mut writer = FastqWriter::new(&path).unwrap();
        writer.write_record(&record).unwrap();
        writer.flush().unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("@read1"));
        assert!(content.contains("ACGTACGT"));
    }
}
