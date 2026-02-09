//! FASTQ file parser with parallel processing support

use super::FastqRecord;
use crate::{Error, Result};
use needletail::{parse_fastx_file, FastxReader};
use rayon::prelude::*;
use std::path::Path;

/// Parallel FASTQ parser using needletail
pub struct FastqParser {
    reader: Box<dyn FastxReader>,
}

impl FastqParser {
    /// Open a FASTQ file (supports .gz and .zst compression)
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = parse_fastx_file(path.as_ref())
            .map_err(|e| Error::FastqParse(format!("Failed to open FASTQ: {}", e)))?;
        Ok(Self { reader })
    }

    /// Read all records into memory
    pub fn read_all(&mut self) -> Result<Vec<FastqRecord>> {
        let mut records = Vec::new();
        while let Some(record) = self.reader.next() {
            let record =
                record.map_err(|e| Error::FastqParse(format!("Failed to read record: {}", e)))?;
            records.push(FastqRecord::new(
                String::from_utf8_lossy(record.id()).to_string(),
                record.seq().to_vec(),
                record.qual().map(|q| q.to_vec()).unwrap_or_default(),
            ));
        }
        Ok(records)
    }

    /// Process records in parallel with a given function
    pub fn process_parallel<F, T>(&mut self, chunk_size: usize, f: F) -> Result<Vec<T>>
    where
        F: Fn(&FastqRecord) -> T + Send + Sync,
        T: Send,
    {
        let records = self.read_all()?;
        let results: Vec<T> = records
            .par_chunks(chunk_size)
            .flat_map(|chunk| chunk.iter().map(&f).collect::<Vec<_>>())
            .collect();
        Ok(results)
    }
}

impl Iterator for FastqParser {
    type Item = Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        self.reader.next().map(|result| {
            result
                .map(|record| {
                    FastqRecord::new(
                        String::from_utf8_lossy(record.id()).to_string(),
                        record.seq().to_vec(),
                        record.qual().map(|q| q.to_vec()).unwrap_or_default(),
                    )
                })
                .map_err(|e| Error::FastqParse(format!("Failed to read record: {}", e)))
        })
    }
}

/// Parse paired-end FASTQ files together
pub struct PairedFastqParser {
    r1_parser: FastqParser,
    r2_parser: FastqParser,
}

impl PairedFastqParser {
    pub fn open<P: AsRef<Path>>(r1_path: P, r2_path: P) -> Result<Self> {
        Ok(Self {
            r1_parser: FastqParser::open(r1_path)?,
            r2_parser: FastqParser::open(r2_path)?,
        })
    }
}

impl Iterator for PairedFastqParser {
    type Item = Result<(FastqRecord, FastqRecord)>;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.r1_parser.next(), self.r2_parser.next()) {
            (Some(Ok(r1)), Some(Ok(r2))) => Some(Ok((r1, r2))),
            (Some(Err(e)), _) | (_, Some(Err(e))) => Some(Err(e)),
            (None, None) => None,
            _ => Some(Err(Error::FastqParse(
                "Paired FASTQ files have different lengths".to_string(),
            ))),
        }
    }
}
