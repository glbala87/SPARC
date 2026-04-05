//! Streaming/chunked processing for large datasets

use crate::fastq::{FastqParser, FastqRecord};
use crate::Result;
use std::path::Path;

/// Configuration for streaming processing
pub struct StreamConfig {
    /// Number of records per chunk
    pub chunk_size: usize,
    /// Maximum memory usage in bytes (advisory)
    pub max_memory: usize,
}

impl Default for StreamConfig {
    fn default() -> Self {
        Self {
            chunk_size: 1_000_000,
            max_memory: 4 * 1024 * 1024 * 1024, // 4GB
        }
    }
}

/// Statistics from streaming processing
#[derive(Debug, Default)]
pub struct StreamStats {
    pub total_records: u64,
    pub chunks_processed: u64,
}

/// Streaming FASTQ processor that processes records in fixed-size chunks
pub struct StreamingProcessor {
    config: StreamConfig,
}

impl StreamingProcessor {
    pub fn new(config: StreamConfig) -> Self {
        Self { config }
    }

    pub fn with_chunk_size(chunk_size: usize) -> Self {
        Self {
            config: StreamConfig {
                chunk_size,
                ..Default::default()
            },
        }
    }

    /// Process FASTQ file in chunks, calling the callback for each chunk
    pub fn process_fastq<P, F>(&self, path: P, mut callback: F) -> Result<StreamStats>
    where
        P: AsRef<Path>,
        F: FnMut(&[FastqRecord]) -> Result<()>,
    {
        let mut parser = FastqParser::open(path)?;
        let mut stats = StreamStats::default();
        let mut chunk = Vec::with_capacity(self.config.chunk_size);

        loop {
            chunk.clear();
            for _ in 0..self.config.chunk_size {
                match parser.next() {
                    Some(Ok(record)) => chunk.push(record),
                    Some(Err(e)) => return Err(e),
                    None => break,
                }
            }
            if chunk.is_empty() {
                break;
            }
            stats.total_records += chunk.len() as u64;
            stats.chunks_processed += 1;
            callback(&chunk)?;
        }
        Ok(stats)
    }

    /// Process FASTQ file in parallel chunks using rayon
    pub fn process_fastq_parallel<P, F, T>(&self, path: P, f: F) -> Result<Vec<T>>
    where
        P: AsRef<Path>,
        F: Fn(&FastqRecord) -> T + Send + Sync,
        T: Send,
    {
        use rayon::prelude::*;

        let mut parser = FastqParser::open(path)?;
        let records = parser.read_all()?;

        let results: Vec<T> = records
            .par_chunks(self.config.chunk_size)
            .flat_map(|chunk| chunk.iter().map(&f).collect::<Vec<_>>())
            .collect();

        Ok(results)
    }
}
