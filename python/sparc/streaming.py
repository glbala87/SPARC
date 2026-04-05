"""
Streaming / memory-efficient processing for large FASTQ files.
"""

from pathlib import Path
from typing import Callable, Iterator, Optional, Union
from dataclasses import dataclass

try:
    from sparc._sparc_py import FastqParser, FastqRecord
    _RUST_AVAILABLE = True
except ImportError:
    _RUST_AVAILABLE = False


@dataclass
class StreamStats:
    """Statistics from streaming processing."""
    total_records: int = 0
    processed_records: int = 0
    chunks_processed: int = 0


class StreamingProcessor:
    """
    Memory-efficient streaming processor for large FASTQ files.

    Processes data in configurable chunk sizes to limit memory usage.

    Parameters
    ----------
    chunk_size : int
        Number of records per chunk (default: 100000)
    max_memory_mb : int
        Approximate memory limit in MB (default: 4096)
    """

    def __init__(self, chunk_size: int = 100_000, max_memory_mb: int = 4096):
        self.chunk_size = chunk_size
        self.max_memory_mb = max_memory_mb

    def process_fastq(
        self,
        path: Union[str, Path],
        callback: Callable[[list], None],
    ) -> StreamStats:
        """
        Process a FASTQ file in chunks.

        Parameters
        ----------
        path : str or Path
            Path to FASTQ file
        callback : callable
            Function called with each chunk of records

        Returns
        -------
        StreamStats
            Processing statistics
        """
        if not _RUST_AVAILABLE:
            raise ImportError("Rust bindings not available. Install with: pip install sparc")

        stats = StreamStats()
        parser = FastqParser(str(path))
        chunk = []

        for record in parser:
            chunk.append(record)
            stats.total_records += 1

            if len(chunk) >= self.chunk_size:
                callback(chunk)
                stats.chunks_processed += 1
                stats.processed_records += len(chunk)
                chunk = []

        if chunk:
            callback(chunk)
            stats.chunks_processed += 1
            stats.processed_records += len(chunk)

        return stats

    def iter_chunks(
        self,
        path: Union[str, Path],
    ) -> Iterator[list]:
        """
        Iterate over chunks of FASTQ records.

        Parameters
        ----------
        path : str or Path
            Path to FASTQ file

        Yields
        ------
        list
            Chunk of FastqRecord objects
        """
        if not _RUST_AVAILABLE:
            raise ImportError("Rust bindings not available. Install with: pip install sparc")

        parser = FastqParser(str(path))
        chunk = []

        for record in parser:
            chunk.append(record)
            if len(chunk) >= self.chunk_size:
                yield chunk
                chunk = []

        if chunk:
            yield chunk
