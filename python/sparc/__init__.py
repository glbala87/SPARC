"""
SPARC - Single-cell Pipeline Accelerated in Rust Core

This package provides Python bindings for the Rust-based SPARC library,
with seamless integration with Scanpy and AnnData.
"""

from typing import TYPE_CHECKING

__version__ = "0.1.0"

# Import Rust bindings
try:
    from sparc._sparc_py import (
        FastqParser,
        FastqRecord,
        FastqWriter,
        BamParser,
        BamRecord,
        Whitelist,
        BarcodeCorrector,
        CountMatrix,
        GeneCounter,
    )

    # Also import new bindings when available
    from sparc._sparc_py import (
        QcMetrics,
        QcReport,
        py_normalize_total as rust_normalize_total,
        py_pca as rust_pca,
        py_run_analysis as rust_run_analysis,
        py_label_propagation as rust_label_propagation,
        adjusted_rand_index as rust_adjusted_rand_index,
        normalized_mutual_info as rust_normalized_mutual_info,
    )

    _RUST_AVAILABLE = True
except ImportError:
    _RUST_AVAILABLE = False

# Import Python modules
from sparc.io import read_fastq, read_bam, read_matrix, write_matrix, write_h5ad, read_h5ad
from sparc.preprocessing import extract_barcodes, correct_barcodes, deduplicate_umis
from sparc.analysis import to_anndata, from_anndata, run_pipeline, normalize_and_analyze, find_marker_genes
from sparc.streaming import StreamingProcessor, StreamStats
from sparc.batch import BatchProcessor, BatchSummary, parse_manifest
from sparc import validation

if TYPE_CHECKING:
    from sparc._sparc_py import (
        FastqParser,
        FastqRecord,
        FastqWriter,
        BamParser,
        BamRecord,
        Whitelist,
        BarcodeCorrector,
        CountMatrix,
        GeneCounter,
    )

__all__ = [
    # Version
    "__version__",
    # Rust classes
    "FastqParser",
    "FastqRecord",
    "FastqWriter",
    "BamParser",
    "BamRecord",
    "Whitelist",
    "BarcodeCorrector",
    "CountMatrix",
    "GeneCounter",
    # I/O functions
    "read_fastq",
    "read_bam",
    "read_matrix",
    "write_matrix",
    "write_h5ad",
    "read_h5ad",
    # Preprocessing
    "extract_barcodes",
    "correct_barcodes",
    "deduplicate_umis",
    # Analysis
    "to_anndata",
    "from_anndata",
    "run_pipeline",
    "normalize_and_analyze",
    "find_marker_genes",
    # Streaming
    "StreamingProcessor",
    "StreamStats",
    # Batch
    "BatchProcessor",
    "BatchSummary",
    "parse_manifest",
    # Validation
    "validation",
    # Utilities
    "check_rust_bindings",
]


def check_rust_bindings() -> bool:
    """Check if Rust bindings are available."""
    return _RUST_AVAILABLE
