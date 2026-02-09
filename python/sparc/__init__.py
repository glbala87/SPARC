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

    _RUST_AVAILABLE = True
except ImportError:
    _RUST_AVAILABLE = False

# Import Python modules
from sparc.io import read_fastq, read_bam, read_matrix, write_matrix
from sparc.preprocessing import extract_barcodes, correct_barcodes, deduplicate_umis
from sparc.analysis import to_anndata, from_anndata, run_pipeline

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
    # Preprocessing
    "extract_barcodes",
    "correct_barcodes",
    "deduplicate_umis",
    # Analysis
    "to_anndata",
    "from_anndata",
    "run_pipeline",
]


def check_rust_bindings() -> bool:
    """Check if Rust bindings are available."""
    return _RUST_AVAILABLE
