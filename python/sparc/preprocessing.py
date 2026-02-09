"""
Preprocessing utilities for single-cell data.
"""

from pathlib import Path
from typing import Optional, Union
from dataclasses import dataclass

import numpy as np

try:
    from sparc._sparc_py import (
        FastqParser,
        Whitelist,
        BarcodeCorrector,
        GeneCounter,
    )
    _RUST_AVAILABLE = True
except ImportError:
    _RUST_AVAILABLE = False


@dataclass
class ExtractionResult:
    """Result of barcode/UMI extraction."""
    total_reads: int
    valid_barcodes: int
    corrected_barcodes: int
    valid_umis: int
    barcodes: list[str]
    umis: list[str]


def extract_barcodes(
    r1_path: Union[str, Path],
    whitelist_path: Union[str, Path],
    barcode_start: int = 0,
    barcode_len: int = 16,
    umi_start: int = 16,
    umi_len: int = 12,
    max_mismatch: int = 1,
    min_quality: int = 10,
) -> ExtractionResult:
    """
    Extract barcodes and UMIs from FASTQ reads.

    Parameters
    ----------
    r1_path : str or Path
        Path to R1 FASTQ file
    whitelist_path : str or Path
        Path to barcode whitelist file
    barcode_start : int
        Start position of barcode in read (default: 0)
    barcode_len : int
        Length of barcode (default: 16)
    umi_start : int
        Start position of UMI in read (default: 16)
    umi_len : int
        Length of UMI (default: 12)
    max_mismatch : int
        Maximum Hamming distance for barcode correction (default: 1)
    min_quality : int
        Minimum mean quality score for barcode region (default: 10)

    Returns
    -------
    ExtractionResult
        Extraction statistics and lists of valid barcodes/UMIs
    """
    if not _RUST_AVAILABLE:
        raise ImportError("Rust bindings not available")

    whitelist = Whitelist(str(whitelist_path))
    corrector = BarcodeCorrector(whitelist, max_mismatch)

    parser = FastqParser(str(r1_path))

    total_reads = 0
    valid_barcodes = 0
    corrected_barcodes = 0
    valid_umis = 0
    barcodes = []
    umis = []

    for record in parser:
        total_reads += 1

        # Check read length
        if len(record.seq) < barcode_start + barcode_len + umi_len:
            continue

        # Extract barcode and UMI
        barcode_seq = record.subsequence(barcode_start, barcode_len)
        umi_seq = record.subsequence(umi_start, umi_len)

        if barcode_seq is None or umi_seq is None:
            continue

        barcode_str = bytes(barcode_seq).decode("ascii", errors="ignore")
        umi_str = bytes(umi_seq).decode("ascii", errors="ignore")

        # Match barcode
        status, corrected, _ = corrector.match_barcode(barcode_str)

        if status == "exact":
            valid_barcodes += 1
            barcodes.append(corrected)
            umis.append(umi_str)
            valid_umis += 1
        elif status == "corrected":
            valid_barcodes += 1
            corrected_barcodes += 1
            barcodes.append(corrected)
            umis.append(umi_str)
            valid_umis += 1

    return ExtractionResult(
        total_reads=total_reads,
        valid_barcodes=valid_barcodes,
        corrected_barcodes=corrected_barcodes,
        valid_umis=valid_umis,
        barcodes=barcodes,
        umis=umis,
    )


def correct_barcodes(
    barcodes: list[str],
    whitelist_path: Union[str, Path],
    max_mismatch: int = 1,
) -> list[Optional[str]]:
    """
    Correct barcodes against a whitelist.

    Parameters
    ----------
    barcodes : list of str
        Barcodes to correct
    whitelist_path : str or Path
        Path to barcode whitelist file
    max_mismatch : int
        Maximum Hamming distance for correction (default: 1)

    Returns
    -------
    list of str or None
        Corrected barcodes (None if no match found)
    """
    if not _RUST_AVAILABLE:
        raise ImportError("Rust bindings not available")

    whitelist = Whitelist(str(whitelist_path))
    corrector = BarcodeCorrector(whitelist, max_mismatch)

    return corrector.correct_batch(barcodes)


def deduplicate_umis(
    barcodes: list[str],
    umis: list[str],
    genes: list[str],
    max_distance: int = 1,
) -> tuple[list[str], list[str], list[str], list[int]]:
    """
    Deduplicate UMIs within each cell-gene group.

    Parameters
    ----------
    barcodes : list of str
        Cell barcodes
    umis : list of str
        UMI sequences
    genes : list of str
        Gene assignments
    max_distance : int
        Maximum edit distance for UMI clustering (default: 1)

    Returns
    -------
    tuple
        (unique_barcodes, unique_genes, representative_umis, counts)
    """
    from collections import defaultdict

    # Group by barcode-gene
    groups = defaultdict(list)
    for bc, umi, gene in zip(barcodes, umis, genes):
        groups[(bc, gene)].append(umi)

    unique_barcodes = []
    unique_genes = []
    representative_umis = []
    counts = []

    for (bc, gene), umi_list in groups.items():
        # Simple deduplication: count unique UMIs
        # For proper clustering, use the Rust implementation
        unique_umis = set(umi_list)
        unique_barcodes.append(bc)
        unique_genes.append(gene)
        representative_umis.append(list(unique_umis)[0])
        counts.append(len(unique_umis))

    return unique_barcodes, unique_genes, representative_umis, counts
