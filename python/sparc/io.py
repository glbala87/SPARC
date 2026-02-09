"""
I/O functions for reading and writing single-cell data.
"""

from pathlib import Path
from typing import Iterator, Optional, Union

import numpy as np
import scipy.sparse as sp

try:
    from sparc._sparc_py import (
        FastqParser,
        FastqRecord,
        BamParser,
        BamRecord,
        CountMatrix,
    )
    _RUST_AVAILABLE = True
except ImportError:
    _RUST_AVAILABLE = False


def read_fastq(path: Union[str, Path]) -> Iterator["FastqRecord"]:
    """
    Read a FASTQ file and iterate over records.

    Parameters
    ----------
    path : str or Path
        Path to FASTQ file (supports .gz and .zst compression)

    Yields
    ------
    FastqRecord
        FASTQ record with id, seq, and qual attributes
    """
    if not _RUST_AVAILABLE:
        raise ImportError("Rust bindings not available. Install with: pip install sparc")

    parser = FastqParser(str(path))
    for record in parser:
        yield record


def read_bam(
    path: Union[str, Path],
    min_mapq: int = 0,
) -> Iterator["BamRecord"]:
    """
    Read a BAM file and iterate over records.

    Parameters
    ----------
    path : str or Path
        Path to BAM file
    min_mapq : int, optional
        Minimum mapping quality filter (default: 0)

    Yields
    ------
    BamRecord
        BAM record with alignment information and tags
    """
    if not _RUST_AVAILABLE:
        raise ImportError("Rust bindings not available. Install with: pip install sparc")

    parser = BamParser(str(path))
    for record in parser:
        if record.mapq >= min_mapq:
            yield record


def read_matrix(
    path: Union[str, Path],
    genome: Optional[str] = None,
) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """
    Read a count matrix from 10x-style directory.

    Parameters
    ----------
    path : str or Path
        Path to matrix directory containing matrix.mtx, barcodes.tsv, genes.tsv
    genome : str, optional
        Genome name for multi-genome matrices

    Returns
    -------
    tuple
        (sparse_matrix, barcodes, genes)
    """
    path = Path(path)

    # Read barcodes
    barcodes_file = path / "barcodes.tsv"
    if not barcodes_file.exists():
        barcodes_file = path / "barcodes.tsv.gz"

    with open(barcodes_file) as f:
        barcodes = [line.strip() for line in f]

    # Read genes
    genes_file = path / "genes.tsv"
    if not genes_file.exists():
        genes_file = path / "features.tsv"
    if not genes_file.exists():
        genes_file = path / "genes.tsv.gz"

    genes = []
    with open(genes_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            genes.append(parts[0] if parts else "")

    # Read matrix
    mtx_file = path / "matrix.mtx"
    if not mtx_file.exists():
        mtx_file = path / "matrix.mtx.gz"

    from scipy.io import mmread
    matrix = mmread(str(mtx_file)).T.tocsr()

    return matrix, barcodes, genes


def write_matrix(
    matrix: sp.spmatrix,
    barcodes: list[str],
    genes: list[str],
    path: Union[str, Path],
    compress: bool = False,
) -> None:
    """
    Write a count matrix to 10x-style directory.

    Parameters
    ----------
    matrix : sparse matrix
        Count matrix (cells x genes)
    barcodes : list of str
        Cell barcodes
    genes : list of str
        Gene names/IDs
    path : str or Path
        Output directory
    compress : bool, optional
        Whether to gzip output files (default: False)
    """
    from scipy.io import mmwrite
    import gzip

    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    # Write matrix (transpose to genes x cells for 10x format)
    mtx_path = path / "matrix.mtx"
    mmwrite(str(mtx_path), matrix.T.tocoo())

    # Write barcodes
    barcodes_path = path / "barcodes.tsv"
    with open(barcodes_path, "w") as f:
        for bc in barcodes:
            f.write(f"{bc}\n")

    # Write genes
    genes_path = path / "genes.tsv"
    with open(genes_path, "w") as f:
        for gene in genes:
            f.write(f"{gene}\t{gene}\n")

    if compress:
        for file_path in [mtx_path, barcodes_path, genes_path]:
            with open(file_path, "rb") as f_in:
                with gzip.open(f"{file_path}.gz", "wb") as f_out:
                    f_out.write(f_in.read())
            file_path.unlink()
