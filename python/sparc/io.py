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

    if not path.exists():
        raise FileNotFoundError(f"Matrix directory not found: {path}")
    if not path.is_dir():
        raise NotADirectoryError(f"Expected a directory, got a file: {path}")

    # Read barcodes
    barcodes_file = path / "barcodes.tsv"
    if not barcodes_file.exists():
        barcodes_file = path / "barcodes.tsv.gz"
    if not barcodes_file.exists():
        raise FileNotFoundError(
            f"barcodes.tsv (or .gz) not found in {path}. "
            f"Available files: {[f.name for f in path.iterdir()]}"
        )

    opener = _open_maybe_gz(barcodes_file)
    with opener as f:
        barcodes = [line.strip() for line in f if line.strip()]

    # Read genes
    genes_file = path / "genes.tsv"
    if not genes_file.exists():
        genes_file = path / "features.tsv"
    if not genes_file.exists():
        genes_file = path / "features.tsv.gz"
    if not genes_file.exists():
        genes_file = path / "genes.tsv.gz"
    if not genes_file.exists():
        raise FileNotFoundError(
            f"genes.tsv/features.tsv (or .gz) not found in {path}. "
            f"Available files: {[f.name for f in path.iterdir()]}"
        )

    genes = []
    opener = _open_maybe_gz(genes_file)
    with opener as f:
        for line in f:
            parts = line.strip().split("\t")
            genes.append(parts[0] if parts else "")

    # Read matrix
    mtx_file = path / "matrix.mtx"
    if not mtx_file.exists():
        mtx_file = path / "matrix.mtx.gz"
    if not mtx_file.exists():
        raise FileNotFoundError(
            f"matrix.mtx (or .gz) not found in {path}. "
            f"Available files: {[f.name for f in path.iterdir()]}"
        )

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


def write_h5ad(
    matrix: sp.spmatrix,
    barcodes: list[str],
    genes: list[str],
    path: Union[str, Path],
    gene_ids: Optional[list[str]] = None,
) -> None:
    """
    Write a count matrix to H5AD (AnnData) format.

    Parameters
    ----------
    matrix : sparse matrix
        Count matrix (cells x genes)
    barcodes : list of str
        Cell barcodes
    genes : list of str
        Gene names
    path : str or Path
        Output .h5ad file path
    gene_ids : list of str, optional
        Gene IDs (defaults to gene names)
    """
    try:
        import anndata as ad
    except ImportError:
        raise ImportError("anndata is required for H5AD export. Install with: pip install anndata")

    import pandas as pd

    obs = pd.DataFrame(index=barcodes)
    var_data = {"gene_name": genes}
    if gene_ids:
        var_data["gene_id"] = gene_ids
    var = pd.DataFrame(var_data, index=genes if not gene_ids else gene_ids)

    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)
    else:
        matrix = matrix.tocsr()

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.write_h5ad(str(path))


def read_h5ad(
    path: Union[str, Path],
) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """
    Read a count matrix from H5AD format.

    Parameters
    ----------
    path : str or Path
        Path to .h5ad file

    Returns
    -------
    tuple
        (sparse_matrix, barcodes, genes)
    """
    try:
        import anndata as ad
    except ImportError:
        raise ImportError("anndata is required for H5AD import. Install with: pip install anndata")

    adata = ad.read_h5ad(str(path))
    matrix = adata.X
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)

    barcodes = list(adata.obs_names)
    genes = list(adata.var_names)

    return matrix, barcodes, genes


def _open_maybe_gz(path: Path):
    """Open a file, auto-detecting gzip compression."""
    import gzip
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)
