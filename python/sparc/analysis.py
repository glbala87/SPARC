"""
Analysis functions and Scanpy integration.
"""

from pathlib import Path
from typing import Optional, Union
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

try:
    import anndata as ad
    import scanpy as sc
    _SCANPY_AVAILABLE = True
except ImportError:
    _SCANPY_AVAILABLE = False

try:
    from sparc._sparc_py import CountMatrix, GeneCounter
    _RUST_AVAILABLE = True
except ImportError:
    _RUST_AVAILABLE = False


def to_anndata(
    matrix: Union["CountMatrix", sp.spmatrix],
    barcodes: Optional[list[str]] = None,
    genes: Optional[list[str]] = None,
) -> "ad.AnnData":
    """
    Convert a count matrix to AnnData format.

    Parameters
    ----------
    matrix : CountMatrix or sparse matrix
        Count matrix (genes x cells or cells x genes)
    barcodes : list of str, optional
        Cell barcodes (required if matrix is sparse)
    genes : list of str, optional
        Gene names (required if matrix is sparse)

    Returns
    -------
    AnnData
        AnnData object ready for Scanpy analysis
    """
    if not _SCANPY_AVAILABLE:
        raise ImportError("scanpy and anndata are required. Install with: pip install scanpy")

    if _RUST_AVAILABLE and isinstance(matrix, CountMatrix):
        # Convert from Rust CountMatrix
        barcodes = matrix.barcodes
        genes = matrix.genes

        # Build sparse matrix
        data = np.array(matrix.values(), dtype=np.float32)
        row = np.array(matrix.row_indices())
        col = np.array(matrix.col_indices())

        sparse_mat = sp.csr_matrix(
            (data, (col, row)),  # cells x genes
            shape=(matrix.n_cols, matrix.n_rows),
        )
    elif sp.issparse(matrix):
        sparse_mat = matrix
        if barcodes is None or genes is None:
            raise ValueError("barcodes and genes are required for sparse matrix input")
    else:
        raise TypeError(f"Unsupported matrix type: {type(matrix)}")

    # Create AnnData
    adata = ad.AnnData(
        X=sparse_mat,
        obs={"barcode": barcodes},
        var={"gene": genes},
    )
    adata.obs_names = barcodes
    adata.var_names = genes

    return adata


def from_anndata(adata: "ad.AnnData") -> "CountMatrix":
    """
    Convert AnnData to CountMatrix format.

    Parameters
    ----------
    adata : AnnData
        AnnData object

    Returns
    -------
    CountMatrix
        Rust CountMatrix object
    """
    if not _RUST_AVAILABLE:
        raise ImportError("Rust bindings not available")

    counter = GeneCounter()

    # Get sparse matrix
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    barcodes = list(adata.obs_names)
    genes = list(adata.var_names)

    # Populate counter
    cx = X.tocoo()
    for i, j, v in zip(cx.row, cx.col, cx.data):
        if v > 0:
            counter.add_count(barcodes[i], genes[j], int(v))

    return counter.build()


@dataclass
class PipelineResult:
    """Result from running the analysis pipeline."""
    adata: "ad.AnnData"
    qc_passed: int
    qc_failed: int
    total_genes: int
    total_cells: int


def run_pipeline(
    matrix_path: Union[str, Path],
    min_genes: int = 200,
    max_genes: int = 10000,
    max_mito: float = 20.0,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    resolution: float = 1.0,
) -> PipelineResult:
    """
    Run a standard single-cell analysis pipeline.

    Parameters
    ----------
    matrix_path : str or Path
        Path to matrix directory
    min_genes : int
        Minimum genes per cell (default: 200)
    max_genes : int
        Maximum genes per cell (default: 10000)
    max_mito : float
        Maximum mitochondrial percentage (default: 20.0)
    n_pcs : int
        Number of principal components (default: 50)
    n_neighbors : int
        Number of neighbors for graph construction (default: 15)
    resolution : float
        Clustering resolution (default: 1.0)

    Returns
    -------
    PipelineResult
        Analysis results including processed AnnData
    """
    if not _SCANPY_AVAILABLE:
        raise ImportError("scanpy is required. Install with: pip install scanpy")

    from sparc.io import read_matrix

    # Load data
    matrix, barcodes, genes = read_matrix(matrix_path)
    adata = to_anndata(matrix, barcodes, genes)

    initial_cells = adata.n_obs

    # QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    adata = adata[adata.obs.pct_counts_mt < max_mito, :]

    qc_passed = adata.n_obs
    qc_failed = initial_cells - qc_passed

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=3)

    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    # Scale
    sc.pp.scale(adata, max_value=10)

    # PCA
    sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1))

    # Neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))

    # UMAP
    sc.tl.umap(adata)

    # Clustering
    sc.tl.leiden(adata, resolution=resolution)

    return PipelineResult(
        adata=adata,
        qc_passed=qc_passed,
        qc_failed=qc_failed,
        total_genes=adata.n_vars,
        total_cells=adata.n_obs,
    )


def normalize_and_analyze(
    adata: "ad.AnnData",
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    n_pcs: int = 50,
    n_neighbors: int = 15,
    resolution: float = 1.0,
) -> "ad.AnnData":
    """
    Run normalization and downstream analysis on an AnnData object.

    Performs: normalize_total -> log1p -> HVG selection -> scale -> PCA ->
    neighbors -> UMAP -> Leiden clustering.

    Parameters
    ----------
    adata : AnnData
        Input AnnData (raw counts)
    target_sum : float
        Target sum for normalization (default: 1e4)
    n_top_genes : int
        Number of highly variable genes to select (default: 2000)
    n_pcs : int
        Number of principal components (default: 50)
    n_neighbors : int
        Number of neighbors for graph construction (default: 15)
    resolution : float
        Leiden clustering resolution (default: 1.0)

    Returns
    -------
    AnnData
        Processed AnnData with embeddings and clusters
    """
    if not _SCANPY_AVAILABLE:
        raise ImportError("scanpy is required. Install with: pip install scanpy")

    adata = adata.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    return adata


def find_marker_genes(
    adata: "ad.AnnData",
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 25,
) -> "ad.AnnData":
    """
    Find marker genes for each cluster.

    Parameters
    ----------
    adata : AnnData
        Processed AnnData with clusters
    groupby : str
        Column in obs to group by (default: "leiden")
    method : str
        Statistical method (default: "wilcoxon")
    n_genes : int
        Number of top genes per group (default: 25)

    Returns
    -------
    AnnData
        AnnData with rank_genes_groups result
    """
    if not _SCANPY_AVAILABLE:
        raise ImportError("scanpy is required. Install with: pip install scanpy")

    use_raw = adata.raw is not None
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method, n_genes=n_genes, use_raw=use_raw)

    return adata
