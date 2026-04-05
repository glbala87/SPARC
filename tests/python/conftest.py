"""Shared test fixtures for SPARC Python tests."""

import tempfile
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp


@pytest.fixture
def tmp_dir():
    """Provide a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as d:
        yield Path(d)


@pytest.fixture
def sample_matrix():
    """Create a sample sparse count matrix (cells x genes)."""
    data = np.array([10, 3, 2, 5, 8, 7], dtype=np.float32)
    row = np.array([0, 0, 0, 1, 1, 1])  # cells
    col = np.array([0, 1, 2, 0, 1, 2])  # genes
    matrix = sp.csr_matrix((data, (row, col)), shape=(2, 3))
    barcodes = ["CELL1", "CELL2"]
    genes = ["GENE1", "GENE2", "GENE3"]
    return matrix, barcodes, genes


@pytest.fixture
def sample_mtx_dir(tmp_dir, sample_matrix):
    """Create a sample 10x-style matrix directory."""
    from scipy.io import mmwrite

    matrix, barcodes, genes = sample_matrix
    mtx_dir = tmp_dir / "matrix"
    mtx_dir.mkdir()

    mmwrite(str(mtx_dir / "matrix.mtx"), matrix.T.tocoo())

    with open(mtx_dir / "barcodes.tsv", "w") as f:
        for bc in barcodes:
            f.write(f"{bc}\n")

    with open(mtx_dir / "genes.tsv", "w") as f:
        for gene in genes:
            f.write(f"{gene}\t{gene}\n")

    return mtx_dir, barcodes, genes
