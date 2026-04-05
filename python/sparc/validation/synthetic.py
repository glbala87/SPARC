"""
Synthetic ground-truth dataset generation for validation.

Generates scRNA-seq data with known cell types, expression profiles,
and barcode sequences using Poisson-distributed counts.
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import scipy.sparse as sp

try:
    import anndata as ad
    _ANNDATA_AVAILABLE = True
except ImportError:
    _ANNDATA_AVAILABLE = False


BASES = list("ACGT")


@dataclass
class SyntheticConfig:
    """Configuration for synthetic dataset generation."""
    n_cells: int = 500
    n_genes: int = 200
    n_cell_types: int = 5
    n_markers_per_type: int = 10
    marker_fold_change: float = 8.0
    base_expression: float = 2.0
    mutation_rate: float = 0.1
    invalid_barcode_rate: float = 0.05
    barcode_len: int = 16
    umi_len: int = 12
    seed: int = 42
    protocol: str = "10x-3prime-v3"


@dataclass
class SyntheticDataset:
    """A complete synthetic dataset with ground truth."""
    # Ground truth count matrix (cells x genes) as dense numpy array
    count_matrix: np.ndarray
    # Cell barcodes
    barcodes: list[str]
    # Gene names
    genes: list[str]
    # Cell type assignment per cell (index)
    cell_type_labels: np.ndarray
    # Cell type names
    cell_type_names: list[str]
    # Mean expression profiles per cell type (n_types x n_genes)
    expression_profiles: np.ndarray
    # Mutated barcodes: list of (original, mutated, distance)
    mutated_barcodes: list[tuple[str, str, int]] = field(default_factory=list)
    # Invalid barcodes
    invalid_barcodes: list[str] = field(default_factory=list)
    # Config used
    config: Optional[SyntheticConfig] = None

    def to_anndata(self) -> "ad.AnnData":
        """Convert to AnnData with ground-truth cell type labels in obs."""
        if not _ANNDATA_AVAILABLE:
            raise ImportError("anndata is required: pip install anndata")

        adata = ad.AnnData(
            X=sp.csr_matrix(self.count_matrix.astype(np.float32)),
            obs={
                "barcode": self.barcodes,
                "cell_type": [self.cell_type_names[i] for i in self.cell_type_labels],
                "cell_type_id": self.cell_type_labels.astype(int),
            },
            var={"gene": self.genes},
        )
        adata.obs_names = self.barcodes
        adata.var_names = self.genes
        return adata


def generate_synthetic_dataset(config: Optional[SyntheticConfig] = None) -> SyntheticDataset:
    """
    Generate a synthetic scRNA-seq dataset with known ground truth.

    Parameters
    ----------
    config : SyntheticConfig, optional
        Generation parameters. Uses defaults if not provided.

    Returns
    -------
    SyntheticDataset
        Dataset with ground-truth labels and expression matrix.
    """
    if config is None:
        config = SyntheticConfig()

    rng = np.random.default_rng(config.seed)

    # Generate unique barcodes
    barcodes = _generate_barcodes(rng, config.n_cells, config.barcode_len)

    # Gene names
    genes = [f"GENE_{i:04d}" for i in range(config.n_genes)]

    # Cell type assignments (round-robin)
    cell_type_labels = np.array([i % config.n_cell_types for i in range(config.n_cells)])
    cell_type_names = [f"CellType_{i}" for i in range(config.n_cell_types)]

    # Build expression profiles
    expression_profiles = _build_expression_profiles(config, rng)

    # Sample counts from Poisson
    count_matrix = np.zeros((config.n_cells, config.n_genes), dtype=np.float64)
    for cell_idx in range(config.n_cells):
        ct = cell_type_labels[cell_idx]
        means = expression_profiles[ct]
        count_matrix[cell_idx] = rng.poisson(lam=np.maximum(means, 0.01))

    # Generate mutated barcodes
    n_mutate = int(config.n_cells * config.mutation_rate)
    mutated_barcodes = []
    for i in range(min(n_mutate, config.n_cells)):
        original = barcodes[i]
        mutated = list(original)
        pos = rng.integers(0, len(mutated))
        orig_base = mutated[pos]
        new_base = orig_base
        while new_base == orig_base:
            new_base = BASES[rng.integers(0, 4)]
        mutated[pos] = new_base
        mutated_barcodes.append((original, "".join(mutated), 1))

    # Generate invalid barcodes
    n_invalid = int(config.n_cells * config.invalid_barcode_rate)
    barcode_set = set(barcodes)
    invalid_barcodes = []
    while len(invalid_barcodes) < n_invalid:
        bc = "".join(
            rng.choice(list("ACGTN"), size=config.barcode_len)
        )
        if bc not in barcode_set:
            invalid_barcodes.append(bc)

    return SyntheticDataset(
        count_matrix=count_matrix,
        barcodes=barcodes,
        genes=genes,
        cell_type_labels=cell_type_labels,
        cell_type_names=cell_type_names,
        expression_profiles=expression_profiles,
        mutated_barcodes=mutated_barcodes,
        invalid_barcodes=invalid_barcodes,
        config=config,
    )


def _generate_barcodes(rng: np.random.Generator, n: int, length: int) -> list[str]:
    """Generate n unique random DNA barcodes."""
    seen = set()
    barcodes = []
    while len(barcodes) < n:
        bc = "".join(rng.choice(BASES, size=length))
        if bc not in seen:
            seen.add(bc)
            barcodes.append(bc)
    return barcodes


def _build_expression_profiles(
    config: SyntheticConfig, rng: np.random.Generator
) -> np.ndarray:
    """Build mean expression profile for each cell type."""
    profiles = np.full(
        (config.n_cell_types, config.n_genes),
        config.base_expression,
        dtype=np.float64,
    )

    for ct in range(config.n_cell_types):
        marker_start = ct * config.n_markers_per_type
        marker_end = min(marker_start + config.n_markers_per_type, config.n_genes)
        profiles[ct, marker_start:marker_end] = (
            config.base_expression * config.marker_fold_change
        )

    # Add random noise
    noise = 1.0 + rng.uniform(-0.3, 0.3, size=profiles.shape)
    profiles = np.maximum(profiles * noise, 0.1)

    return profiles
