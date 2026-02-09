"""
Visualization helpers for single-cell data.
"""

from typing import Optional, Union

import numpy as np

try:
    import matplotlib.pyplot as plt
    import plotly.express as px
    import plotly.graph_objects as go
    _PLOTTING_AVAILABLE = True
except ImportError:
    _PLOTTING_AVAILABLE = False

try:
    import anndata as ad
    import scanpy as sc
    _SCANPY_AVAILABLE = True
except ImportError:
    _SCANPY_AVAILABLE = False


def plot_qc_violin(
    adata: "ad.AnnData",
    keys: Optional[list[str]] = None,
    groupby: Optional[str] = None,
    save: Optional[str] = None,
) -> None:
    """
    Plot QC metrics as violin plots.

    Parameters
    ----------
    adata : AnnData
        AnnData object with QC metrics
    keys : list of str, optional
        Metrics to plot (default: n_genes_by_counts, total_counts, pct_counts_mt)
    groupby : str, optional
        Column to group by
    save : str, optional
        Path to save figure
    """
    if not _SCANPY_AVAILABLE:
        raise ImportError("scanpy is required")

    if keys is None:
        keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

    # Filter to existing keys
    keys = [k for k in keys if k in adata.obs.columns]

    sc.pl.violin(adata, keys, groupby=groupby, save=save)


def plot_qc_scatter(
    adata: "ad.AnnData",
    x: str = "total_counts",
    y: str = "n_genes_by_counts",
    color: str = "pct_counts_mt",
    save: Optional[str] = None,
) -> None:
    """
    Plot QC metrics as scatter plot.

    Parameters
    ----------
    adata : AnnData
        AnnData object with QC metrics
    x : str
        X-axis metric
    y : str
        Y-axis metric
    color : str
        Color metric
    save : str, optional
        Path to save figure
    """
    if not _PLOTTING_AVAILABLE:
        raise ImportError("matplotlib and plotly are required")

    fig, ax = plt.subplots(figsize=(8, 6))

    scatter = ax.scatter(
        adata.obs[x],
        adata.obs[y],
        c=adata.obs[color] if color in adata.obs.columns else None,
        cmap="viridis",
        s=1,
        alpha=0.5,
    )

    ax.set_xlabel(x)
    ax.set_ylabel(y)

    if color in adata.obs.columns:
        plt.colorbar(scatter, label=color)

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    else:
        plt.show()


def plot_umap_interactive(
    adata: "ad.AnnData",
    color: str = "leiden",
    hover_data: Optional[list[str]] = None,
    title: Optional[str] = None,
) -> "go.Figure":
    """
    Create interactive UMAP plot with Plotly.

    Parameters
    ----------
    adata : AnnData
        AnnData object with UMAP coordinates
    color : str
        Column to color by
    hover_data : list of str, optional
        Additional columns to show on hover
    title : str, optional
        Plot title

    Returns
    -------
    plotly.graph_objects.Figure
        Interactive Plotly figure
    """
    if not _PLOTTING_AVAILABLE:
        raise ImportError("plotly is required")

    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP not found. Run sc.tl.umap first.")

    umap_coords = adata.obsm["X_umap"]

    # Build hover data
    hover_dict = {"UMAP1": umap_coords[:, 0], "UMAP2": umap_coords[:, 1]}

    if color in adata.obs.columns:
        hover_dict[color] = adata.obs[color].values

    if hover_data:
        for col in hover_data:
            if col in adata.obs.columns:
                hover_dict[col] = adata.obs[col].values

    import pandas as pd
    df = pd.DataFrame(hover_dict)

    fig = px.scatter(
        df,
        x="UMAP1",
        y="UMAP2",
        color=color if color in df.columns else None,
        hover_data=list(hover_dict.keys()),
        title=title or f"UMAP colored by {color}",
    )

    fig.update_traces(marker=dict(size=3, opacity=0.7))
    fig.update_layout(
        template="plotly_white",
        width=800,
        height=600,
    )

    return fig


def plot_gene_expression(
    adata: "ad.AnnData",
    genes: Union[str, list[str]],
    use_raw: bool = True,
) -> None:
    """
    Plot gene expression on UMAP.

    Parameters
    ----------
    adata : AnnData
        AnnData object
    genes : str or list of str
        Gene(s) to plot
    use_raw : bool
        Whether to use raw data (default: True)
    """
    if not _SCANPY_AVAILABLE:
        raise ImportError("scanpy is required")

    if isinstance(genes, str):
        genes = [genes]

    sc.pl.umap(adata, color=genes, use_raw=use_raw)


def plot_knee(
    counts_per_cell: np.ndarray,
    expected_cells: Optional[int] = None,
    save: Optional[str] = None,
) -> None:
    """
    Plot knee plot for cell calling.

    Parameters
    ----------
    counts_per_cell : array
        UMI counts per cell barcode
    expected_cells : int, optional
        Expected number of cells
    save : str, optional
        Path to save figure
    """
    if not _PLOTTING_AVAILABLE:
        raise ImportError("matplotlib is required")

    sorted_counts = np.sort(counts_per_cell)[::-1]
    ranks = np.arange(1, len(sorted_counts) + 1)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.loglog(ranks, sorted_counts, linewidth=1)
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel("UMI counts")
    ax.set_title("Knee plot")

    if expected_cells is not None and expected_cells < len(sorted_counts):
        ax.axvline(expected_cells, color="red", linestyle="--", label=f"Expected: {expected_cells}")
        ax.legend()

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")
    else:
        plt.show()
