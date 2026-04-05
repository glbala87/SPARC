"""Tests for the SPARC Python package."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp


class TestIO:
    """Tests for sparc.io module."""

    def test_read_matrix(self, sample_mtx_dir):
        from sparc.io import read_matrix

        mtx_dir, expected_barcodes, expected_genes = sample_mtx_dir
        matrix, barcodes, genes = read_matrix(mtx_dir)

        assert isinstance(matrix, sp.csr_matrix)
        assert matrix.shape[0] == len(expected_barcodes)
        assert len(barcodes) == len(expected_barcodes)
        assert len(genes) == len(expected_genes)

    def test_write_matrix(self, tmp_dir, sample_matrix):
        from sparc.io import read_matrix, write_matrix

        matrix, barcodes, genes = sample_matrix
        out_dir = tmp_dir / "output"

        write_matrix(matrix, barcodes, genes, out_dir)

        assert (out_dir / "matrix.mtx").exists()
        assert (out_dir / "barcodes.tsv").exists()
        assert (out_dir / "genes.tsv").exists()

        # Roundtrip
        matrix2, barcodes2, genes2 = read_matrix(out_dir)
        assert matrix2.shape == matrix.shape
        assert barcodes2 == barcodes

    def test_write_h5ad(self, tmp_dir, sample_matrix):
        pytest.importorskip("anndata")
        from sparc.io import write_h5ad, read_h5ad

        matrix, barcodes, genes = sample_matrix
        h5ad_path = tmp_dir / "test.h5ad"

        write_h5ad(matrix, barcodes, genes, h5ad_path)
        assert h5ad_path.exists()

        matrix2, barcodes2, genes2 = read_h5ad(h5ad_path)
        assert matrix2.shape == matrix.shape
        assert barcodes2 == barcodes
        assert genes2 == genes

    def test_write_matrix_compressed(self, tmp_dir, sample_matrix):
        from sparc.io import write_matrix

        matrix, barcodes, genes = sample_matrix
        out_dir = tmp_dir / "compressed"

        write_matrix(matrix, barcodes, genes, out_dir, compress=True)

        assert (out_dir / "matrix.mtx.gz").exists()
        assert (out_dir / "barcodes.tsv.gz").exists()
        assert (out_dir / "genes.tsv.gz").exists()
        assert not (out_dir / "matrix.mtx").exists()


class TestAnalysis:
    """Tests for sparc.analysis module."""

    def test_to_anndata(self, sample_matrix):
        pytest.importorskip("anndata")
        pytest.importorskip("scanpy")
        from sparc.analysis import to_anndata

        matrix, barcodes, genes = sample_matrix
        adata = to_anndata(matrix, barcodes, genes)

        assert adata.n_obs == len(barcodes)
        assert adata.n_vars == len(genes)

    def test_from_anndata_roundtrip(self, sample_matrix):
        pytest.importorskip("anndata")
        pytest.importorskip("scanpy")
        pytest.importorskip("sparc._sparc_py")
        from sparc.analysis import to_anndata, from_anndata

        matrix, barcodes, genes = sample_matrix
        adata = to_anndata(matrix, barcodes, genes)
        count_matrix = from_anndata(adata)

        assert count_matrix.n_cols == len(barcodes)
        assert count_matrix.n_rows == len(genes)

    def test_normalize_and_analyze(self, sample_matrix):
        pytest.importorskip("anndata")
        sc = pytest.importorskip("scanpy")
        leidenalg = pytest.importorskip("leidenalg")
        from sparc.analysis import to_anndata, normalize_and_analyze

        # Need larger matrix for HVG detection
        np.random.seed(42)
        n_cells, n_genes = 100, 200
        data = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)
        matrix = sp.csr_matrix(data)
        barcodes = [f"CELL{i}" for i in range(n_cells)]
        genes = [f"GENE{i}" for i in range(n_genes)]

        adata = to_anndata(matrix, barcodes, genes)
        result = normalize_and_analyze(adata, n_top_genes=50, n_pcs=10)

        assert "X_pca" in result.obsm
        assert "X_umap" in result.obsm
        assert "leiden" in result.obs.columns

    def test_find_marker_genes(self, sample_matrix):
        pytest.importorskip("anndata")
        sc = pytest.importorskip("scanpy")
        leidenalg = pytest.importorskip("leidenalg")
        from sparc.analysis import to_anndata, normalize_and_analyze, find_marker_genes

        np.random.seed(42)
        n_cells, n_genes = 100, 200
        data = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)
        matrix = sp.csr_matrix(data)
        barcodes = [f"CELL{i}" for i in range(n_cells)]
        genes = [f"GENE{i}" for i in range(n_genes)]

        adata = to_anndata(matrix, barcodes, genes)
        adata = normalize_and_analyze(adata, n_top_genes=50, n_pcs=10)
        adata = find_marker_genes(adata, n_genes=5)

        assert adata.uns.get("rank_genes_groups") is not None


class TestStreaming:
    """Tests for sparc.streaming module."""

    def test_streaming_processor_init(self):
        from sparc.streaming import StreamingProcessor

        processor = StreamingProcessor(chunk_size=50000, max_memory_mb=2048)
        assert processor.chunk_size == 50000
        assert processor.max_memory_mb == 2048

    def test_stream_stats_defaults(self):
        from sparc.streaming import StreamStats

        stats = StreamStats()
        assert stats.total_records == 0
        assert stats.processed_records == 0
        assert stats.chunks_processed == 0


class TestBatch:
    """Tests for sparc.batch module."""

    def test_parse_manifest(self, tmp_dir):
        from sparc.batch import parse_manifest

        manifest_path = tmp_dir / "manifest.csv"
        manifest_path.write_text(
            "sample_name,r1,r2,whitelist\n"
            "sample1,/data/s1_R1.fastq.gz,/data/s1_R2.fastq.gz,/data/whitelist.txt\n"
            "sample2,/data/s2_R1.fastq.gz,/data/s2_R2.fastq.gz,/data/whitelist.txt\n"
        )

        samples = parse_manifest(manifest_path)
        assert len(samples) == 2
        assert samples[0].name == "sample1"
        assert samples[1].name == "sample2"
        assert samples[0].r1 == Path("/data/s1_R1.fastq.gz")

    def test_parse_manifest_with_protocol(self, tmp_dir):
        from sparc.batch import parse_manifest

        manifest_path = tmp_dir / "manifest.csv"
        manifest_path.write_text(
            "sample_name,r1,r2,whitelist,protocol\n"
            "sample1,/data/s1_R1.fq.gz,/data/s1_R2.fq.gz,/data/wl.txt,drop-seq\n"
        )

        samples = parse_manifest(manifest_path)
        assert len(samples) == 1
        assert samples[0].protocol == "drop-seq"

    def test_batch_processor_init(self, tmp_dir):
        from sparc.batch import BatchProcessor

        manifest = tmp_dir / "manifest.csv"
        manifest.write_text("sample_name,r1,r2,whitelist\n")

        processor = BatchProcessor(
            manifest=manifest,
            reference="/ref",
            output=tmp_dir / "output",
            parallel=4,
        )
        assert processor.parallel == 4
