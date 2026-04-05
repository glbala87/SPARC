"""
Background pipeline worker using Celery.
"""

import json
import logging
import os
from collections import defaultdict
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Celery configuration
REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")

try:
    from celery import Celery

    celery_app = Celery(
        "sparc",
        broker=REDIS_URL,
        backend=REDIS_URL,
    )

    celery_app.conf.update(
        task_serializer="json",
        accept_content=["json"],
        result_serializer="json",
        timezone="UTC",
        enable_utc=True,
        task_track_started=True,
        task_time_limit=3600,  # 1 hour
        worker_prefetch_multiplier=1,
    )

    CELERY_AVAILABLE = True
except ImportError:
    CELERY_AVAILABLE = False
    celery_app = None


# Protocol configurations
PROTOCOL_CONFIG = {
    "10x-3prime-v3": {"barcode_start": 0, "barcode_len": 16, "umi_start": 16, "umi_len": 12},
    "10x-3prime-v2": {"barcode_start": 0, "barcode_len": 16, "umi_start": 16, "umi_len": 10},
    "10x-5prime-v2": {"barcode_start": 0, "barcode_len": 16, "umi_start": 16, "umi_len": 10},
    "drop-seq": {"barcode_start": 0, "barcode_len": 12, "umi_start": 12, "umi_len": 8},
    "indrop": {"barcode_start": 0, "barcode_len": 16, "umi_start": 16, "umi_len": 6},
    "sci-rna-seq": {"barcode_start": 0, "barcode_len": 10, "umi_start": 10, "umi_len": 8},
}


def run_pipeline(
    job_id: str,
    r1_path: str,
    r2_path: Optional[str],
    whitelist_path: Optional[str],
    output_dir: str,
    config: dict,
) -> dict:
    """
    Run the analysis pipeline.

    This is a synchronous function that can be called directly
    or wrapped in a Celery task.
    """
    try:
        result = {
            "total_reads": 0,
            "valid_barcodes": 0,
            "corrected_barcodes": 0,
            "invalid_barcodes": 0,
            "cells": 0,
            "genes": 0,
            "median_genes_per_cell": 0,
            "median_umis_per_cell": 0,
        }

        protocol = config.get("protocol", "10x-3prime-v3")
        proto_cfg = PROTOCOL_CONFIG.get(protocol, PROTOCOL_CONFIG["10x-3prime-v3"])

        barcode_start = proto_cfg["barcode_start"]
        barcode_len = proto_cfg["barcode_len"]
        umi_start = proto_cfg["umi_start"]
        umi_len = proto_cfg["umi_len"]

        # Try to use Rust bindings for performance
        try:
            from sparc import Whitelist, BarcodeCorrector, FastqParser, GeneCounter
            return _run_with_rust(
                job_id, r1_path, r2_path, whitelist_path, output_dir,
                config, result, barcode_start, barcode_len, umi_start, umi_len,
            )
        except ImportError:
            logger.info("Rust bindings not available, using pure Python pipeline")

        # Pure Python fallback
        return _run_pure_python(
            job_id, r1_path, r2_path, whitelist_path, output_dir,
            config, result, barcode_start, barcode_len, umi_start, umi_len,
        )

    except Exception as e:
        logger.exception("Pipeline failed for job %s", job_id)
        return {"error": str(e)}


def _run_with_rust(
    job_id, r1_path, r2_path, whitelist_path, output_dir,
    config, result, barcode_start, barcode_len, umi_start, umi_len,
):
    """Run pipeline using Rust bindings."""
    from sparc import Whitelist, BarcodeCorrector, FastqParser, GeneCounter

    max_mismatch = config.get("max_mismatch", 1)

    # Load whitelist
    corrector = None
    if whitelist_path and Path(whitelist_path).exists():
        whitelist = Whitelist(whitelist_path)
        corrector = BarcodeCorrector(whitelist, max_mismatch)

    # Step 1: Extract barcodes from R1 and build count matrix
    parser = FastqParser(r1_path)
    counter = GeneCounter()

    # Track barcode-UMI-gene associations
    barcode_umi_gene = defaultdict(lambda: defaultdict(set))

    for record in parser:
        result["total_reads"] += 1

        if len(record.seq) < barcode_start + barcode_len + umi_len:
            continue

        barcode_seq = record.subsequence(barcode_start, barcode_len)
        umi_seq = record.subsequence(umi_start, umi_len)

        if barcode_seq is None or umi_seq is None:
            continue

        barcode_str = bytes(barcode_seq).decode("ascii", errors="ignore")
        umi_str = bytes(umi_seq).decode("ascii", errors="ignore")

        # Match barcode
        final_barcode = None
        if corrector:
            status, corrected, _ = corrector.match_barcode(barcode_str)
            if status == "exact":
                result["valid_barcodes"] += 1
                final_barcode = corrected
            elif status == "corrected":
                result["valid_barcodes"] += 1
                result["corrected_barcodes"] += 1
                final_barcode = corrected
            else:
                result["invalid_barcodes"] += 1
        else:
            final_barcode = barcode_str
            result["valid_barcodes"] += 1

        if final_barcode:
            # Extract gene info from read name if available
            gene = _extract_gene_from_read(record.id)
            if gene:
                barcode_umi_gene[final_barcode][gene].add(umi_str)

    # Step 2: UMI deduplication and count matrix building
    for barcode, gene_umis in barcode_umi_gene.items():
        for gene, umis in gene_umis.items():
            # Count unique UMIs per barcode-gene pair
            counter.add_count(barcode, gene, len(umis))

    # Step 3: Process R2 if available (gene assignment from alignment)
    if r2_path and Path(r2_path).exists():
        logger.info("R2 file available at %s", r2_path)
        # R2 processing would typically require alignment; counts already built from R1 tags

    # Step 4: Build and write matrix
    matrix = counter.build()
    result["cells"] = matrix.n_cols
    result["genes"] = matrix.n_rows

    # Compute QC stats
    if matrix.n_cols > 0:
        import numpy as np
        genes_per_cell = np.array(matrix.genes_per_cell())
        counts_per_cell = np.array(matrix.counts_per_cell())
        result["median_genes_per_cell"] = int(np.median(genes_per_cell)) if len(genes_per_cell) > 0 else 0
        result["median_umis_per_cell"] = int(np.median(counts_per_cell)) if len(counts_per_cell) > 0 else 0

    # Save output
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    matrix.write_mtx(str(output_path / "matrix.mtx"))
    matrix.write_barcodes(str(output_path / "barcodes.tsv"))
    matrix.write_genes(str(output_path / "genes.tsv"))

    # Write QC report
    qc_report = {
        "job_id": job_id,
        "total_reads": result["total_reads"],
        "valid_barcodes": result["valid_barcodes"],
        "corrected_barcodes": result["corrected_barcodes"],
        "invalid_barcodes": result["invalid_barcodes"],
        "cells": result["cells"],
        "genes": result["genes"],
        "median_genes_per_cell": result["median_genes_per_cell"],
        "median_umis_per_cell": result["median_umis_per_cell"],
    }
    with open(output_path / "qc_report.json", "w") as f:
        json.dump(qc_report, f, indent=2)

    return result


def _run_pure_python(
    job_id, r1_path, r2_path, whitelist_path, output_dir,
    config, result, barcode_start, barcode_len, umi_start, umi_len,
):
    """Pure Python fallback pipeline when Rust bindings are unavailable."""
    import gzip
    from scipy.io import mmwrite
    import scipy.sparse as sp
    import numpy as np

    max_mismatch = config.get("max_mismatch", 1)

    # Load whitelist
    whitelist_set = set()
    if whitelist_path and Path(whitelist_path).exists():
        with open(whitelist_path) as f:
            for line in f:
                bc = line.strip()
                if bc and not bc.startswith("#"):
                    whitelist_set.add(bc)

    # Parse FASTQ
    barcode_gene_counts = defaultdict(lambda: defaultdict(set))

    open_fn = gzip.open if r1_path.endswith(".gz") else open
    mode = "rt" if r1_path.endswith(".gz") else "r"

    with open_fn(r1_path, mode) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # +
            f.readline()  # qual

            result["total_reads"] += 1

            if len(seq) < barcode_start + barcode_len + umi_len:
                continue

            barcode = seq[barcode_start:barcode_start + barcode_len]
            umi = seq[umi_start:umi_start + umi_len]

            if whitelist_set and barcode not in whitelist_set:
                result["invalid_barcodes"] += 1
                continue

            result["valid_barcodes"] += 1

            gene = _extract_gene_from_read(header)
            if gene:
                barcode_gene_counts[barcode][gene].add(umi)

    # Build count matrix
    all_barcodes = sorted(barcode_gene_counts.keys())
    all_genes = sorted({g for bc_genes in barcode_gene_counts.values() for g in bc_genes})

    bc_idx = {bc: i for i, bc in enumerate(all_barcodes)}
    gene_idx = {g: i for i, g in enumerate(all_genes)}

    rows, cols, vals = [], [], []
    for bc, gene_umis in barcode_gene_counts.items():
        for gene, umis in gene_umis.items():
            rows.append(gene_idx[gene])
            cols.append(bc_idx[bc])
            vals.append(len(umis))

    result["cells"] = len(all_barcodes)
    result["genes"] = len(all_genes)

    # Save
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if rows:
        matrix = sp.coo_matrix(
            (np.array(vals, dtype=np.int32), (np.array(rows), np.array(cols))),
            shape=(len(all_genes), len(all_barcodes)),
        )
        mmwrite(str(output_path / "matrix.mtx"), matrix)

        # Compute QC
        genes_per_cell = np.diff(matrix.tocsc().indptr)
        counts_per_cell = np.array(matrix.tocsc().sum(axis=0)).ravel()
        result["median_genes_per_cell"] = int(np.median(genes_per_cell)) if len(genes_per_cell) > 0 else 0
        result["median_umis_per_cell"] = int(np.median(counts_per_cell)) if len(counts_per_cell) > 0 else 0

    with open(output_path / "barcodes.tsv", "w") as f:
        for bc in all_barcodes:
            f.write(f"{bc}\n")

    with open(output_path / "genes.tsv", "w") as f:
        for gene in all_genes:
            f.write(f"{gene}\t{gene}\n")

    # QC report
    qc_report = {
        "job_id": job_id,
        **{k: v for k, v in result.items()},
    }
    with open(output_path / "qc_report.json", "w") as f:
        json.dump(qc_report, f, indent=2)

    return result


def _extract_gene_from_read(read_id: str) -> Optional[str]:
    """Extract gene name from read ID if tagged (e.g., READ_00000001:BARCODE:GENE)."""
    parts = read_id.split(":")
    if len(parts) >= 3:
        gene = parts[-1].strip()
        if gene and gene != "NONE":
            return gene
    return None


if CELERY_AVAILABLE:

    @celery_app.task(bind=True)
    def run_pipeline_task(
        self,
        job_id: str,
        r1_path: str,
        r2_path: Optional[str],
        whitelist_path: Optional[str],
        output_dir: str,
        config: dict,
    ):
        """Celery task wrapper for pipeline."""

        def update_progress(progress: float, message: str):
            self.update_state(
                state="PROGRESS",
                meta={"progress": progress, "message": message},
            )

        update_progress(0.1, "Starting pipeline...")

        result = run_pipeline(
            job_id=job_id,
            r1_path=r1_path,
            r2_path=r2_path,
            whitelist_path=whitelist_path,
            output_dir=output_dir,
            config=config,
        )

        if "error" not in result:
            update_progress(1.0, "Pipeline completed")

        return result
