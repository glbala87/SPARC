"""
Background pipeline worker using Celery.
"""

import os
from pathlib import Path
from typing import Optional

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
        from sparc import Whitelist, BarcodeCorrector, FastqParser, GeneCounter

        result = {
            "total_reads": 0,
            "valid_barcodes": 0,
            "corrected_barcodes": 0,
            "cells": 0,
            "genes": 0,
        }

        # Load whitelist if provided
        if whitelist_path:
            whitelist = Whitelist(whitelist_path)
            corrector = BarcodeCorrector(whitelist, config.get("max_mismatch", 1))
        else:
            corrector = None

        # Process FASTQ
        parser = FastqParser(r1_path)

        barcode_start = 0
        barcode_len = 16
        umi_start = 16
        umi_len = 12

        if config.get("protocol") == "10x-3prime-v2":
            umi_len = 10
        elif config.get("protocol") == "10x-5prime-v2":
            umi_len = 10

        counter = GeneCounter()

        for record in parser:
            result["total_reads"] += 1

            if len(record.seq) < barcode_start + barcode_len + umi_len:
                continue

            barcode_seq = record.subsequence(barcode_start, barcode_len)
            if barcode_seq is None:
                continue

            barcode_str = bytes(barcode_seq).decode("ascii", errors="ignore")

            if corrector:
                status, corrected, _ = corrector.match_barcode(barcode_str)
                if status == "exact":
                    result["valid_barcodes"] += 1
                elif status == "corrected":
                    result["valid_barcodes"] += 1
                    result["corrected_barcodes"] += 1

        # Build matrix (placeholder)
        matrix = counter.build()
        result["cells"] = matrix.n_cols
        result["genes"] = matrix.n_rows

        # Save output
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        matrix.write_mtx(str(output_path / "matrix.mtx"))
        matrix.write_barcodes(str(output_path / "barcodes.tsv"))
        matrix.write_genes(str(output_path / "genes.tsv"))

        return result

    except Exception as e:
        return {"error": str(e)}


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

        return result
