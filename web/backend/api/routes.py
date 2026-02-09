"""
API routes for sparc web interface.
"""

import os
import uuid
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, File, UploadFile, HTTPException, BackgroundTasks
from pydantic import BaseModel

router = APIRouter()

# Configuration
UPLOAD_DIR = Path(os.getenv("SCTOOLS_UPLOAD_DIR", "/tmp/sparc/uploads"))
OUTPUT_DIR = Path(os.getenv("SCTOOLS_OUTPUT_DIR", "/tmp/sparc/outputs"))

# In-memory job store (use Redis in production)
jobs: dict[str, dict] = {}


class PipelineConfig(BaseModel):
    """Pipeline configuration."""
    sample_name: str = "sample"
    protocol: str = "10x-3prime-v3"
    max_mismatch: int = 1
    min_genes: int = 200
    max_genes: int = 10000
    max_mito: float = 20.0
    n_pcs: int = 50
    resolution: float = 1.0


class JobStatus(BaseModel):
    """Job status response."""
    job_id: str
    status: str
    progress: float
    message: str
    result: Optional[dict] = None


@router.post("/upload")
async def upload_files(
    r1: UploadFile = File(...),
    r2: Optional[UploadFile] = File(None),
    whitelist: Optional[UploadFile] = File(None),
):
    """Upload FASTQ and whitelist files."""
    job_id = str(uuid.uuid4())
    job_dir = UPLOAD_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    # Save R1
    r1_path = job_dir / r1.filename
    with open(r1_path, "wb") as f:
        content = await r1.read()
        f.write(content)

    # Save R2 if provided
    r2_path = None
    if r2:
        r2_path = job_dir / r2.filename
        with open(r2_path, "wb") as f:
            content = await r2.read()
            f.write(content)

    # Save whitelist if provided
    whitelist_path = None
    if whitelist:
        whitelist_path = job_dir / whitelist.filename
        with open(whitelist_path, "wb") as f:
            content = await whitelist.read()
            f.write(content)

    return {
        "job_id": job_id,
        "files": {
            "r1": str(r1_path),
            "r2": str(r2_path) if r2_path else None,
            "whitelist": str(whitelist_path) if whitelist_path else None,
        },
    }


@router.post("/pipeline/{job_id}")
async def start_pipeline(
    job_id: str,
    config: PipelineConfig,
    background_tasks: BackgroundTasks,
):
    """Start the analysis pipeline."""
    job_dir = UPLOAD_DIR / job_id

    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found")

    # Initialize job status
    jobs[job_id] = {
        "status": "queued",
        "progress": 0.0,
        "message": "Pipeline queued",
        "config": config.model_dump(),
    }

    # Start background task
    background_tasks.add_task(run_pipeline_task, job_id, config)

    return {"job_id": job_id, "status": "queued"}


async def run_pipeline_task(job_id: str, config: PipelineConfig):
    """Run the pipeline as a background task."""
    try:
        jobs[job_id]["status"] = "running"
        jobs[job_id]["progress"] = 0.1
        jobs[job_id]["message"] = "Starting pipeline..."

        # Simulate pipeline steps
        import asyncio

        # Step 1: Extract barcodes
        jobs[job_id]["progress"] = 0.2
        jobs[job_id]["message"] = "Extracting barcodes and UMIs..."
        await asyncio.sleep(1)

        # Step 2: Align reads
        jobs[job_id]["progress"] = 0.4
        jobs[job_id]["message"] = "Aligning reads..."
        await asyncio.sleep(1)

        # Step 3: Count matrix
        jobs[job_id]["progress"] = 0.6
        jobs[job_id]["message"] = "Generating count matrix..."
        await asyncio.sleep(1)

        # Step 4: QC
        jobs[job_id]["progress"] = 0.8
        jobs[job_id]["message"] = "Running QC..."
        await asyncio.sleep(1)

        # Complete
        jobs[job_id]["status"] = "completed"
        jobs[job_id]["progress"] = 1.0
        jobs[job_id]["message"] = "Pipeline completed successfully"
        jobs[job_id]["result"] = {
            "total_reads": 1000000,
            "valid_barcodes": 950000,
            "cells": 5000,
            "genes": 20000,
            "median_genes_per_cell": 2500,
        }

    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["message"] = str(e)


@router.get("/pipeline/{job_id}/status")
async def get_pipeline_status(job_id: str) -> JobStatus:
    """Get pipeline status."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]
    return JobStatus(
        job_id=job_id,
        status=job["status"],
        progress=job["progress"],
        message=job["message"],
        result=job.get("result"),
    )


@router.get("/pipeline/{job_id}/results")
async def get_pipeline_results(job_id: str):
    """Get pipeline results."""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs[job_id]

    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Pipeline not completed")

    output_dir = OUTPUT_DIR / job_id

    return {
        "job_id": job_id,
        "result": job.get("result"),
        "files": {
            "matrix": str(output_dir / "matrix.mtx") if output_dir.exists() else None,
            "barcodes": str(output_dir / "barcodes.tsv") if output_dir.exists() else None,
            "genes": str(output_dir / "genes.tsv") if output_dir.exists() else None,
            "qc_report": str(output_dir / "qc_report.json") if output_dir.exists() else None,
        },
    }


@router.get("/whitelists")
async def list_whitelists():
    """List available barcode whitelists."""
    return {
        "whitelists": [
            {"name": "10x 3M v3", "file": "3M-february-2018.txt", "barcodes": 6794880},
            {"name": "10x 737K v2", "file": "737K-august-2016.txt", "barcodes": 737280},
        ]
    }


@router.get("/protocols")
async def list_protocols():
    """List available protocols."""
    return {
        "protocols": [
            {
                "id": "10x-3prime-v3",
                "name": "10x Genomics 3' v3",
                "barcode_len": 16,
                "umi_len": 12,
            },
            {
                "id": "10x-3prime-v2",
                "name": "10x Genomics 3' v2",
                "barcode_len": 16,
                "umi_len": 10,
            },
            {
                "id": "10x-5prime-v2",
                "name": "10x Genomics 5' v2",
                "barcode_len": 16,
                "umi_len": 10,
            },
        ]
    }
