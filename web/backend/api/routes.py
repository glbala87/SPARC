"""
API routes for sparc web interface.
"""

import asyncio
import json
import logging
import os
import re
import shutil
import time
import uuid
from pathlib import Path, PurePosixPath
from typing import Optional

from fastapi import APIRouter, File, UploadFile, HTTPException, BackgroundTasks, Depends, Header, Request
from pydantic import BaseModel, field_validator

logger = logging.getLogger("sparc.api")

router = APIRouter()

# ─── Authentication ───────────────────────────────────────────────────

_raw_tokens = os.getenv("SPARC_API_TOKENS", "")
if not _raw_tokens:
    logger.warning(
        "SPARC_API_TOKENS not set — all authenticated endpoints will reject requests. "
        "Set SPARC_API_TOKENS=<token> or SPARC_API_TOKENS= (empty) to disable auth."
    )
    API_TOKENS: set[str] = set()
    _AUTH_DISABLED = False
elif _raw_tokens.strip() == "":
    API_TOKENS = set()
    _AUTH_DISABLED = True
    logger.warning("SPARC_API_TOKENS is empty — authentication is DISABLED")
else:
    API_TOKENS = set(t.strip() for t in _raw_tokens.split(",") if t.strip())
    _AUTH_DISABLED = False


def verify_token(authorization: str = Header(None)):
    """Verify API token from Authorization header."""
    if _AUTH_DISABLED:
        return
    if not authorization:
        raise HTTPException(status_code=401, detail="Missing Authorization header")
    token = authorization.replace("Bearer ", "").strip()
    if token not in API_TOKENS:
        raise HTTPException(status_code=403, detail="Invalid API token")


# ─── Configuration ────────────────────────────────────────────────────

UPLOAD_DIR = Path(os.getenv("SPARC_UPLOAD_DIR", os.getenv("SCTOOLS_UPLOAD_DIR", "/tmp/sparc/uploads")))
OUTPUT_DIR = Path(os.getenv("SPARC_OUTPUT_DIR", os.getenv("SCTOOLS_OUTPUT_DIR", "/tmp/sparc/outputs")))
MAX_UPLOAD_SIZE = int(os.getenv("SPARC_MAX_UPLOAD_MB", "5000")) * 1024 * 1024
MAX_WHITELIST_SIZE = int(os.getenv("SPARC_MAX_WHITELIST_MB", "500")) * 1024 * 1024
MAX_JOBS = int(os.getenv("SPARC_MAX_JOBS", "100"))
MAX_CONCURRENT_PIPELINES = int(os.getenv("SPARC_MAX_CONCURRENT", "4"))

# Track running pipelines
_running_pipelines = 0
_pipeline_lock = asyncio.Lock()

# ─── Job Store (thread-safe) ─────────────────────────────────────────

try:
    import redis as redis_lib
    _REDIS_URL = os.getenv("REDIS_URL", "")
    if _REDIS_URL:
        _redis = redis_lib.Redis.from_url(_REDIS_URL, decode_responses=True)
        _redis.ping()
        _USE_REDIS = True
        logger.info("Using Redis for job storage at %s", _REDIS_URL.split("@")[-1])
    else:
        _USE_REDIS = False
except Exception:
    _USE_REDIS = False

_jobs_mem: dict[str, dict] = {}
_jobs_lock = asyncio.Lock()


async def _get_job(job_id: str) -> Optional[dict]:
    if _USE_REDIS:
        data = _redis.get(f"sparc:job:{job_id}")
        return json.loads(data) if data else None
    async with _jobs_lock:
        return _jobs_mem.get(job_id)


async def _set_job(job_id: str, job: dict):
    if _USE_REDIS:
        _redis.set(f"sparc:job:{job_id}", json.dumps(job), ex=86400)
        return
    async with _jobs_lock:
        if job_id not in _jobs_mem and len(_jobs_mem) >= MAX_JOBS:
            oldest = next(iter(_jobs_mem))
            logger.warning("Job store full, evicting oldest job %s", oldest)
            del _jobs_mem[oldest]
        _jobs_mem[job_id] = job


async def _delete_job(job_id: str):
    if _USE_REDIS:
        _redis.delete(f"sparc:job:{job_id}")
        return
    async with _jobs_lock:
        _jobs_mem.pop(job_id, None)


async def _list_jobs() -> list[dict]:
    if _USE_REDIS:
        keys = _redis.keys("sparc:job:*")
        jobs = []
        for key in keys:
            data = _redis.get(key)
            if data:
                jobs.append(json.loads(data))
        return jobs
    async with _jobs_lock:
        return list(_jobs_mem.values())


# ─── Security helpers ─────────────────────────────────────────────────

_SAFE_FILENAME_RE = re.compile(r"^[a-zA-Z0-9][a-zA-Z0-9._\-]{0,254}$")


def _sanitize_filename(filename: str) -> str:
    """Sanitize an upload filename to prevent path traversal."""
    # Take only the basename (strip any directory components)
    name = PurePosixPath(filename).name
    # Remove any remaining path separators
    name = name.replace("/", "_").replace("\\", "_").replace("\x00", "")
    if not name or name.startswith("."):
        name = f"upload_{uuid.uuid4().hex[:8]}"
    if not _SAFE_FILENAME_RE.match(name):
        ext = Path(name).suffix if "." in name else ""
        name = f"upload_{uuid.uuid4().hex[:8]}{ext}"
    return name


def _relative_path(path: Path, base: Path) -> str:
    """Return a relative path string, never exposing absolute filesystem paths."""
    try:
        return str(path.relative_to(base))
    except ValueError:
        return path.name


# ─── Validation models ────────────────────────────────────────────────

VALID_PROTOCOLS = {
    "10x-3prime-v3", "10x-3prime-v2", "10x-5prime-v2",
    "drop-seq", "indrop", "sci-rna-seq", "smart-seq2",
}


class PipelineConfig(BaseModel):
    """Pipeline configuration with validated constraints."""
    sample_name: str = "sample"
    protocol: str = "10x-3prime-v3"
    max_mismatch: int = 1
    min_genes: int = 200
    max_genes: int = 10000
    max_mito: float = 20.0
    n_pcs: int = 50
    resolution: float = 1.0

    @field_validator("protocol")
    @classmethod
    def validate_protocol(cls, v):
        if v not in VALID_PROTOCOLS:
            raise ValueError(f"Invalid protocol: {v}. Must be one of: {VALID_PROTOCOLS}")
        return v

    @field_validator("max_mismatch")
    @classmethod
    def validate_max_mismatch(cls, v):
        if not 0 <= v <= 3:
            raise ValueError("max_mismatch must be 0-3")
        return v

    @field_validator("min_genes")
    @classmethod
    def validate_min_genes(cls, v):
        if not 0 <= v <= 100000:
            raise ValueError("min_genes must be 0-100000")
        return v

    @field_validator("max_genes")
    @classmethod
    def validate_max_genes(cls, v):
        if not 1 <= v <= 1000000:
            raise ValueError("max_genes must be 1-1000000")
        return v

    @field_validator("max_mito")
    @classmethod
    def validate_max_mito(cls, v):
        if not 0.0 <= v <= 100.0:
            raise ValueError("max_mito must be 0-100")
        return v

    @field_validator("n_pcs")
    @classmethod
    def validate_n_pcs(cls, v):
        if not 2 <= v <= 500:
            raise ValueError("n_pcs must be 2-500")
        return v

    @field_validator("resolution")
    @classmethod
    def validate_resolution(cls, v):
        if not 0.01 <= v <= 50.0:
            raise ValueError("resolution must be 0.01-50.0")
        return v


class JobStatus(BaseModel):
    """Job status response."""
    job_id: str
    status: str
    progress: float
    message: str
    result: Optional[dict] = None


# ─── FASTQ validation ─────────────────────────────────────────────────

def _validate_fastq(path: Path) -> bool:
    """Validate that file looks like a FASTQ (checks header + first record structure)."""
    try:
        with open(path, "rb") as f:
            header = f.read(4)
            if not header:
                return False
            # Gzip-compressed
            if header[:2] == b"\x1f\x8b":
                return True
            # Zstd-compressed
            if header[:4] == b"\x28\xb5\x2f\xfd":
                return True
            # Plain FASTQ: first char must be '@'
            if header[0:1] != b"@":
                return False
        # Read first 4 lines to validate structure
        with open(path, "r") as f:
            line1 = f.readline()  # @header
            line2 = f.readline()  # sequence
            line3 = f.readline()  # +
            line4 = f.readline()  # quality
            if not (line1.startswith("@") and line3.startswith("+") and len(line2.strip()) == len(line4.strip())):
                return False
        return True
    except Exception:
        return False


# ─── Endpoints ────────────────────────────────────────────────────────

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

    try:
        # Save R1
        r1_name = _sanitize_filename(r1.filename or "r1.fastq")
        r1_path = job_dir / r1_name
        content = await r1.read()
        if len(content) > MAX_UPLOAD_SIZE:
            raise HTTPException(status_code=413, detail="R1 file exceeds max upload size")
        r1_path.write_bytes(content)
        logger.info("Job %s: uploaded R1 (%s, %d bytes)", job_id, r1_name, len(content))

        if not _validate_fastq(r1_path):
            raise HTTPException(status_code=400, detail="R1 does not appear to be a valid FASTQ file")

        # Save R2 if provided
        r2_name = None
        if r2:
            r2_name = _sanitize_filename(r2.filename or "r2.fastq")
            r2_path = job_dir / r2_name
            content = await r2.read()
            if len(content) > MAX_UPLOAD_SIZE:
                raise HTTPException(status_code=413, detail="R2 file exceeds max upload size")
            r2_path.write_bytes(content)
            logger.info("Job %s: uploaded R2 (%s, %d bytes)", job_id, r2_name, len(content))

            if not _validate_fastq(r2_path):
                raise HTTPException(status_code=400, detail="R2 does not appear to be a valid FASTQ file")

        # Save whitelist if provided
        wl_name = None
        if whitelist:
            wl_name = _sanitize_filename(whitelist.filename or "whitelist.txt")
            wl_path = job_dir / wl_name
            content = await whitelist.read()
            if len(content) > MAX_WHITELIST_SIZE:
                raise HTTPException(status_code=413, detail="Whitelist file exceeds max size")
            wl_path.write_bytes(content)
            logger.info("Job %s: uploaded whitelist (%s, %d bytes)", job_id, wl_name, len(content))

    except HTTPException:
        shutil.rmtree(job_dir, ignore_errors=True)
        raise
    except Exception as e:
        shutil.rmtree(job_dir, ignore_errors=True)
        logger.exception("Upload failed for job %s", job_id)
        raise HTTPException(status_code=500, detail="Upload failed")

    return {
        "job_id": job_id,
        "files": {
            "r1": r1_name,
            "r2": r2_name,
            "whitelist": wl_name,
        },
    }


@router.post("/pipeline/{job_id}")
async def start_pipeline(
    job_id: str,
    config: PipelineConfig,
    background_tasks: BackgroundTasks,
):
    """Start the analysis pipeline."""
    # Validate job_id format (UUID)
    try:
        uuid.UUID(job_id)
    except ValueError:
        raise HTTPException(status_code=400, detail="Invalid job ID format")

    job_dir = UPLOAD_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found. Upload files first.")

    r1_files = list(job_dir.glob("*.fastq*")) + list(job_dir.glob("*.fq*"))
    if not r1_files:
        raise HTTPException(status_code=400, detail="No FASTQ files found in upload directory")

    existing = await _get_job(job_id)
    if existing and existing.get("status") == "running":
        raise HTTPException(status_code=409, detail="Pipeline already running for this job")

    # Check concurrent pipeline limit
    global _running_pipelines
    async with _pipeline_lock:
        if _running_pipelines >= MAX_CONCURRENT_PIPELINES:
            raise HTTPException(
                status_code=429,
                detail=f"Too many concurrent pipelines ({MAX_CONCURRENT_PIPELINES} max). Try again later.",
            )

    job = {
        "job_id": job_id,
        "status": "queued",
        "progress": 0.0,
        "message": "Pipeline queued",
        "config": config.model_dump(),
        "result": None,
        "started_at": time.time(),
    }
    await _set_job(job_id, job)
    logger.info("Job %s: pipeline queued (protocol=%s)", job_id, config.protocol)

    background_tasks.add_task(run_pipeline_task, job_id, config)

    return {"job_id": job_id, "status": "queued"}


async def run_pipeline_task(job_id: str, config: PipelineConfig):
    """Run the actual SPARC pipeline as a background task."""
    global _running_pipelines
    async with _pipeline_lock:
        _running_pipelines += 1
    logger.info("Job %s: pipeline started (concurrent=%d)", job_id, _running_pipelines)

    try:
        job = await _get_job(job_id) or {}
        job.update({"status": "running", "progress": 0.05, "message": "Starting pipeline..."})
        await _set_job(job_id, job)

        job_dir = UPLOAD_DIR / job_id
        output_dir = OUTPUT_DIR / job_id
        output_dir.mkdir(parents=True, exist_ok=True)

        r1_path = None
        r2_path = None
        whitelist_path = None

        for f in job_dir.iterdir():
            name = f.name.lower()
            if "whitelist" in name or name.endswith(".txt"):
                whitelist_path = str(f)
            elif "r2" in name or "_2." in name or "_R2" in f.name:
                r2_path = str(f)
            elif r1_path is None:
                r1_path = str(f)

        if not r1_path:
            raise FileNotFoundError("No R1 FASTQ file found")

        from web.backend.workers.pipeline import run_pipeline

        job.update({"progress": 0.1, "message": "Extracting barcodes and UMIs..."})
        await _set_job(job_id, job)

        result = await asyncio.get_event_loop().run_in_executor(
            None, run_pipeline,
            job_id, r1_path, r2_path, whitelist_path, str(output_dir), config.model_dump(),
        )

        if "error" in result:
            logger.error("Job %s: pipeline error: %s", job_id, result["error"])
            job.update({"status": "failed", "message": f"Pipeline error: {result['error']}"})
            await _set_job(job_id, job)
            return

        elapsed = time.time() - job.get("started_at", time.time())
        logger.info("Job %s: pipeline completed in %.1fs (cells=%s, genes=%s)",
                     job_id, elapsed, result.get("cells"), result.get("genes"))
        job.update({
            "status": "completed", "progress": 1.0,
            "message": "Pipeline completed successfully", "result": result,
        })
        await _set_job(job_id, job)

    except Exception as e:
        logger.exception("Job %s: pipeline failed", job_id)
        job = await _get_job(job_id) or {"job_id": job_id}
        job.update({"status": "failed", "message": str(e)})
        await _set_job(job_id, job)
    finally:
        async with _pipeline_lock:
            _running_pipelines -= 1
        logger.info("Job %s: pipeline slot released (concurrent=%d)", job_id, _running_pipelines)


@router.get("/pipeline/{job_id}/status")
async def get_pipeline_status(job_id: str) -> JobStatus:
    """Get pipeline status."""
    job = await _get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return JobStatus(
        job_id=job_id, status=job["status"], progress=job["progress"],
        message=job["message"], result=job.get("result"),
    )


@router.get("/pipeline/{job_id}/results")
async def get_pipeline_results(job_id: str):
    """Get pipeline results (no absolute paths exposed)."""
    job = await _get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Pipeline not completed")

    output_dir = OUTPUT_DIR / job_id
    files = {}
    for name in ["matrix.mtx", "barcodes.tsv", "genes.tsv", "qc_report.json"]:
        files[name.split(".")[0]] = name if (output_dir / name).exists() else None

    return {"job_id": job_id, "result": job.get("result"), "files": files}


@router.delete("/pipeline/{job_id}", dependencies=[Depends(verify_token)])
async def delete_job(job_id: str):
    """Delete a job and its data."""
    job = await _get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.get("status") == "running":
        raise HTTPException(status_code=409, detail="Cannot delete a running job")

    shutil.rmtree(UPLOAD_DIR / job_id, ignore_errors=True)
    shutil.rmtree(OUTPUT_DIR / job_id, ignore_errors=True)
    await _delete_job(job_id)
    logger.info("Job %s: deleted", job_id)
    return {"job_id": job_id, "deleted": True}


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
            {"id": p, "name": p} for p in sorted(VALID_PROTOCOLS)
        ]
    }


@router.post("/whitelists/upload", dependencies=[Depends(verify_token)])
async def upload_custom_whitelist(
    whitelist: UploadFile = File(...),
    name: str = "custom",
):
    """Upload a custom barcode whitelist."""
    # Sanitize name
    safe_name = re.sub(r"[^a-zA-Z0-9_\-]", "_", name)[:64]
    whitelist_dir = UPLOAD_DIR / "whitelists"
    whitelist_dir.mkdir(parents=True, exist_ok=True)

    whitelist_path = whitelist_dir / f"{safe_name}.txt"
    content = await whitelist.read()
    if len(content) > MAX_WHITELIST_SIZE:
        raise HTTPException(status_code=413, detail="Whitelist file exceeds max size")
    whitelist_path.write_bytes(content)
    line_count = content.decode(errors="replace").count("\n")
    logger.info("Uploaded custom whitelist '%s' (%d barcodes)", safe_name, line_count)

    return {"name": safe_name, "barcodes": line_count}


@router.get("/jobs")
async def list_jobs():
    """List all pipeline jobs with their status."""
    all_jobs = await _list_jobs()
    return {
        "jobs": [
            {
                "job_id": job.get("job_id", "unknown"),
                "status": job.get("status", "unknown"),
                "progress": job.get("progress", 0),
                "message": job.get("message", ""),
            }
            for job in all_jobs
        ]
    }


@router.post("/compare", dependencies=[Depends(verify_token)])
async def compare_samples(job_ids: list[str]):
    """Compare results across multiple completed pipeline jobs."""
    results = []
    for job_id in job_ids:
        job = await _get_job(job_id)
        if not job or job["status"] != "completed" or not job.get("result"):
            continue
        results.append({
            "job_id": job_id,
            "sample_name": job.get("config", {}).get("sample_name", job_id),
            **job["result"],
        })
    if not results:
        raise HTTPException(status_code=400, detail="No completed jobs found for comparison")
    return {"samples": results}


@router.get("/pipeline/{job_id}/notebook")
async def export_notebook(job_id: str):
    """Export pipeline results as a Jupyter notebook."""
    job = await _get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    sample_name = job.get("config", {}).get("sample_name", "sample")
    notebook = {
        "nbformat": 4, "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        },
        "cells": [
            {"cell_type": "markdown", "metadata": {},
             "source": [f"# SPARC Analysis: {sample_name}\n"]},
            {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
             "source": ["import sparc\nimport scanpy as sc\n"]},
            {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
             "source": [f"adata = sparc.run_pipeline('{job_id}')\n",
                        "adata = sparc.normalize_and_analyze(adata)\n",
                        "sc.pl.umap(adata, color='leiden')\n"]},
        ],
    }
    return notebook
