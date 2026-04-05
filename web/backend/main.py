"""
FastAPI backend for sparc web interface.
"""

import asyncio
import logging
import os
import signal
import time
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from api.routes import router, _USE_REDIS
from api.websocket import websocket_router

# ─── Logging setup ────────────────────────────────────────────────────

LOG_LEVEL = os.getenv("SPARC_LOG_LEVEL", "INFO").upper()
logging.basicConfig(
    level=getattr(logging, LOG_LEVEL, logging.INFO),
    format="%(asctime)s %(levelname)s [%(name)s] %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("sparc")

# ─── Configuration ────────────────────────────────────────────────────

UPLOAD_DIR = Path(os.getenv("SPARC_UPLOAD_DIR", os.getenv("SCTOOLS_UPLOAD_DIR", "/tmp/sparc/uploads")))
OUTPUT_DIR = Path(os.getenv("SPARC_OUTPUT_DIR", os.getenv("SCTOOLS_OUTPUT_DIR", "/tmp/sparc/outputs")))
ALLOWED_ORIGINS = [
    o.strip() for o in os.getenv("SPARC_CORS_ORIGINS", "http://localhost:3000,http://localhost:5173").split(",")
    if o.strip()
]

# ─── Metrics (simple counters) ────────────────────────────────────────

_metrics = {
    "requests_total": 0,
    "requests_failed": 0,
    "jobs_completed": 0,
    "jobs_failed": 0,
    "uptime_start": time.time(),
}

# ─── Graceful shutdown ────────────────────────────────────────────────

_shutdown_event = asyncio.Event()


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize and cleanup application state."""
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    logger.info("SPARC API starting (upload_dir=%s, output_dir=%s)", UPLOAD_DIR, OUTPUT_DIR)

    # Register signal handlers for graceful shutdown
    loop = asyncio.get_event_loop()
    for sig in (signal.SIGTERM, signal.SIGINT):
        loop.add_signal_handler(sig, lambda s=sig: _handle_signal(s))

    yield

    logger.info("SPARC API shutting down — waiting for running pipelines...")
    _shutdown_event.set()
    # Give running pipelines up to 30s to complete
    for _ in range(30):
        from api.routes import _running_pipelines
        if _running_pipelines <= 0:
            break
        logger.info("Waiting for %d running pipeline(s)...", _running_pipelines)
        await asyncio.sleep(1)
    logger.info("SPARC API shutdown complete")


def _handle_signal(sig):
    logger.info("Received signal %s, initiating graceful shutdown", sig.name)
    _shutdown_event.set()


# ─── App creation ─────────────────────────────────────────────────────

app = FastAPI(
    title="SPARC API",
    description="Single-cell RNA-seq analysis API powered by SPARC",
    version="0.1.0",
    lifespan=lifespan,
)

# CORS — configurable origins, restricted methods
app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["GET", "POST", "DELETE"],
    allow_headers=["Authorization", "Content-Type"],
)


# ─── Request logging middleware ───────────────────────────────────────

@app.middleware("http")
async def log_requests(request: Request, call_next):
    """Log every request with timing and status."""
    start = time.time()
    _metrics["requests_total"] += 1

    response = await call_next(request)

    elapsed = (time.time() - start) * 1000
    status = response.status_code

    if status >= 400:
        _metrics["requests_failed"] += 1

    logger.info(
        "%s %s %d (%.1fms)",
        request.method, request.url.path, status, elapsed,
    )
    return response


# ─── Request size limit middleware ────────────────────────────────────

MAX_REQUEST_BODY = int(os.getenv("SPARC_MAX_REQUEST_MB", "6000")) * 1024 * 1024  # 6GB default


@app.middleware("http")
async def limit_request_size(request: Request, call_next):
    """Reject requests with body larger than limit."""
    content_length = request.headers.get("content-length")
    if content_length and int(content_length) > MAX_REQUEST_BODY:
        return JSONResponse(
            status_code=413,
            content={"detail": "Request body too large"},
        )
    return await call_next(request)


# ─── Include routers ──────────────────────────────────────────────────

app.include_router(router, prefix="/api")
app.include_router(websocket_router, prefix="/ws")


# ─── Health check (real) ──────────────────────────────────────────────

@app.get("/health")
async def health_check():
    """Health check with dependency verification."""
    health = {
        "status": "healthy",
        "version": "0.1.0",
        "uptime_seconds": round(time.time() - _metrics["uptime_start"], 1),
        "checks": {},
    }

    # Check Redis
    if _USE_REDIS:
        try:
            from api.routes import _redis
            _redis.ping()
            health["checks"]["redis"] = "ok"
        except Exception as e:
            health["checks"]["redis"] = f"error: {e}"
            health["status"] = "degraded"

    # Check disk space
    try:
        stat = os.statvfs(str(UPLOAD_DIR))
        free_gb = (stat.f_bavail * stat.f_frsize) / (1024 ** 3)
        health["checks"]["disk_free_gb"] = round(free_gb, 2)
        if free_gb < 1.0:
            health["status"] = "degraded"
            health["checks"]["disk"] = "low"
        else:
            health["checks"]["disk"] = "ok"
    except Exception:
        health["checks"]["disk"] = "unknown"

    # Running pipelines
    from api.routes import _running_pipelines, MAX_CONCURRENT_PIPELINES
    health["checks"]["running_pipelines"] = _running_pipelines
    health["checks"]["max_concurrent"] = MAX_CONCURRENT_PIPELINES

    return health


# ─── Metrics endpoint ─────────────────────────────────────────────────

@app.get("/metrics")
async def metrics():
    """Basic metrics for monitoring."""
    from api.routes import _running_pipelines

    return {
        "requests_total": _metrics["requests_total"],
        "requests_failed": _metrics["requests_failed"],
        "running_pipelines": _running_pipelines,
        "uptime_seconds": round(time.time() - _metrics["uptime_start"], 1),
    }


# ─── Main entry point ────────────────────────────────────────────────

def main():
    """Run the server (development mode)."""
    import uvicorn

    host = os.getenv("SPARC_HOST", "0.0.0.0")
    port = int(os.getenv("SPARC_PORT", "8000"))
    reload = os.getenv("SPARC_ENV", "production") == "development"

    logger.info("Starting SPARC API on %s:%d (reload=%s)", host, port, reload)

    uvicorn.run(
        "main:app",
        host=host,
        port=port,
        reload=reload,
    )


if __name__ == "__main__":
    main()
