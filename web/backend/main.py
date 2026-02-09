"""
FastAPI backend for sparc web interface.
"""

import os
from pathlib import Path
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from api.routes import router
from api.websocket import websocket_router


# Configuration
UPLOAD_DIR = Path(os.getenv("SCTOOLS_UPLOAD_DIR", "/tmp/sparc/uploads"))
OUTPUT_DIR = Path(os.getenv("SCTOOLS_OUTPUT_DIR", "/tmp/sparc/outputs"))


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize application state."""
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    yield


app = FastAPI(
    title="sparc API",
    description="Single-cell RNA-seq analysis API",
    version="0.1.0",
    lifespan=lifespan,
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(router, prefix="/api")
app.include_router(websocket_router, prefix="/ws")


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "version": "0.1.0"}


def main():
    """Run the server."""
    import uvicorn
    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
    )


if __name__ == "__main__":
    main()
