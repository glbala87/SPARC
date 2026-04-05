"""
SPARC Web UI launcher.

Provides the `sparc-web` console script entry point
that starts the FastAPI server with uvicorn.
"""

import os
import sys


def main():
    """Launch the SPARC web API server."""
    try:
        import uvicorn
    except ImportError:
        print(
            "Error: uvicorn is not installed. "
            "Install web dependencies with: pip install sparc-sc[web]",
            file=sys.stderr,
        )
        sys.exit(1)

    host = os.getenv("SPARC_HOST", "127.0.0.1")
    port = int(os.getenv("SPARC_PORT", "8000"))
    reload = os.getenv("SPARC_RELOAD", "false").lower() == "true"
    workers = int(os.getenv("SPARC_WORKERS", "1"))

    print(f"Starting SPARC Web UI at http://{host}:{port}")
    print("Press Ctrl+C to stop\n")

    uvicorn.run(
        "web.backend.main:app",
        host=host,
        port=port,
        reload=reload,
        workers=workers,
    )


if __name__ == "__main__":
    main()
