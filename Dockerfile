# ── Stage 1: Build Rust binaries ──────────────────────────────────────
FROM rust:1.77-bookworm AS rust-builder

RUN apt-get update && apt-get install -y \
    libhts-dev libdeflate-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    libssl-dev zlib1g-dev cmake pkg-config && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY crates/ crates/

# Build sparc CLI
RUN cargo build --release -p sparc-cli && \
    cp target/release/sparc /usr/local/bin/sparc

# ── Stage 2: Build Python wheel with maturin ─────────────────────────
FROM rust:1.77-bookworm AS py-builder

RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv python3-dev \
    libhts-dev libdeflate-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    libssl-dev zlib1g-dev cmake pkg-config && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --break-system-packages maturin

WORKDIR /build
COPY . .

RUN maturin build --release -o /wheels

# ── Stage 3: Production image ────────────────────────────────────────
FROM python:3.12-slim-bookworm AS production

RUN apt-get update && apt-get install -y --no-install-recommends \
    libhts3 libdeflate0 libbz2-1.0 liblzma5 libcurl4 zlib1g && \
    rm -rf /var/lib/apt/lists/*

# Copy sparc CLI binary
COPY --from=rust-builder /usr/local/bin/sparc /usr/local/bin/sparc

# Install Python wheel
COPY --from=py-builder /wheels/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl && rm -f /tmp/*.whl

# Install optional web dependencies
RUN pip install --no-cache-dir fastapi uvicorn celery redis python-multipart websockets

# Create directories
RUN mkdir -p /data/uploads /data/outputs /data/whitelists

ENV SCTOOLS_UPLOAD_DIR=/data/uploads
ENV SCTOOLS_OUTPUT_DIR=/data/outputs

WORKDIR /app
COPY web/ web/
COPY python/ python/

EXPOSE 8000

# Default: run the web API
CMD ["uvicorn", "web.backend.main:app", "--host", "0.0.0.0", "--port", "8000"]
