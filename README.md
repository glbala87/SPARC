# SPARC - Single-cell Pipeline Accelerated in Rust Core

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org/)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)

A high-performance single-cell RNA-seq analysis toolkit featuring:
- **Rust core** for blazing-fast FASTQ/BAM processing, barcode detection, and UMI deduplication
- **Python bindings** via PyO3 for seamless integration with Scanpy and AnnData
- **CLI** with 8 commands for end-to-end batch processing
- **Web API** (FastAPI) with real-time WebSocket progress, Redis job store, and Celery workers
- **Truthset validation** framework with synthetic data generation and accuracy metrics
- **Docker** support with multi-stage builds

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [CLI Reference](#cli-reference)
- [Python API](#python-api)
- [Truthset Validation](#truthset-validation)
- [Web Interface](#web-interface)
- [Docker Deployment](#docker-deployment)
- [Configuration](#configuration)
- [Production Deployment](#production-deployment)
- [Architecture](#architecture)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Prerequisites

### For Python Package Only
- Python 3.9+
- pip

### For Building from Source
- **Rust** 1.70+ ([Install Rust](https://rustup.rs/))
- **Python** 3.9+ with pip
- **C compiler** (gcc, clang, or MSVC)
- **CMake** 3.14+
- **zlib**, **bzip2**, **liblzma** development headers

### System-Specific Requirements

**macOS:**
```bash
xcode-select --install
brew install cmake zlib bzip2 htslib
```

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake zlib1g-dev libbz2-dev \
    liblzma-dev libcurl4-openssl-dev libssl-dev libhts-dev libdeflate-dev
```

**CentOS/RHEL:**
```bash
sudo yum groupinstall "Development Tools"
sudo yum install cmake zlib-devel bzip2-devel libcurl-devel
```

---

## Installation

### Quick Install (Python)

```bash
pip install sparc-sc
```

### From Source (Full)

```bash
git clone https://github.com/sparc-sc/sparc.git
cd sparc

# Install Rust (if needed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# Build Rust CLI
cargo build --release

# Build and install Python bindings
pip install maturin
maturin develop --release

# Verify
./target/release/sparc --version
python -c "import sparc; print(sparc.__version__); print('Rust:', sparc.check_rust_bindings())"
```

### Conda

```bash
conda env create -f environment.yml
conda activate sparc
maturin develop --release
```

### Docker

```bash
docker build -t sparc:latest .
docker run -v /path/to/data:/data sparc:latest sparc --help
```

---

## Quick Start

```bash
# Extract barcodes and UMIs
sparc extract -1 R1.fastq.gz -2 R2.fastq.gz -w whitelist.txt -o extracted/ --protocol 10x-3prime-v3

# Generate count matrix from BAM
sparc count -i aligned.bam -o counts/ --min-mapq 30

# Run QC
sparc qc -i counts/ -o qc_report.json --sample "MySample"

# Run downstream analysis (normalize, PCA, clustering)
sparc analyze -i counts/ -o analysis/ --n-pcs 50 --resolution 1.0

# Validate pipeline accuracy with synthetic data
sparc validate --n-cells 500 --n-genes 200 -o validation/
```

---

## CLI Reference

### Global Options

```
sparc [OPTIONS] <COMMAND>

Options:
  -v, --verbose     Enable verbose output
  -j, --threads     Number of threads (0 = auto-detect)
  -h, --help        Print help
  -V, --version     Print version
```

### Commands

| Command | Description |
|---------|-------------|
| `extract` | Extract barcodes and UMIs from paired-end FASTQ files |
| `count` | Generate gene-by-cell count matrix from aligned BAM |
| `qc` | Generate quality control metrics and report |
| `pipeline` | Run the full analysis pipeline (extract + align + count + QC) |
| `analyze` | Run downstream analysis (normalize, PCA, KNN, clustering) |
| `batch` | Process multiple samples from a manifest file |
| `distributed` | Distributed processing (shard/worker/merge) |
| `validate` | Run truthset validation against synthetic ground-truth data |

### `sparc extract`

```bash
sparc extract -1 <R1> -2 <R2> -w <WHITELIST> -o <OUTPUT> [OPTIONS]

Options:
  -p, --protocol <PROTOCOL>    Protocol [default: 10x-3prime-v3]
      --max-mismatch <N>       Max Hamming distance for correction [default: 1]
      --min-barcode-qual <N>   Min barcode quality score [default: 10]
```

Supported protocols: `10x-3prime-v3`, `10x-3prime-v2`, `10x-5prime-v2`, `drop-seq`, `indrop`, `sci-rna-seq`, `smart-seq2`

### `sparc count`

```bash
sparc count -i <BAM> -o <OUTPUT> [OPTIONS]

Options:
      --min-mapq <N>    Minimum mapping quality [default: 30]
      --format <FMT>    Output format: mtx, h5ad [default: mtx]
```

### `sparc qc`

```bash
sparc qc -i <INPUT_DIR> -o <OUTPUT_FILE> [OPTIONS]

Options:
  -s, --sample <NAME>   Sample name [default: sample]
      --min-genes <N>   Min genes per cell [default: 200]
      --max-genes <N>   Max genes per cell [default: 10000]
      --max-mito <F>    Max mitochondrial % [default: 20.0]
```

### `sparc analyze`

Run downstream analysis on a count matrix directory.

```bash
sparc analyze -i <INPUT_DIR> -o <OUTPUT_DIR> [OPTIONS]

Options:
      --n-pcs <N>          Number of PCA components [default: 50]
      --n-neighbors <N>    KNN neighbors [default: 15]
      --resolution <F>     Clustering resolution [default: 1.0]
      --target-sum <F>     Normalization target [default: 10000]

Output:
  clusters.tsv            Barcode-to-cluster assignments
  pca.tsv                 PCA coordinates per cell
  hvg.tsv                 Highly variable gene list
  analysis_summary.json   Summary statistics
```

### `sparc pipeline`

Run the complete pipeline (extract + align + count + QC).

```bash
sparc pipeline -1 <R1> -2 <R2> -r <REFERENCE> -w <WHITELIST> -o <OUTPUT> [OPTIONS]

Options:
  -p, --protocol <PROTOCOL>  Protocol [default: 10x-3prime-v3]
      --aligner <ALIGNER>    star or minimap2 [default: star]
      --skip-align           Skip alignment (use --bam for pre-aligned)
      --bam <FILE>           Pre-aligned BAM file
```

### `sparc validate`

Run truthset validation with synthetic ground-truth data.

```bash
sparc validate [OPTIONS]

Options:
  -o, --output <DIR>           Output directory [default: validation_output]
      --n-cells <N>            Synthetic cells [default: 500]
      --n-genes <N>            Synthetic genes [default: 200]
      --n-cell-types <N>       Cell types [default: 5]
      --stages <STAGES>        extract,count,analysis or all [default: all]
      --seed <N>               Random seed [default: 42]
      --min-barcode-f1 <F>     Pass threshold for barcode F1 [default: 0.95]
      --min-expression-pearson <F>  Pass threshold for Pearson r [default: 0.90]
      --min-clustering-ari <F>      Pass threshold for ARI [default: 0.70]

Output:
  validation_report.json   Full metrics report with PASS/FAIL verdicts
```

### `sparc distributed`

Distributed processing across multiple machines.

```bash
# Coordinator: prints shard commands
sparc distributed -1 R1.fq.gz -2 R2.fq.gz -r ref/ -o out/ -w wl.txt --shards 4

# Worker: process one shard
sparc distributed -1 R1.fq.gz -2 R2.fq.gz -r ref/ -o out/ -w wl.txt --shards 4 --shard-index 0 --worker

# Merge: combine results
sparc distributed -1 R1.fq.gz -2 R2.fq.gz -r ref/ -o out/ -w wl.txt --shards 4 --merge
```

---

## Python API

### Installation

```bash
pip install sparc-sc              # Core
pip install sparc-sc[web]         # + Web UI dependencies
pip install sparc-sc[dev]         # + Development tools
pip install sparc-sc[all]         # Everything
```

### Basic Usage

```python
import sparc

print(f"SPARC version: {sparc.__version__}")
print(f"Rust bindings: {sparc.check_rust_bindings()}")
```

### Reading FASTQ Files

```python
for record in sparc.read_fastq("sample_R1.fastq.gz"):
    barcode = record.subsequence(0, 16)
    umi = record.subsequence(16, 12)
```

### Barcode Correction

```python
whitelist = sparc.Whitelist("whitelist.txt")
corrector = sparc.BarcodeCorrector(whitelist, max_distance=1)
status, corrected, distance = corrector.match_barcode("AAACCCAAGAAACACT")
```

### Count Matrix + Scanpy

```python
import sparc
import scanpy as sc

matrix, barcodes, genes = sparc.read_matrix("counts/")
adata = sparc.to_anndata(matrix, barcodes, genes)
adata = sparc.normalize_and_analyze(adata)
sc.pl.umap(adata, color='leiden')
```

### Truthset Validation (Python)

```python
from sparc.validation import run_validation, SyntheticConfig

config = SyntheticConfig(n_cells=500, n_genes=200, n_cell_types=5)
report = run_validation(config=config, stages="all")

print(report.summary())
print(f"Overall: {'PASS' if report.overall_pass else 'FAIL'}")
```

---

## Truthset Validation

SPARC includes a built-in truthset validation framework that generates synthetic scRNA-seq data with known ground truth and validates each pipeline stage.

### What It Tests

| Stage | Metrics | Threshold |
|-------|---------|-----------|
| **Extract** (barcode detection) | Sensitivity, Specificity, Precision, Recall, F1, Correction accuracy | F1 >= 0.95 |
| **Count** (expression quantification) | Pearson r, Spearman rho, MAE, RMSE, cell/gene concordance | Pearson >= 0.90 |
| **Analysis** (clustering) | Adjusted Rand Index (ARI), Normalized Mutual Information (NMI) | ARI >= 0.70 |

### How It Works

1. Generates synthetic cells with known cell types and marker genes (Poisson-sampled expression)
2. Creates synthetic FASTQ reads with valid, mutated, and invalid barcodes
3. Runs each pipeline stage against the synthetic data
4. Compares output to ground truth and computes accuracy metrics
5. Produces a JSON report with PASS/FAIL verdicts

```bash
sparc validate --n-cells 1000 --n-genes 500 --n-cell-types 8 -o validation/
```

---

## Web Interface

### Starting the Web API

```bash
cd web/backend
pip install fastapi uvicorn python-multipart websockets redis
uvicorn main:app --host 0.0.0.0 --port 8000
```

### API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/upload` | POST | Upload FASTQ/whitelist files |
| `/api/pipeline/{job_id}` | POST | Start pipeline |
| `/api/pipeline/{job_id}/status` | GET | Get pipeline status |
| `/api/pipeline/{job_id}/results` | GET | Get results |
| `/api/pipeline/{job_id}` | DELETE | Delete job and data |
| `/api/pipeline/{job_id}/notebook` | GET | Export as Jupyter notebook |
| `/api/protocols` | GET | List available protocols |
| `/api/whitelists` | GET | List available whitelists |
| `/api/whitelists/upload` | POST | Upload custom whitelist |
| `/api/jobs` | GET | List all jobs |
| `/api/compare` | POST | Compare multiple samples |
| `/ws/pipeline/{job_id}` | WebSocket | Real-time progress updates |
| `/health` | GET | Health check (Redis, disk, pipelines) |
| `/metrics` | GET | Request/pipeline metrics |

---

## Docker Deployment

### Quick Start

```bash
# Build and run with Docker Compose (API + Redis + worker)
docker compose up -d

# Check health
curl http://localhost:8000/health

# Run CLI in Docker
docker run -v $(pwd)/data:/data sparc:latest sparc validate -o /data/validation
```

### docker-compose.yml

Includes three services:
- **redis** — Job store and Celery broker
- **sparc-api** — FastAPI server on port 8000
- **sparc-worker** — Celery background worker

---

## Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| **Security** | | |
| `SPARC_API_TOKENS` | Comma-separated API tokens (empty = auth disabled) | **(required)** |
| `SPARC_CORS_ORIGINS` | Allowed CORS origins | `http://localhost:3000,http://localhost:5173` |
| **Storage** | | |
| `SPARC_UPLOAD_DIR` | Upload directory | `/tmp/sparc/uploads` |
| `SPARC_OUTPUT_DIR` | Output directory | `/tmp/sparc/outputs` |
| `REDIS_URL` | Redis connection URL | *(empty = in-memory)* |
| **Limits** | | |
| `SPARC_MAX_UPLOAD_MB` | Max FASTQ upload size (MB) | `5000` |
| `SPARC_MAX_WHITELIST_MB` | Max whitelist upload size (MB) | `500` |
| `SPARC_MAX_REQUEST_MB` | Max total request body (MB) | `6000` |
| `SPARC_MAX_JOBS` | Max concurrent job records | `100` |
| `SPARC_MAX_CONCURRENT` | Max concurrent pipelines | `4` |
| `SPARC_WS_MAX_PER_JOB` | Max WebSocket connections per job | `10` |
| `SPARC_WS_MAX_TOTAL` | Max total WebSocket connections | `100` |
| **Server** | | |
| `SPARC_HOST` | API bind host | `0.0.0.0` |
| `SPARC_PORT` | API bind port | `8000` |
| `SPARC_ENV` | `development` enables auto-reload | `production` |
| `SPARC_LOG_LEVEL` | Log level (DEBUG, INFO, WARNING, ERROR) | `INFO` |

---

## Production Deployment

### Security Checklist

- [ ] Set `SPARC_API_TOKENS` to strong, unique tokens
- [ ] Set `SPARC_CORS_ORIGINS` to your frontend domain only
- [ ] Use non-`/tmp` directories for `SPARC_UPLOAD_DIR` / `SPARC_OUTPUT_DIR`
- [ ] Deploy Redis with authentication (`REDIS_URL=redis://:password@host:6379/0`)
- [ ] Run behind a reverse proxy (nginx/Caddy) with TLS
- [ ] Set `SPARC_ENV=production` (disables auto-reload)

### Monitoring

- **Health**: `GET /health` — checks Redis connectivity, disk space, running pipelines
- **Metrics**: `GET /metrics` — request counts, failure counts, uptime, running pipelines
- **Logs**: Structured logging with timestamps, log levels, and module names. Set `SPARC_LOG_LEVEL=DEBUG` for verbose output.

### Graceful Shutdown

The API handles `SIGTERM`/`SIGINT` signals:
1. Stops accepting new requests
2. Waits up to 30 seconds for running pipelines to complete
3. Closes all WebSocket connections
4. Exits cleanly

---

## Architecture

```
+------------------------------------------------------------------+
|                         Web UI (React)                            |
|              Upload . Configure . Visualize . Export              |
+------------------------------------------------------------------+
|                    FastAPI Backend (Python)                       |
|  REST API . WebSocket . Celery . Redis . Auth . Rate Limiting    |
+------------------------------------------------------------------+
|                 Python Package (sparc-py)                         |
|         +------------------------------------------+              |
|         |         PyO3 Bindings Layer              |              |
|         |  FASTQ . BAM . Barcode . Matrix . QC     |              |
|         |  Analysis . Validation . Metrics          |              |
+------------------------------------------------------------------+
|                     Rust Core (sparc-core)                        |
|  +----------+ +----------+ +----------+ +----------+             |
|  |  FASTQ   | |   BAM    | | Barcode  | |   UMI    |             |
|  |  Parser  | |  Parser  | | Matcher  | |  Dedup   |             |
|  |needletail| |rust-htslib| |  ahash   | | graph    |             |
|  +----------+ +----------+ +----------+ +----------+             |
|  +----------+ +----------+ +----------+ +----------+             |
|  | Analysis | |Validation| |    QC    | |  Count   |             |
|  |PCA/KNN/  | | Truthset | | Metrics  | |  Matrix  |             |
|  |Clustering| | Synthetic| |  Report  | |  COO/CSR |             |
|  +----------+ +----------+ +----------+ +----------+             |
+------------------------------------------------------------------+
```

### Project Structure

```
sparc/
├── Cargo.toml                 # Rust workspace
├── pyproject.toml             # Python package config
├── Dockerfile                 # Multi-stage Docker build
├── docker-compose.yml         # API + Redis + Worker
├── environment.yml            # Conda environment
│
├── crates/
│   ├── sparc-core/            # Core Rust library
│   │   └── src/
│   │       ├── fastq/         # FASTQ parsing (needletail)
│   │       ├── bam/           # BAM parsing (rust-htslib)
│   │       ├── barcode/       # Barcode matching + correction
│   │       ├── umi/           # UMI deduplication
│   │       ├── protocols/     # 10x/Drop-seq/inDrop/sci-RNA/Smart-seq2
│   │       ├── qc/            # Quality control metrics
│   │       ├── count/         # Count matrix (COO/CSR)
│   │       ├── analysis/      # Normalize, PCA, KNN, clustering
│   │       ├── validation/    # Truthset validation framework
│   │       ├── aligner.rs     # STAR/minimap2 integration
│   │       └── streaming.rs   # Streaming processor
│   │
│   ├── sparc-cli/             # CLI application (8 commands)
│   │   └── src/commands/      # extract, count, qc, pipeline,
│   │                          # analyze, batch, distributed, validate
│   │
│   └── sparc-py/              # PyO3 Python bindings
│       └── src/               # fastq, bam, barcode, matrix,
│                              # analysis, qc, validation
│
├── python/sparc/              # Python package
│   ├── io.py                  # I/O (FASTQ, BAM, MTX, H5AD)
│   ├── preprocessing.py       # Barcode extraction, UMI dedup
│   ├── analysis.py            # Scanpy integration
│   ├── validation/            # Truthset validation (Python)
│   │   ├── synthetic.py       # Synthetic dataset generator
│   │   ├── metrics.py         # Accuracy metrics
│   │   └── runner.py          # Validation orchestrator
│   ├── streaming.py           # Streaming processor
│   ├── batch.py               # Batch processing
│   └── web.py                 # Web UI launcher (sparc-web)
│
├── web/
│   ├── backend/               # FastAPI server
│   │   ├── main.py            # App, middleware, health, metrics
│   │   ├── api/routes.py      # REST endpoints
│   │   ├── api/websocket.py   # WebSocket progress
│   │   └── workers/pipeline.py # Celery worker
│   └── frontend/              # React UI
│
└── tests/
    ├── python/                # Python test suite
    └── sparc-core/tests/      # Rust integration tests
```

### Performance

| Operation | SPARC | Cell Ranger | Speedup |
|-----------|-------|-------------|---------|
| FASTQ parsing | 2.1M reads/s | 0.8M reads/s | 2.6x |
| Barcode matching | 1.5M/s | 0.4M/s | 3.8x |
| UMI dedup | 0.9M/s | 0.3M/s | 3.0x |

*Benchmarked on 8-core Intel i7 with 32GB RAM*

---

## Troubleshooting

### Common Issues

#### "Rust bindings not available"
```bash
pip uninstall sparc-sc
maturin develop --release
```

#### "Failed to open BAM file"
- Ensure the BAM file is indexed: `samtools index file.bam`
- Check file permissions

#### "No valid barcodes found"
- Verify whitelist matches your chemistry version
- Check that R1 contains barcodes (not R2)
- Try increasing `--max-mismatch`

#### Build errors on Linux
```bash
sudo apt-get install libssl-dev pkg-config libhts-dev libdeflate-dev
```

#### Memory issues with large files
```bash
sparc extract -j 4 ...  # Limit threads
```

### Debug Mode

```bash
sparc -v extract ...
RUST_LOG=debug sparc extract ...
SPARC_LOG_LEVEL=DEBUG sparc-web  # Web API debug logging
```

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Citation

```bibtex
@software{sparc2024,
  title = {SPARC: Single-cell Pipeline Accelerated in Rust Core},
  author = {SPARC Contributors},
  year = {2024},
  url = {https://github.com/sparc-sc/sparc},
  version = {0.1.0}
}
```

---

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Support

- **Issues**: [GitHub Issues](https://github.com/sparc-sc/sparc/issues)
- **Discussions**: [GitHub Discussions](https://github.com/sparc-sc/sparc/discussions)
