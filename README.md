# SPARC - Single-cell Pipeline Accelerated in Rust Core

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org/)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)

A high-performance single-cell RNA-seq analysis toolkit featuring:
- **Rust core** for blazing-fast FASTQ/BAM processing, barcode detection, and UMI deduplication
- **Python bindings** via PyO3 for seamless integration with Scanpy and AnnData
- **CLI** for command-line batch processing
- **Web UI** for interactive analysis and visualization

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
  - [Quick Install (Python)](#quick-install-python)
  - [From Source (Full)](#from-source-full)
  - [Docker](#docker)
- [Quick Start](#quick-start)
- [CLI Reference](#cli-reference)
- [Python API](#python-api)
- [Web Interface](#web-interface)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [Architecture](#architecture)
- [License](#license)

---

## Prerequisites

### For Python Package Only
- Python 3.9 or higher
- pip

### For Building from Source
- **Rust** 1.70+ ([Install Rust](https://rustup.rs/))
- **Python** 3.9+ with pip
- **C compiler** (gcc, clang, or MSVC)
- **CMake** 3.14+ (for rust-htslib)
- **zlib** development headers
- **bzip2** development headers (optional, for BAM support)

### System-Specific Requirements

**macOS:**
```bash
# Install Xcode command line tools
xcode-select --install

# Install dependencies via Homebrew
brew install cmake zlib bzip2
```

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake zlib1g-dev libbz2-dev libcurl4-openssl-dev
```

**CentOS/RHEL:**
```bash
sudo yum groupinstall "Development Tools"
sudo yum install cmake zlib-devel bzip2-devel libcurl-devel
```

**Windows:**
- Install [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
- Install [CMake](https://cmake.org/download/)

---

## Installation

### Quick Install (Python)

```bash
pip install sparc-sc
```

This installs pre-built wheels with all Rust components included.

### From Source (Full)

#### Step 1: Clone the Repository
```bash
git clone https://github.com/sparc-sc/sparc.git
cd sparc
```

#### Step 2: Install Rust (if not already installed)
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

#### Step 3: Build Rust Components
```bash
# Build in release mode for best performance
cargo build --release

# The CLI binary will be at: target/release/sparc
```

#### Step 4: Install Python Package
```bash
# Install maturin (Python-Rust build tool)
pip install maturin

# Build and install Python bindings
maturin develop --release

# Or install with optional dependencies
pip install -e ".[all]"
```

#### Step 5: Verify Installation
```bash
# Check CLI
./target/release/sparc --version

# Check Python bindings
python -c "import sparc; print(sparc.__version__)"
```

### Docker

```bash
# Build the Docker image
docker build -t sparc:latest .

# Run SPARC in Docker
docker run -v /path/to/data:/data sparc:latest sparc extract -1 /data/R1.fq.gz -2 /data/R2.fq.gz -o /data/output
```

---

## Quick Start

### Example Workflow

```bash
# 1. Extract barcodes and UMIs from FASTQ files
sparc extract \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -w data/whitelists/3M-february-2018.txt \
    -o extracted/ \
    --protocol 10x-3prime-v3

# 2. Align reads with STAR (external tool)
STAR --genomeDir /path/to/star_index \
     --readFilesIn extracted/R2_filtered.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix aligned/

# 3. Generate count matrix from aligned BAM
sparc count \
    -i aligned/Aligned.sortedByCoord.out.bam \
    -o counts/ \
    --min-mapq 30

# 4. Generate QC report
sparc qc \
    -i counts/ \
    -o qc_report.json \
    --sample "MySample" \
    --min-genes 200 \
    --max-genes 10000 \
    --max-mito 20
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

#### `sparc extract` - Extract Barcodes and UMIs

Extract cell barcodes and UMIs from paired-end FASTQ files.

```bash
sparc extract [OPTIONS] -1 <R1> -2 <R2> -w <WHITELIST> -o <OUTPUT>

Required Arguments:
  -1, --r1 <R1>              R1 FASTQ file (barcode/UMI read)
  -2, --r2 <R2>              R2 FASTQ file (cDNA read)
  -w, --whitelist <FILE>     Barcode whitelist file
  -o, --output <DIR>         Output directory

Options:
  -p, --protocol <PROTOCOL>  Protocol [default: 10x-3prime-v3]
                             Values: 10x-3prime-v3, 10x-3prime-v2, 10x-5prime-v2
      --max-mismatch <N>     Max Hamming distance for correction [default: 1]
      --min-barcode-qual <N> Min barcode quality score [default: 10]

Examples:
  sparc extract -1 R1.fq.gz -2 R2.fq.gz -w whitelist.txt -o out/
  sparc extract -1 R1.fq.gz -2 R2.fq.gz -w whitelist.txt -o out/ --protocol 10x-3prime-v2
```

#### `sparc count` - Generate Count Matrix

Generate a gene-by-cell count matrix from aligned BAM file.

```bash
sparc count [OPTIONS] -i <INPUT> -o <OUTPUT>

Required Arguments:
  -i, --input <BAM>     Input BAM file (with CB, UB, GN/GX tags)
  -o, --output <DIR>    Output directory

Options:
      --min-mapq <N>    Minimum mapping quality [default: 30]
      --format <FMT>    Output format: mtx, h5ad [default: mtx]

Examples:
  sparc count -i aligned.bam -o counts/
  sparc count -i aligned.bam -o counts/ --min-mapq 20 --format mtx
```

#### `sparc qc` - Generate QC Report

Generate quality control metrics and report.

```bash
sparc qc [OPTIONS] -i <INPUT> -o <OUTPUT>

Required Arguments:
  -i, --input <DIR>     Input matrix directory
  -o, --output <FILE>   Output JSON report file

Options:
  -s, --sample <NAME>   Sample name [default: sample]
      --min-genes <N>   Min genes per cell [default: 200]
      --max-genes <N>   Max genes per cell [default: 10000]
      --max-mito <F>    Max mitochondrial % [default: 20.0]

Examples:
  sparc qc -i counts/ -o qc_report.json
  sparc qc -i counts/ -o qc.json --sample "Patient1" --min-genes 500
```

#### `sparc pipeline` - Run Full Pipeline

Run the complete analysis pipeline (experimental).

```bash
sparc pipeline [OPTIONS] -1 <R1> -2 <R2> -r <REFERENCE> -w <WHITELIST> -o <OUTPUT>

Required Arguments:
  -1, --r1 <R1>              R1 FASTQ file
  -2, --r2 <R2>              R2 FASTQ file
  -r, --reference <DIR>      Reference genome directory
  -w, --whitelist <FILE>     Barcode whitelist
  -o, --output <DIR>         Output directory

Options:
  -p, --protocol <PROTOCOL>  Protocol [default: 10x-3prime-v3]
  -s, --sample <NAME>        Sample name [default: sample]
      --expect-cells <N>     Expected number of cells
      --force-cells <N>      Force number of cells
```

---

## Python API

### Installation

```bash
pip install sparc-sc

# With optional dependencies
pip install sparc-sc[scanpy]  # Scanpy integration
pip install sparc-sc[web]     # Web UI dependencies
pip install sparc-sc[all]     # Everything
```

### Basic Usage

```python
import sparc

# Check if Rust bindings are available
print(f"SPARC version: {sparc.__version__}")
print(f"Rust bindings: {sparc.check_rust_bindings()}")
```

### Reading FASTQ Files

```python
import sparc

# Iterate over FASTQ records (supports .gz compression)
for record in sparc.read_fastq("sample_R1.fastq.gz"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq_str()}")
    print(f"Quality: {record.mean_quality():.2f}")

    # Extract subsequences
    barcode = record.subsequence(0, 16)  # First 16bp
    umi = record.subsequence(16, 12)     # Next 12bp
```

### Barcode Correction

```python
import sparc

# Load barcode whitelist
whitelist = sparc.Whitelist("data/whitelists/3M-february-2018.txt")
print(f"Loaded {len(whitelist)} barcodes")

# Create corrector with max 1 mismatch
corrector = sparc.BarcodeCorrector(whitelist, max_distance=1)

# Correct a barcode
status, corrected, distance = corrector.match_barcode("AAACCCAAGAAACACT")
print(f"Status: {status}, Corrected: {corrected}, Distance: {distance}")

# Batch correction
barcodes = ["AAACCCAAGAAACACT", "AAACCCAAGAAACCAT", "INVALIDBARCODE"]
corrected = corrector.correct_batch(barcodes)
```

### Reading BAM Files

```python
import sparc

# Iterate over BAM records
for record in sparc.read_bam("aligned.bam", min_mapq=30):
    print(f"Name: {record.name}")
    print(f"Cell Barcode: {record.cell_barcode}")
    print(f"UMI: {record.umi}")
    print(f"Gene: {record.gene_name}")
    print(f"Mapped: {record.is_mapped}, MAPQ: {record.mapq}")
```

### Building Count Matrix

```python
import sparc

# Create gene counter
counter = sparc.GeneCounter()

# Add counts
counter.increment("CELL1-1", "TP53")
counter.increment("CELL1-1", "TP53")  # Duplicate UMI
counter.increment("CELL1-1", "BRCA1")
counter.increment("CELL2-1", "TP53")

# Build matrix
matrix = counter.build()
print(f"Matrix: {matrix.n_rows} genes x {matrix.n_cols} cells")
print(f"Non-zero entries: {matrix.nnz}")

# Get statistics
print(f"Counts per cell: {matrix.counts_per_cell()}")
print(f"Genes per cell: {matrix.genes_per_cell()}")

# Export to Matrix Market format
matrix.write_mtx("output/matrix.mtx")
matrix.write_barcodes("output/barcodes.tsv")
matrix.write_genes("output/genes.tsv")
```

### Scanpy Integration

```python
import sparc
import scanpy as sc

# Read matrix and convert to AnnData
matrix, barcodes, genes = sparc.read_matrix("counts/")
adata = sparc.to_anndata(matrix, barcodes, genes)

# Or run the full pipeline
result = sparc.run_pipeline(
    "counts/",
    min_genes=200,
    max_genes=10000,
    max_mito=20.0,
    n_pcs=50,
    resolution=1.0,
)

print(f"Cells passing QC: {result.qc_passed}")
print(f"Final cells: {result.total_cells}")
print(f"Genes: {result.total_genes}")

# Access the processed AnnData
adata = result.adata
sc.pl.umap(adata, color='leiden')
```

### Complete Example

```python
import sparc

# Step 1: Extract barcodes
result = sparc.extract_barcodes(
    r1_path="sample_R1.fastq.gz",
    whitelist_path="whitelist.txt",
    barcode_start=0,
    barcode_len=16,
    umi_start=16,
    umi_len=12,
    max_mismatch=1,
)
print(f"Total reads: {result.total_reads}")
print(f"Valid barcodes: {result.valid_barcodes} ({result.valid_barcodes/result.total_reads*100:.1f}%)")
print(f"Corrected: {result.corrected_barcodes}")

# Step 2: Build count matrix from extracted data
counter = sparc.GeneCounter()
for barcode, gene in zip(result.barcodes, genes_from_alignment):
    counter.increment(barcode, gene)

matrix = counter.build()

# Step 3: Convert to AnnData and analyze
adata = sparc.to_anndata(matrix)
```

---

## Web Interface

### Starting the Web UI

#### Backend (FastAPI)

```bash
cd web/backend

# Install dependencies
pip install fastapi uvicorn python-multipart websockets

# Start server
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
```

#### Frontend (React)

```bash
cd web/frontend

# Install dependencies
npm install

# Start development server
npm run dev
```

Open http://localhost:3000 in your browser.

### Features

1. **Upload**: Drag-and-drop FASTQ files and whitelist
2. **Configure**: Select protocol and set parameters
3. **Pipeline**: Run analysis with real-time progress
4. **QC Report**: Interactive quality control metrics
5. **Visualization**: UMAP plots with cluster exploration

### API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/upload` | POST | Upload FASTQ/whitelist files |
| `/api/pipeline/{job_id}` | POST | Start pipeline |
| `/api/pipeline/{job_id}/status` | GET | Get pipeline status |
| `/api/pipeline/{job_id}/results` | GET | Get results |
| `/api/protocols` | GET | List available protocols |
| `/api/whitelists` | GET | List available whitelists |
| `/ws/pipeline/{job_id}` | WebSocket | Real-time progress |

---

## Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `SPARC_UPLOAD_DIR` | Upload directory | `/tmp/sparc/uploads` |
| `SPARC_OUTPUT_DIR` | Output directory | `/tmp/sparc/outputs` |
| `SPARC_THREADS` | Number of threads | Auto-detect |
| `REDIS_URL` | Redis URL (for Celery) | `redis://localhost:6379/0` |

### Barcode Whitelists

Place whitelist files in `data/whitelists/`:
- `3M-february-2018.txt` - 10x Genomics v3 (6.7M barcodes)
- `737K-august-2016.txt` - 10x Genomics v2 (737K barcodes)

Download from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads):
```bash
# v3 whitelist
curl -O https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
gunzip 3M-february-2018.txt.gz

# v2 whitelist
curl -O https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt
```

---

## Troubleshooting

### Common Issues

#### "Rust bindings not available"
```bash
# Rebuild Python bindings
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
# Install missing dependencies
sudo apt-get install libssl-dev pkg-config
```

#### Memory issues with large files
```bash
# Limit threads to reduce memory
sparc extract -j 4 ...
```

### Debug Mode

```bash
# Enable verbose logging
sparc -v extract ...

# Set Rust log level
RUST_LOG=debug sparc extract ...
```

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         Web UI (React)                          │
│              Upload • Configure • Visualize • Export            │
├─────────────────────────────────────────────────────────────────┤
│                    FastAPI Backend (Python)                     │
│                REST API • WebSocket • Celery                    │
├─────────────────────────────────────────────────────────────────┤
│                 Python Package (sparc-py)                       │
│         ┌─────────────────────────────────────────┐             │
│         │         PyO3 Bindings Layer             │             │
│         │   NumPy arrays • Iterators • Types      │             │
├─────────┴─────────────────────────────────────────┴─────────────┤
│                     Rust Core (sparc-core)                      │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐           │
│  │  FASTQ   │ │   BAM    │ │ Barcode  │ │   UMI    │           │
│  │  Parser  │ │  Parser  │ │ Matcher  │ │  Dedup   │           │
│  │needletail│ │rust-htslib│ │  ahash   │ │ graph    │           │
│  └──────────┘ └──────────┘ └──────────┘ └──────────┘           │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐           │
│  │ 10x 3'   │ │  10x 5'  │ │    QC    │ │  Count   │           │
│  │ Protocol │ │ Protocol │ │ Metrics  │ │  Matrix  │           │
│  └──────────┘ └──────────┘ └──────────┘ └──────────┘           │
└─────────────────────────────────────────────────────────────────┘
```

### Performance

| Operation | SPARC | Cell Ranger | Speedup |
|-----------|-------|-------------|---------|
| FASTQ parsing | 2.1M reads/s | 0.8M reads/s | 2.6x |
| Barcode matching | 1.5M/s | 0.4M/s | 3.8x |
| UMI dedup | 0.9M/s | 0.3M/s | 3.0x |

*Benchmarked on 8-core Intel i7 with 32GB RAM*

---

## Project Structure

```
sparc/
├── Cargo.toml                 # Rust workspace
├── pyproject.toml             # Python package config
├── README.md                  # This file
│
├── crates/
│   ├── sparc-core/            # Core Rust library
│   │   └── src/
│   │       ├── fastq/         # FASTQ parsing
│   │       ├── bam/           # BAM parsing
│   │       ├── barcode/       # Barcode matching
│   │       ├── umi/           # UMI deduplication
│   │       ├── protocols/     # 10x kit implementations
│   │       ├── qc/            # Quality metrics
│   │       └── count/         # Count matrix
│   │
│   ├── sparc-cli/             # CLI application
│   │   └── src/commands/      # extract, count, qc, pipeline
│   │
│   └── sparc-py/              # PyO3 Python bindings
│
├── python/sparc/              # Python wrapper
│   ├── io.py                  # I/O functions
│   ├── preprocessing.py       # Preprocessing
│   ├── analysis.py            # Scanpy integration
│   └── plotting.py            # Visualization
│
├── web/
│   ├── backend/               # FastAPI server
│   └── frontend/              # React UI
│
├── data/whitelists/           # Barcode whitelists
└── tests/                     # Test suites
```

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Citation

If you use SPARC in your research, please cite:

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
