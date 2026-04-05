"""
Multi-sample batch processing.
"""

import csv
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union


@dataclass
class SampleConfig:
    """Configuration for a single sample."""
    name: str
    r1: Path
    r2: Path
    whitelist: Path
    protocol: str = "10x-3prime-v3"


@dataclass
class BatchResult:
    """Result from batch processing."""
    sample: str
    success: bool
    message: str = ""
    output_dir: Optional[Path] = None


@dataclass
class BatchSummary:
    """Summary of batch processing results."""
    total: int = 0
    succeeded: int = 0
    failed: int = 0
    results: list[BatchResult] = field(default_factory=list)


def parse_manifest(path: Union[str, Path]) -> list[SampleConfig]:
    """
    Parse a sample manifest CSV file.

    Expected columns: sample_name, r1, r2, whitelist[, protocol]

    Parameters
    ----------
    path : str or Path
        Path to manifest CSV file

    Returns
    -------
    list of SampleConfig
        Parsed sample configurations
    """
    samples = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("sample_name", row.get("name", ""))
            if not name:
                continue
            samples.append(SampleConfig(
                name=name.strip(),
                r1=Path(row["r1"].strip()),
                r2=Path(row["r2"].strip()),
                whitelist=Path(row["whitelist"].strip()),
                protocol=row.get("protocol", "10x-3prime-v3").strip(),
            ))
    return samples


def _run_sample(
    sample: SampleConfig,
    reference: Path,
    output_dir: Path,
    aligner: str = "star",
    max_mismatch: int = 1,
) -> BatchResult:
    """Run the SPARC pipeline on a single sample."""
    sample_output = output_dir / sample.name
    try:
        cmd = [
            "sparc", "pipeline",
            "-1", str(sample.r1),
            "-2", str(sample.r2),
            "-r", str(reference),
            "-o", str(sample_output),
            "-w", str(sample.whitelist),
            "-p", sample.protocol,
            "--aligner", aligner,
            "--max-mismatch", str(max_mismatch),
            "-s", sample.name,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return BatchResult(
            sample=sample.name,
            success=True,
            message="Pipeline completed successfully",
            output_dir=sample_output,
        )
    except subprocess.CalledProcessError as e:
        return BatchResult(
            sample=sample.name,
            success=False,
            message=f"Pipeline failed: {e.stderr}",
        )
    except FileNotFoundError:
        return BatchResult(
            sample=sample.name,
            success=False,
            message="sparc CLI not found in PATH",
        )


class BatchProcessor:
    """
    Process multiple samples from a manifest file.

    Parameters
    ----------
    manifest : str or Path
        Path to sample manifest CSV
    reference : str or Path
        Path to reference genome directory
    output : str or Path
        Output directory
    parallel : int
        Number of samples to process in parallel (default: 1)
    aligner : str
        Aligner to use (default: "star")
    max_mismatch : int
        Maximum barcode mismatch (default: 1)
    """

    def __init__(
        self,
        manifest: Union[str, Path],
        reference: Union[str, Path],
        output: Union[str, Path],
        parallel: int = 1,
        aligner: str = "star",
        max_mismatch: int = 1,
    ):
        self.manifest = Path(manifest)
        self.reference = Path(reference)
        self.output = Path(output)
        self.parallel = parallel
        self.aligner = aligner
        self.max_mismatch = max_mismatch

    def run(self) -> BatchSummary:
        """
        Run the batch processing pipeline.

        Returns
        -------
        BatchSummary
            Summary of all sample results
        """
        samples = parse_manifest(self.manifest)
        self.output.mkdir(parents=True, exist_ok=True)

        summary = BatchSummary(total=len(samples))

        if self.parallel > 1:
            with ProcessPoolExecutor(max_workers=self.parallel) as executor:
                futures = {
                    executor.submit(
                        _run_sample, sample, self.reference, self.output,
                        self.aligner, self.max_mismatch
                    ): sample
                    for sample in samples
                }
                for future in as_completed(futures):
                    result = future.result()
                    summary.results.append(result)
                    if result.success:
                        summary.succeeded += 1
                    else:
                        summary.failed += 1
        else:
            for sample in samples:
                result = _run_sample(
                    sample, self.reference, self.output,
                    self.aligner, self.max_mismatch
                )
                summary.results.append(result)
                if result.success:
                    summary.succeeded += 1
                else:
                    summary.failed += 1

        return summary
