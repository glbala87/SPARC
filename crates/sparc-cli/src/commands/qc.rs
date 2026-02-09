//! Generate QC report

use anyhow::{Context, Result};
use clap::Args;
use sparc_core::qc::{CellMetrics, QcMetrics, QcReport};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

#[derive(Args)]
pub struct QcArgs {
    /// Input matrix directory (with matrix.mtx, barcodes.tsv, genes.tsv)
    #[arg(short, long)]
    input: PathBuf,

    /// Output QC report file (JSON)
    #[arg(short, long)]
    output: PathBuf,

    /// Sample name
    #[arg(short, long, default_value = "sample")]
    sample: String,

    /// Minimum genes per cell
    #[arg(long, default_value = "200")]
    min_genes: u64,

    /// Maximum genes per cell
    #[arg(long, default_value = "10000")]
    max_genes: u64,

    /// Maximum mitochondrial percentage
    #[arg(long, default_value = "20.0")]
    max_mito: f64,
}

pub fn run(args: QcArgs) -> Result<()> {
    log::info!("Reading count matrix from {:?}", args.input);

    // Read barcodes
    let barcodes_path = args.input.join("barcodes.tsv");
    let barcodes: Vec<String> = BufReader::new(File::open(&barcodes_path)?)
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?;

    // Read genes
    let genes_path = args.input.join("genes.tsv");
    let genes: Vec<String> = BufReader::new(File::open(&genes_path)?)
        .lines()
        .map(|l| l.map(|s| s.split('\t').next().unwrap_or("").to_string()))
        .collect::<std::io::Result<Vec<_>>>()?;

    // Read matrix
    let mtx_path = args.input.join("matrix.mtx");
    let mtx_file = BufReader::new(File::open(&mtx_path)?);

    let mut n_rows = 0;
    let mut n_cols = 0;
    let mut counts_per_cell: Vec<u64> = Vec::new();
    let mut genes_per_cell: Vec<ahash::AHashSet<usize>> = Vec::new();

    for (i, line) in mtx_file.lines().enumerate() {
        let line = line?;
        if line.starts_with('%') {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 3 && i <= 2 {
            // Header line
            n_rows = parts[0].parse().unwrap_or(0);
            n_cols = parts[1].parse().unwrap_or(0);
            counts_per_cell = vec![0; n_cols];
            genes_per_cell = (0..n_cols).map(|_| ahash::AHashSet::new()).collect();
        } else if parts.len() == 3 {
            // Data line (1-indexed)
            let row: usize = parts[0].parse().context("Invalid row index")?;
            let col: usize = parts[1].parse().context("Invalid col index")?;
            let val: u64 = parts[2].parse().context("Invalid value")?;

            if col > 0 && col <= n_cols && row > 0 && row <= n_rows {
                counts_per_cell[col - 1] += val;
                genes_per_cell[col - 1].insert(row - 1);
            }
        }
    }

    let genes_per_cell_count: Vec<u64> = genes_per_cell.iter().map(|s| s.len() as u64).collect();

    // Calculate metrics
    let mut metrics = QcMetrics::new();
    metrics.num_cells = n_cols as u64;
    metrics.total_genes = n_rows as u64;
    metrics.update_from_cells(&counts_per_cell, &genes_per_cell_count, &counts_per_cell);

    // Build report
    let mut report = QcReport::new(args.sample.clone());
    report.metrics = metrics;

    // Per-cell metrics
    for (i, barcode) in barcodes.iter().enumerate() {
        let cell_metrics = CellMetrics {
            barcode: barcode.clone(),
            reads: counts_per_cell.get(i).copied().unwrap_or(0),
            genes: genes_per_cell_count.get(i).copied().unwrap_or(0),
            umis: counts_per_cell.get(i).copied().unwrap_or(0),
            mito_percent: 0.0, // Would need gene annotations
        };
        report.per_cell_metrics.push(cell_metrics);
    }

    // Generate warnings
    report.generate_warnings();

    // Apply filters
    let filtered_cells: Vec<_> = report
        .per_cell_metrics
        .iter()
        .filter(|c| c.genes >= args.min_genes && c.genes <= args.max_genes)
        .count();

    // Write report
    let json = report.to_json()?;
    std::fs::write(&args.output, &json)?;

    // Print summary
    println!("\n=== QC Summary ===");
    println!("Sample:              {}", args.sample);
    println!("Total cells:         {}", n_cols);
    println!("Total genes:         {}", n_rows);
    println!("Median genes/cell:   {:.0}", report.metrics.median_genes_per_cell);
    println!("Median UMIs/cell:    {:.0}", report.metrics.median_umi_per_cell);
    println!("Cells passing QC:    {} ({:.1}%)",
        filtered_cells,
        filtered_cells as f64 / n_cols as f64 * 100.0
    );

    if !report.warnings.is_empty() {
        println!("\nWarnings:");
        for warning in &report.warnings {
            println!("  - {}", warning);
        }
    }

    println!("\nQC report written to {:?}", args.output);

    Ok(())
}
