//! Generate gene count matrix from BAM file

use anyhow::{Context, Result};
use clap::Args;
use indicatif::{ProgressBar, ProgressStyle};
use sparc_core::{
    bam::BamParser,
    count::GeneCounter,
};
use std::path::PathBuf;

#[derive(Args)]
pub struct CountArgs {
    /// Input BAM file (with CB, UB, GN/GX tags)
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory for matrix files
    #[arg(short, long)]
    output: PathBuf,

    /// Minimum mapping quality
    #[arg(long, default_value = "30")]
    min_mapq: u8,

    /// Output format (mtx, h5ad)
    #[arg(long, default_value = "mtx")]
    format: String,
}

pub fn run(args: CountArgs) -> Result<()> {
    log::info!("Opening BAM file: {:?}", args.input);
    let mut parser = BamParser::open(&args.input)
        .context("Failed to open BAM file")?;

    // Create output directory
    std::fs::create_dir_all(&args.output)?;

    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    let mut counter = GeneCounter::new();
    let mut total_reads = 0u64;
    let mut assigned_reads = 0u64;

    // Process BAM records
    for result in &mut parser {
        let record = result?;
        total_reads += 1;

        if total_reads % 100000 == 0 {
            progress.set_message(format!(
                "Processed {} reads, {} assigned ({:.1}%)",
                total_reads,
                assigned_reads,
                assigned_reads as f64 / total_reads as f64 * 100.0
            ));
        }

        // Skip unmapped or low quality
        if !record.is_mapped || record.mapq < args.min_mapq {
            continue;
        }

        // Need cell barcode and gene
        let (barcode, gene) = match (&record.cell_barcode, &record.gene_name) {
            (Some(bc), Some(gn)) => (bc, gn),
            (Some(bc), None) => {
                // Try gene_id if gene_name not available
                if let Some(gx) = &record.gene_id {
                    (bc, gx)
                } else {
                    continue;
                }
            }
            _ => continue,
        };

        counter.increment(barcode, gene);
        assigned_reads += 1;
    }

    progress.finish_with_message(format!(
        "Done! Processed {} reads",
        total_reads
    ));

    // Build matrix
    log::info!("Building count matrix...");
    let matrix = counter.build();

    log::info!("Matrix dimensions: {} genes x {} cells",
        matrix.n_rows, matrix.n_cols);
    log::info!("Non-zero entries: {}", matrix.values.len());

    // Write output
    match args.format.as_str() {
        "mtx" => {
            let mtx_path = args.output.join("matrix.mtx");
            let barcodes_path = args.output.join("barcodes.tsv");
            let genes_path = args.output.join("genes.tsv");

            log::info!("Writing Matrix Market files...");
            matrix.write_mtx(&mtx_path)?;
            matrix.write_barcodes(&barcodes_path)?;
            matrix.write_genes(&genes_path)?;

            println!("\nOutput files:");
            println!("  {:?}", mtx_path);
            println!("  {:?}", barcodes_path);
            println!("  {:?}", genes_path);
        }
        "h5ad" => {
            anyhow::bail!("H5AD format not yet implemented");
        }
        _ => anyhow::bail!("Unknown format: {}", args.format),
    }

    // Print summary
    println!("\n=== Count Summary ===");
    println!("Total reads:    {}", total_reads);
    println!("Assigned reads: {} ({:.1}%)",
        assigned_reads,
        assigned_reads as f64 / total_reads as f64 * 100.0
    );
    println!("Cells:          {}", matrix.n_cols);
    println!("Genes:          {}", matrix.n_rows);

    Ok(())
}
