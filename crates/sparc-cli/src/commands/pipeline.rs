//! Run full analysis pipeline

use anyhow::{Context, Result};
use clap::Args;
use indicatif::{ProgressBar, ProgressStyle};
use sparc_core::{
    aligner::{Aligner, AlignerConfig},
    bam::BamParser,
    barcode::{BarcodeCorrector, BarcodeMatch, Whitelist},
    count::GeneCounter,
    fastq::FastqParser,
    protocols::{DropSeq, InDrop, Protocol, SciRNA, SmartSeq2, TenX3Prime, TenX5Prime},
    qc::{CellMetrics, QcMetrics, QcReport},
};
use std::path::PathBuf;

#[derive(Args)]
pub struct PipelineArgs {
    /// Input R1 FASTQ file
    #[arg(short = '1', long)]
    pub(crate) r1: PathBuf,

    /// Input R2 FASTQ file
    #[arg(short = '2', long)]
    pub(crate) r2: PathBuf,

    /// Reference genome directory
    #[arg(short = 'r', long)]
    pub(crate) reference: PathBuf,

    /// Output directory
    #[arg(short, long)]
    pub(crate) output: PathBuf,

    /// Barcode whitelist file
    #[arg(short = 'w', long)]
    pub(crate) whitelist: PathBuf,

    /// Protocol (10x-3prime-v3, 10x-3prime-v2, 10x-5prime-v2, drop-seq, indrop, sci-rna-seq, smart-seq2)
    #[arg(short, long, default_value = "10x-3prime-v3")]
    pub(crate) protocol: String,

    /// Sample name
    #[arg(short, long, default_value = "sample")]
    pub(crate) sample: String,

    /// Aligner (star, minimap2)
    #[arg(long, default_value = "star")]
    pub(crate) aligner: String,

    /// Maximum Hamming distance for barcode correction
    #[arg(long, default_value = "1")]
    pub(crate) max_mismatch: u32,

    /// Minimum barcode quality score
    #[arg(long, default_value = "10")]
    pub(crate) min_barcode_qual: u8,

    /// Minimum mapping quality
    #[arg(long, default_value = "30")]
    pub(crate) min_mapq: u8,

    /// Expected number of cells
    #[arg(long)]
    pub(crate) expect_cells: Option<u32>,

    /// Force number of cells
    #[arg(long)]
    pub(crate) force_cells: Option<u32>,

    /// Skip alignment step
    #[arg(long)]
    pub(crate) skip_align: bool,

    /// Pre-aligned BAM file (with --skip-align)
    #[arg(long)]
    pub(crate) bam: Option<PathBuf>,

    /// Minimum genes per cell for QC
    #[arg(long, default_value = "200")]
    pub(crate) min_genes: u64,

    /// Maximum genes per cell for QC
    #[arg(long, default_value = "10000")]
    pub(crate) max_genes: u64,
}

pub(crate) fn get_protocol(name: &str) -> Result<Box<dyn Protocol>> {
    match name {
        "10x-3prime-v3" => Ok(Box::new(TenX3Prime::v3())),
        "10x-3prime-v2" => Ok(Box::new(TenX3Prime::v2())),
        "10x-5prime-v2" => Ok(Box::new(TenX5Prime::v2())),
        "drop-seq" => Ok(Box::new(DropSeq::new())),
        "indrop" => Ok(Box::new(InDrop::new())),
        "sci-rna-seq" => Ok(Box::new(SciRNA::new())),
        "smart-seq2" => Ok(Box::new(SmartSeq2::new("sample".to_string()))),
        _ => anyhow::bail!("Unknown protocol: {}", name),
    }
}

pub fn run(args: PipelineArgs) -> Result<()> {
    println!("=== SPARC Pipeline ===\n");

    // Create output directories
    let extract_dir = args.output.join("extraction");
    let align_dir = args.output.join("alignment");
    let count_dir = args.output.join("counts");
    let qc_dir = args.output.join("qc");

    std::fs::create_dir_all(&extract_dir)?;
    std::fs::create_dir_all(&align_dir)?;
    std::fs::create_dir_all(&count_dir)?;
    std::fs::create_dir_all(&qc_dir)?;

    // ===== Step 1: Extract barcodes =====
    println!("--- Step 1/4: Extracting barcodes and UMIs ---");

    let whitelist =
        Whitelist::from_file(&args.whitelist).context("Failed to load barcode whitelist")?;
    log::info!("Loaded {} barcodes", whitelist.len());

    let corrector = BarcodeCorrector::new(whitelist, args.max_mismatch);
    let protocol = get_protocol(&args.protocol)?;

    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .expect("valid progress template"),
    );

    let mut r1_parser =
        FastqParser::open(&args.r1).context("Failed to open R1 FASTQ")?;

    let mut total_reads = 0u64;
    let mut valid_barcode = 0u64;
    let mut corrected_barcode = 0u64;

    for result in &mut r1_parser {
        let record = result?;
        total_reads += 1;

        if total_reads % 100000 == 0 {
            pb.set_message(format!(
                "Processed {} reads, {} valid barcodes ({:.1}%)",
                total_reads,
                valid_barcode,
                valid_barcode as f64 / total_reads as f64 * 100.0
            ));
        }

        let components = match protocol.extract_r1(&record.seq, &record.qual) {
            Ok(c) => c,
            Err(_) => continue,
        };

        if !components.barcode_quality_ok(args.min_barcode_qual) {
            continue;
        }

        let barcode_str = components.barcode_str();
        match corrector.match_barcode(&barcode_str) {
            BarcodeMatch::Exact(_) => {
                valid_barcode += 1;
            }
            BarcodeMatch::Corrected(_, _, _) => {
                valid_barcode += 1;
                corrected_barcode += 1;
            }
            BarcodeMatch::NoMatch(_) => {}
        }
    }

    pb.finish_with_message(format!(
        "Extraction complete: {} valid from {} reads",
        valid_barcode, total_reads
    ));

    println!("  Total reads:        {}", total_reads);
    println!(
        "  Valid barcodes:     {} ({:.1}%)",
        valid_barcode,
        valid_barcode as f64 / total_reads.max(1) as f64 * 100.0
    );
    println!(
        "  Corrected barcodes: {} ({:.1}%)",
        corrected_barcode,
        corrected_barcode as f64 / total_reads.max(1) as f64 * 100.0
    );

    // ===== Step 2: Alignment =====
    let bam_path = if args.skip_align {
        println!("\n--- Step 2/4: Alignment (skipped) ---");
        args.bam.clone().unwrap_or_else(|| args.output.join("aligned.bam"))
    } else {
        println!("\n--- Step 2/4: Aligning reads ---");

        let aligner_config = match args.aligner.as_str() {
            "star" => AlignerConfig::star(args.reference.clone(), rayon::current_num_threads()),
            "minimap2" => {
                AlignerConfig::minimap2(args.reference.clone(), rayon::current_num_threads())
            }
            _ => anyhow::bail!("Unknown aligner: {}", args.aligner),
        };

        let aligner = Aligner::new(aligner_config);

        if !aligner.is_available() {
            println!("  WARNING: {} not found in PATH", aligner.binary_name());
            println!("  Skipping alignment. Use --skip-align --bam <file> with a pre-aligned BAM.");
            args.output.join("aligned.bam")
        } else {
            let bam = aligner
                .align(&args.r2, Some(&args.r1), &align_dir)
                .context("Alignment failed")?;
            println!("  Alignment complete: {:?}", bam);
            bam
        }
    };

    // ===== Step 3: Count matrix =====
    println!("\n--- Step 3/4: Generating count matrix ---");

    if bam_path.exists() {
        let pb3 = ProgressBar::new_spinner();
        pb3.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .expect("valid progress template"),
        );

        let mut bam_parser =
            BamParser::open(&bam_path).context("Failed to open BAM file")?;

        let mut counter = GeneCounter::new();
        let mut bam_total = 0u64;
        let mut assigned = 0u64;

        for result in &mut bam_parser {
            let record = result?;
            bam_total += 1;

            if bam_total % 100000 == 0 {
                pb3.set_message(format!("Processed {} reads, {} assigned", bam_total, assigned));
            }

            if !record.is_mapped || record.mapq < args.min_mapq {
                continue;
            }

            let (barcode, gene) = match (&record.cell_barcode, &record.gene_name) {
                (Some(bc), Some(gn)) => (bc, gn),
                (Some(bc), None) => {
                    if let Some(gx) = &record.gene_id {
                        (bc, gx)
                    } else {
                        continue;
                    }
                }
                _ => continue,
            };

            counter.increment(barcode, gene);
            assigned += 1;
        }

        pb3.finish_with_message(format!(
            "Count matrix: {} assigned from {} reads",
            assigned, bam_total
        ));

        let matrix = counter.build();

        matrix.write_mtx(count_dir.join("matrix.mtx"))?;
        matrix.write_barcodes(count_dir.join("barcodes.tsv"))?;
        matrix.write_genes(count_dir.join("genes.tsv"))?;

        println!(
            "  Matrix: {} genes x {} cells",
            matrix.n_rows, matrix.n_cols
        );
        println!("  Non-zero entries: {}", matrix.values.len());

        // ===== Step 4: QC =====
        println!("\n--- Step 4/4: Quality control ---");

        let counts_per_cell = matrix.counts_per_cell();
        let genes_per_cell = matrix.genes_per_cell();

        let mut metrics = QcMetrics::new();
        metrics.total_reads = total_reads;
        metrics.valid_barcode_reads = valid_barcode;
        metrics.mapped_reads = bam_total;
        metrics.assigned_reads = assigned;
        metrics.num_cells = matrix.n_cols as u64;
        metrics.total_genes = matrix.n_rows as u64;
        metrics.update_from_cells(&counts_per_cell, &genes_per_cell, &counts_per_cell);

        let mut report = QcReport::new(args.sample.clone());
        report.metrics = metrics;

        for (i, barcode) in matrix.barcodes.iter().enumerate() {
            report.per_cell_metrics.push(CellMetrics {
                barcode: barcode.clone(),
                reads: counts_per_cell.get(i).copied().unwrap_or(0),
                genes: genes_per_cell.get(i).copied().unwrap_or(0),
                umis: counts_per_cell.get(i).copied().unwrap_or(0),
                mito_percent: 0.0,
            });
        }

        report.generate_warnings();

        let filtered_cells = report
            .per_cell_metrics
            .iter()
            .filter(|c| c.genes >= args.min_genes && c.genes <= args.max_genes)
            .count();

        let json = report.to_json()?;
        std::fs::write(qc_dir.join("qc_report.json"), &json)?;

        println!(
            "  Median genes/cell: {:.0}",
            report.metrics.median_genes_per_cell
        );
        println!(
            "  Median UMIs/cell:  {:.0}",
            report.metrics.median_umi_per_cell
        );
        println!(
            "  Cells passing QC:  {} ({:.1}%)",
            filtered_cells,
            filtered_cells as f64 / matrix.n_cols.max(1) as f64 * 100.0
        );

        if !report.warnings.is_empty() {
            println!("\n  Warnings:");
            for w in &report.warnings {
                println!("    - {}", w);
            }
        }
    } else {
        println!(
            "  BAM file not found at {:?}, skipping count matrix and QC",
            bam_path
        );
    }

    println!("\n=== Pipeline Complete ===");
    println!("Output directory: {:?}", args.output);

    Ok(())
}
