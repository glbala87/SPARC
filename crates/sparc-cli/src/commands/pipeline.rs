//! Run full analysis pipeline

use anyhow::Result;
use clap::Args;
use std::path::PathBuf;

#[derive(Args)]
pub struct PipelineArgs {
    /// Input R1 FASTQ file
    #[arg(short = '1', long)]
    r1: PathBuf,

    /// Input R2 FASTQ file
    #[arg(short = '2', long)]
    r2: PathBuf,

    /// Reference genome directory
    #[arg(short = 'r', long)]
    reference: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Barcode whitelist file
    #[arg(short = 'w', long)]
    whitelist: PathBuf,

    /// Protocol
    #[arg(short, long, default_value = "10x-3prime-v3")]
    protocol: String,

    /// Sample name
    #[arg(short, long, default_value = "sample")]
    sample: String,

    /// Expected number of cells
    #[arg(long)]
    expect_cells: Option<u32>,

    /// Force number of cells
    #[arg(long)]
    force_cells: Option<u32>,
}

pub fn run(args: PipelineArgs) -> Result<()> {
    println!("=== sparc Pipeline ===\n");

    println!("Configuration:");
    println!("  R1 FASTQ:   {:?}", args.r1);
    println!("  R2 FASTQ:   {:?}", args.r2);
    println!("  Reference:  {:?}", args.reference);
    println!("  Output:     {:?}", args.output);
    println!("  Whitelist:  {:?}", args.whitelist);
    println!("  Protocol:   {}", args.protocol);
    println!("  Sample:     {}", args.sample);

    if let Some(n) = args.expect_cells {
        println!("  Expected cells: {}", n);
    }
    if let Some(n) = args.force_cells {
        println!("  Force cells: {}", n);
    }

    // Create output directory
    std::fs::create_dir_all(&args.output)?;

    println!("\n--- Step 1: Extract barcodes and UMIs ---");
    println!("(Not yet implemented - run 'sparc extract' separately)\n");

    println!("--- Step 2: Align reads ---");
    println!("(Requires external aligner - STAR or minimap2)\n");

    println!("--- Step 3: Generate count matrix ---");
    println!("(Not yet implemented - run 'sparc count' separately)\n");

    println!("--- Step 4: QC and filtering ---");
    println!("(Not yet implemented - run 'sparc qc' separately)\n");

    println!("Pipeline execution not yet fully implemented.");
    println!("Please run individual commands:");
    println!("  1. sparc extract -1 R1.fq.gz -2 R2.fq.gz -w whitelist.txt -o extracted/");
    println!("  2. STAR --genomeDir ref/ --readFilesIn extracted/*.fq.gz ...");
    println!("  3. sparc count -i aligned.bam -o counts/");
    println!("  4. sparc qc -i counts/ -o qc_report.json");

    Ok(())
}
