//! Extract barcodes and UMIs from FASTQ files

use anyhow::{Context, Result};
use clap::Args;
use indicatif::{ProgressBar, ProgressStyle};
use sparc_core::{
    barcode::{BarcodeCorrector, BarcodeMatch, Whitelist},
    fastq::FastqParser,
    protocols::{Protocol, TenX3Prime, TenX5Prime},
};
use std::path::PathBuf;

#[derive(Args)]
pub struct ExtractArgs {
    /// Input R1 FASTQ file (barcode/UMI read)
    #[arg(short = '1', long)]
    r1: PathBuf,

    /// Input R2 FASTQ file (cDNA read)
    #[arg(short = '2', long)]
    r2: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Barcode whitelist file
    #[arg(short = 'w', long)]
    whitelist: PathBuf,

    /// Protocol (10x-3prime-v3, 10x-3prime-v2, 10x-5prime-v2)
    #[arg(short, long, default_value = "10x-3prime-v3")]
    protocol: String,

    /// Maximum Hamming distance for barcode correction
    #[arg(long, default_value = "1")]
    max_mismatch: u32,

    /// Minimum barcode quality score
    #[arg(long, default_value = "10")]
    min_barcode_qual: u8,
}

pub fn run(args: ExtractArgs) -> Result<()> {
    log::info!("Loading barcode whitelist from {:?}", args.whitelist);
    let whitelist = Whitelist::from_file(&args.whitelist)
        .context("Failed to load barcode whitelist")?;
    log::info!("Loaded {} barcodes", whitelist.len());

    let corrector = BarcodeCorrector::new(whitelist, args.max_mismatch);

    let protocol: Box<dyn Protocol> = match args.protocol.as_str() {
        "10x-3prime-v3" => Box::new(TenX3Prime::v3()),
        "10x-3prime-v2" => Box::new(TenX3Prime::v2()),
        "10x-5prime-v2" => Box::new(TenX5Prime::v2()),
        _ => anyhow::bail!("Unknown protocol: {}", args.protocol),
    };

    log::info!("Using protocol: {} {}", protocol.name(), protocol.version());

    // Create output directory
    std::fs::create_dir_all(&args.output)?;

    // Open input files
    let mut r1_parser = FastqParser::open(&args.r1)
        .context("Failed to open R1 FASTQ")?;

    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );

    let mut total_reads = 0u64;
    let mut valid_barcode = 0u64;
    let mut corrected_barcode = 0u64;

    // Process reads
    for result in &mut r1_parser {
        let record = result?;
        total_reads += 1;

        if total_reads % 100000 == 0 {
            progress.set_message(format!(
                "Processed {} reads, {} valid barcodes ({:.1}%)",
                total_reads,
                valid_barcode,
                valid_barcode as f64 / total_reads as f64 * 100.0
            ));
        }

        // Extract barcode and UMI
        let components = match protocol.extract_r1(&record.seq, &record.qual) {
            Ok(c) => c,
            Err(_) => continue,
        };

        // Check barcode quality
        if !components.barcode_quality_ok(args.min_barcode_qual) {
            continue;
        }

        // Match barcode
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

    progress.finish_with_message(format!(
        "Done! Processed {} reads",
        total_reads
    ));

    // Print summary
    println!("\n=== Extraction Summary ===");
    println!("Total reads:        {}", total_reads);
    println!("Valid barcodes:     {} ({:.1}%)",
        valid_barcode,
        valid_barcode as f64 / total_reads as f64 * 100.0
    );
    println!("Corrected barcodes: {} ({:.1}%)",
        corrected_barcode,
        corrected_barcode as f64 / total_reads as f64 * 100.0
    );

    Ok(())
}
