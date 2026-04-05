//! Multi-sample batch processing

use anyhow::{Context, Result};
use clap::Args;
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

#[derive(Args)]
pub struct BatchArgs {
    /// Sample manifest CSV (columns: sample_name,r1,r2,whitelist)
    #[arg(short, long)]
    manifest: PathBuf,

    /// Reference genome directory
    #[arg(short = 'r', long)]
    reference: PathBuf,

    /// Output directory
    #[arg(short, long)]
    output: PathBuf,

    /// Protocol
    #[arg(short, long, default_value = "10x-3prime-v3")]
    protocol: String,

    /// Aligner
    #[arg(long, default_value = "star")]
    aligner: String,

    /// Maximum barcode mismatch
    #[arg(long, default_value = "1")]
    max_mismatch: u32,

    /// Parallel samples
    #[arg(long, default_value = "1")]
    parallel_samples: usize,
}

#[derive(Debug, Clone)]
struct SampleConfig {
    name: String,
    r1: PathBuf,
    r2: PathBuf,
    whitelist: PathBuf,
}

fn parse_manifest(path: &PathBuf) -> Result<Vec<SampleConfig>> {
    let file = File::open(path).context("Failed to open manifest")?;
    let reader = BufReader::new(file);
    let mut samples = Vec::new();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();

        if i == 0 && line.to_lowercase().contains("sample") {
            continue;
        }
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 4 {
            anyhow::bail!(
                "Line {}: expected 4 columns (sample_name,r1,r2,whitelist), got {}",
                i + 1,
                parts.len()
            );
        }

        samples.push(SampleConfig {
            name: parts[0].trim().to_string(),
            r1: PathBuf::from(parts[1].trim()),
            r2: PathBuf::from(parts[2].trim()),
            whitelist: PathBuf::from(parts[3].trim()),
        });
    }

    Ok(samples)
}

pub fn run(args: BatchArgs) -> Result<()> {
    println!("=== SPARC Batch Processing ===\n");

    let samples = parse_manifest(&args.manifest)?;
    println!("Found {} samples in manifest\n", samples.len());

    std::fs::create_dir_all(&args.output)?;

    let results: Vec<(String, Result<()>)> = if args.parallel_samples > 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(args.parallel_samples)
            .build()
            .context("Failed to build thread pool")?;

        pool.install(|| {
            samples
                .par_iter()
                .map(|sample| {
                    let result = process_sample(sample, &args);
                    (sample.name.clone(), result)
                })
                .collect()
        })
    } else {
        samples
            .iter()
            .map(|sample| {
                println!("Processing sample: {}", sample.name);
                let result = process_sample(sample, &args);
                (sample.name.clone(), result)
            })
            .collect()
    };

    println!("\n=== Batch Summary ===");
    let mut succeeded = 0;
    let mut failed = 0;

    for (name, result) in &results {
        match result {
            Ok(()) => {
                println!("  [OK]   {}", name);
                succeeded += 1;
            }
            Err(e) => {
                println!("  [FAIL] {}: {}", name, e);
                failed += 1;
            }
        }
    }

    println!("\nTotal: {} succeeded, {} failed", succeeded, failed);

    if failed > 0 {
        anyhow::bail!("{} samples failed", failed);
    }

    Ok(())
}

fn process_sample(sample: &SampleConfig, args: &BatchArgs) -> Result<()> {
    let sample_output = args.output.join(&sample.name);

    let pipeline_args = super::pipeline::PipelineArgs {
        r1: sample.r1.clone(),
        r2: sample.r2.clone(),
        reference: args.reference.clone(),
        output: sample_output,
        whitelist: sample.whitelist.clone(),
        protocol: args.protocol.clone(),
        sample: sample.name.clone(),
        aligner: args.aligner.clone(),
        max_mismatch: args.max_mismatch,
        min_barcode_qual: 10,
        min_mapq: 30,
        expect_cells: None,
        force_cells: None,
        skip_align: false,
        bam: None,
        min_genes: 200,
        max_genes: 10000,
    };

    super::pipeline::run(pipeline_args)
}
