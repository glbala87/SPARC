//! Distributed processing support

use anyhow::{Context, Result};
use clap::Args;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

#[derive(Args)]
pub struct DistributedArgs {
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

    /// Number of shards
    #[arg(long, default_value = "4")]
    shards: usize,

    /// Current shard index (0-based, worker mode)
    #[arg(long)]
    shard_index: Option<usize>,

    /// Worker mode: process only one shard
    #[arg(long)]
    worker: bool,

    /// Merge mode: combine shard results
    #[arg(long)]
    merge: bool,
}

pub fn run(args: DistributedArgs) -> Result<()> {
    if args.merge {
        return merge_shards(&args);
    }

    if args.worker {
        let shard_idx = args
            .shard_index
            .context("--shard-index required in worker mode")?;
        return process_shard(&args, shard_idx);
    }

    // Coordinator mode
    println!("=== SPARC Distributed Processing ===\n");
    println!("Splitting work into {} shards\n", args.shards);

    std::fs::create_dir_all(args.output.join("shards"))?;

    println!("Run these commands on separate machines:\n");
    for i in 0..args.shards {
        println!("  sparc distributed \\");
        println!("    -1 {:?} -2 {:?} \\", args.r1, args.r2);
        println!("    -r {:?} -o {:?} \\", args.reference, args.output);
        println!("    -w {:?} \\", args.whitelist);
        println!(
            "    --protocol {} --shards {} --shard-index {} --worker\n",
            args.protocol, args.shards, i
        );
    }

    println!("After all shards complete, merge with:\n");
    println!("  sparc distributed \\");
    println!("    -1 {:?} -2 {:?} \\", args.r1, args.r2);
    println!("    -r {:?} -o {:?} \\", args.reference, args.output);
    println!("    -w {:?} \\", args.whitelist);
    println!("    --shards {} --merge\n", args.shards);

    Ok(())
}

fn process_shard(args: &DistributedArgs, shard_index: usize) -> Result<()> {
    println!(
        "=== Processing shard {}/{} ===\n",
        shard_index + 1,
        args.shards
    );

    let shard_output = args.output.join(format!("shard_{}", shard_index));
    std::fs::create_dir_all(&shard_output)?;

    use sparc_core::{
        barcode::{BarcodeCorrector, Whitelist},
        count::GeneCounter,
        fastq::FastqParser,
    };

    let whitelist = Whitelist::from_file(&args.whitelist)?;
    let corrector = BarcodeCorrector::new(whitelist, 1);
    let protocol = super::pipeline::get_protocol(&args.protocol)?;

    let mut parser = FastqParser::open(&args.r1)?;
    let mut counter = GeneCounter::new();
    let mut record_num = 0u64;
    let mut processed = 0u64;

    for result in &mut parser {
        let record = result?;

        if (record_num as usize % args.shards) == shard_index {
            processed += 1;

            if let Ok(components) = protocol.extract_r1(&record.seq, &record.qual) {
                if components.barcode_quality_ok(10) {
                    let barcode_str = components.barcode_str();
                    if let Some(bc) = corrector.match_barcode(&barcode_str).barcode() {
                        counter.increment(bc, &components.umi_str());
                    }
                }
            }
        }

        record_num += 1;
    }

    let matrix = counter.build();
    matrix.write_mtx(shard_output.join("matrix.mtx"))?;
    matrix.write_barcodes(shard_output.join("barcodes.tsv"))?;
    matrix.write_genes(shard_output.join("genes.tsv"))?;

    println!(
        "Shard {} complete: processed {} of {} reads",
        shard_index, processed, record_num
    );
    println!("Matrix: {} genes x {} cells", matrix.n_rows, matrix.n_cols);

    Ok(())
}

fn merge_shards(args: &DistributedArgs) -> Result<()> {
    use sparc_core::count::GeneCounter;

    println!("=== Merging {} shards ===\n", args.shards);

    let mut counter = GeneCounter::new();
    let mut shards_found = 0usize;
    let mut shards_missing = Vec::new();

    for i in 0..args.shards {
        let shard_dir = args.output.join(format!("shard_{}", i));
        let mtx_path = shard_dir.join("matrix.mtx");
        let barcodes_path = shard_dir.join("barcodes.tsv");
        let genes_path = shard_dir.join("genes.tsv");

        if !mtx_path.exists() {
            println!("  Warning: shard {} not found at {:?}, skipping", i, shard_dir);
            shards_missing.push(i);
            continue;
        }

        if !barcodes_path.exists() {
            anyhow::bail!("Shard {} is incomplete: barcodes.tsv missing at {:?}", i, barcodes_path);
        }
        if !genes_path.exists() {
            anyhow::bail!("Shard {} is incomplete: genes.tsv missing at {:?}", i, genes_path);
        }

        let barcodes: Vec<String> = BufReader::new(File::open(&barcodes_path)?)
            .lines()
            .collect::<std::io::Result<_>>()?;
        let genes: Vec<String> = BufReader::new(File::open(&genes_path)?)
            .lines()
            .map(|l| {
                l.map(|s| {
                    s.split('\t')
                        .next()
                        .unwrap_or("")
                        .to_string()
                })
            })
            .filter(|r| r.as_ref().map_or(true, |s| !s.is_empty()))
            .collect::<std::io::Result<_>>()?;

        if barcodes.is_empty() && genes.is_empty() {
            println!("  Shard {} is empty, skipping", i);
            shards_found += 1;
            continue;
        }

        let mut header_skipped = false;
        let mtx_file = BufReader::new(File::open(&mtx_path)?);
        for line in mtx_file.lines() {
            let line = line?;
            if line.starts_with('%') {
                continue;
            }
            // Skip the dimensions line (first non-comment line)
            if !header_skipped {
                header_skipped = true;
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() != 3 {
                continue;
            }
            let row = match parts[0].parse::<usize>() {
                Ok(r) if r >= 1 => r,
                _ => continue,
            };
            let col = match parts[1].parse::<usize>() {
                Ok(c) if c >= 1 => c,
                _ => continue,
            };
            let val = match parts[2].parse::<u32>() {
                Ok(v) => v,
                _ => continue,
            };

            // Bounds check (1-indexed MTX format)
            if row > genes.len() || col > barcodes.len() {
                log::warn!(
                    "Shard {}: out-of-bounds entry row={} col={} (genes={}, barcodes={})",
                    i, row, col, genes.len(), barcodes.len()
                );
                continue;
            }

            let gene = &genes[row - 1];
            let barcode = &barcodes[col - 1];
            if !gene.is_empty() && !barcode.is_empty() {
                counter.add_count(barcode, gene, val);
            }
        }

        shards_found += 1;
        println!("  Merged shard {}", i);
    }

    if shards_found == 0 {
        anyhow::bail!(
            "No shard data found. Missing shards: {:?}. Ensure all workers completed before merging.",
            shards_missing
        );
    }

    if !shards_missing.is_empty() {
        println!(
            "\n  WARNING: {} of {} shards missing: {:?}",
            shards_missing.len(),
            args.shards,
            shards_missing
        );
        println!("  Merged result is PARTIAL. Re-run missing shards and merge again.");
    }

    let merged_dir = args.output.join("merged");
    std::fs::create_dir_all(&merged_dir)?;

    let matrix = counter.build();
    matrix.write_mtx(merged_dir.join("matrix.mtx"))?;
    matrix.write_barcodes(merged_dir.join("barcodes.tsv"))?;
    matrix.write_genes(merged_dir.join("genes.tsv"))?;

    println!(
        "\nMerged matrix: {} genes x {} cells ({}/{} shards)",
        matrix.n_rows, matrix.n_cols, shards_found, args.shards
    );
    println!("Output: {:?}", merged_dir);

    Ok(())
}
