//! Run downstream analysis on a count matrix (normalize, PCA, clustering)

use anyhow::{Context, Result};
use clap::Args;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use sparc_core::analysis::{
    cluster::label_propagation,
    neighbors::build_knn_graph,
    normalize::{highly_variable_genes, log1p_transform, normalize_total, scale},
    pca::pca,
};

#[derive(Args)]
pub struct AnalyzeArgs {
    /// Input count matrix directory (with matrix.mtx, barcodes.tsv, genes.tsv)
    #[arg(short, long)]
    input: PathBuf,

    /// Output directory for analysis results
    #[arg(short, long)]
    output: PathBuf,

    /// Number of principal components
    #[arg(long, default_value = "50")]
    n_pcs: usize,

    /// Number of nearest neighbors for graph construction
    #[arg(long, default_value = "15")]
    n_neighbors: usize,

    /// Clustering resolution
    #[arg(long, default_value = "1.0")]
    resolution: f64,

    /// Target sum for normalization
    #[arg(long, default_value = "10000")]
    target_sum: f64,

    /// Minimum mean for HVG selection
    #[arg(long, default_value = "0.0125")]
    min_mean: f64,

    /// Maximum mean for HVG selection
    #[arg(long, default_value = "3.0")]
    max_mean: f64,

    /// Minimum dispersion for HVG selection
    #[arg(long, default_value = "0.5")]
    min_disp: f64,
}

pub fn run(args: AnalyzeArgs) -> Result<()> {
    println!("=== SPARC Analysis ===\n");

    std::fs::create_dir_all(&args.output)?;

    // ─── Load count matrix ───────────────────────────────────────
    println!("Loading count matrix from {:?}...", args.input);

    let barcodes_path = args.input.join("barcodes.tsv");
    let genes_path = args.input.join("genes.tsv");
    let mtx_path = args.input.join("matrix.mtx");

    if !mtx_path.exists() {
        anyhow::bail!("matrix.mtx not found in {:?}", args.input);
    }
    if !barcodes_path.exists() {
        anyhow::bail!("barcodes.tsv not found in {:?}", args.input);
    }
    if !genes_path.exists() {
        anyhow::bail!("genes.tsv not found in {:?}", args.input);
    }

    let barcodes: Vec<String> = BufReader::new(File::open(&barcodes_path)?)
        .lines()
        .collect::<std::io::Result<_>>()?;
    let genes: Vec<String> = BufReader::new(File::open(&genes_path)?)
        .lines()
        .map(|l| l.map(|s| s.split('\t').next().unwrap_or("").to_string()))
        .collect::<std::io::Result<_>>()?;

    let n_cells = barcodes.len();
    let n_genes = genes.len();
    println!("  {} cells x {} genes", n_cells, n_genes);

    if n_cells == 0 || n_genes == 0 {
        anyhow::bail!("Empty matrix: {} cells, {} genes", n_cells, n_genes);
    }

    // Parse MTX to dense: cells x genes
    let mut data = vec![vec![0.0f64; n_genes]; n_cells];
    let mut header_skipped = false;
    for line in BufReader::new(File::open(&mtx_path)?).lines() {
        let line = line?;
        if line.starts_with('%') {
            continue;
        }
        if !header_skipped {
            header_skipped = true;
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() == 3 {
            if let (Ok(row), Ok(col), Ok(val)) = (
                parts[0].parse::<usize>(),
                parts[1].parse::<usize>(),
                parts[2].parse::<f64>(),
            ) {
                if row >= 1 && row <= n_genes && col >= 1 && col <= n_cells {
                    data[col - 1][row - 1] = val;
                }
            }
        }
    }

    // ─── Normalize ───────────────────────────────────────────────
    println!("Normalizing (target_sum={})...", args.target_sum);
    normalize_total(&mut data, args.target_sum);
    log1p_transform(&mut data);

    // ─── HVG selection ───────────────────────────────────────────
    println!(
        "Selecting highly variable genes (min_mean={}, max_mean={}, min_disp={})...",
        args.min_mean, args.max_mean, args.min_disp
    );
    let hvg_indices = highly_variable_genes(&data, args.min_mean, args.max_mean, args.min_disp);
    println!("  {} HVGs selected", hvg_indices.len());

    if hvg_indices.is_empty() {
        println!("  WARNING: No HVGs found. Using all genes.");
    }

    // Subset to HVGs if any
    let mut analysis_data = if hvg_indices.is_empty() {
        data.clone()
    } else {
        data.iter()
            .map(|cell| hvg_indices.iter().map(|&g| cell[g]).collect())
            .collect()
    };

    let hvg_names: Vec<&str> = if hvg_indices.is_empty() {
        genes.iter().map(|s| s.as_str()).collect()
    } else {
        hvg_indices.iter().map(|&i| genes[i].as_str()).collect()
    };

    // ─── Scale ───────────────────────────────────────────────────
    println!("Scaling...");
    scale(&mut analysis_data, Some(10.0));

    // ─── PCA ─────────────────────────────────────────────────────
    let n_pcs = args.n_pcs.min(analysis_data[0].len() - 1).min(n_cells - 1);
    println!("Computing PCA ({} components)...", n_pcs);
    let (pca_coords, variances) = pca(&analysis_data, n_pcs);
    println!(
        "  Variance explained (PC1): {:.4}, (PC2): {:.4}",
        variances.first().unwrap_or(&0.0),
        variances.get(1).unwrap_or(&0.0)
    );

    // ─── KNN graph ───────────────────────────────────────────────
    let k = args.n_neighbors.min(n_cells - 1);
    println!("Building KNN graph (k={})...", k);
    let graph = build_knn_graph(&pca_coords, k);

    // ─── Clustering ──────────────────────────────────────────────
    println!("Clustering (resolution={})...", args.resolution);
    let labels = label_propagation(&graph, args.resolution, 100);
    let n_clusters = labels.iter().copied().collect::<std::collections::HashSet<_>>().len();
    println!("  Found {} clusters", n_clusters);

    // ─── Save results ────────────────────────────────────────────
    println!("\nSaving results to {:?}...", args.output);

    // Cluster assignments
    {
        let mut f = File::create(args.output.join("clusters.tsv"))?;
        writeln!(f, "barcode\tcluster")?;
        for (i, bc) in barcodes.iter().enumerate() {
            writeln!(f, "{}\t{}", bc, labels[i])?;
        }
    }

    // PCA coordinates
    {
        let mut f = File::create(args.output.join("pca.tsv"))?;
        let header: Vec<String> = (1..=n_pcs).map(|i| format!("PC{}", i)).collect();
        writeln!(f, "barcode\t{}", header.join("\t"))?;
        for (i, bc) in barcodes.iter().enumerate() {
            let vals: Vec<String> = pca_coords[i].iter().map(|v| format!("{:.6}", v)).collect();
            writeln!(f, "{}\t{}", bc, vals.join("\t"))?;
        }
    }

    // HVG list
    {
        let mut f = File::create(args.output.join("hvg.tsv"))?;
        writeln!(f, "gene")?;
        for gene in &hvg_names {
            writeln!(f, "{}", gene)?;
        }
    }

    // Summary JSON
    {
        let summary = serde_json::json!({
            "n_cells": n_cells,
            "n_genes": n_genes,
            "n_hvgs": hvg_names.len(),
            "n_pcs": n_pcs,
            "n_neighbors": k,
            "resolution": args.resolution,
            "n_clusters": n_clusters,
            "variance_explained": variances,
        });
        let mut f = File::create(args.output.join("analysis_summary.json"))?;
        serde_json::to_writer_pretty(&mut f, &summary)
            .context("Failed to write summary JSON")?;
    }

    // Cluster sizes
    let mut cluster_sizes: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for &l in &labels {
        *cluster_sizes.entry(l).or_insert(0) += 1;
    }
    let mut sizes: Vec<_> = cluster_sizes.iter().collect();
    sizes.sort_by_key(|&(&k, _)| k);

    println!("\nCluster sizes:");
    for (cluster, size) in &sizes {
        println!("  Cluster {}: {} cells", cluster, size);
    }

    println!("\n=== Analysis Complete ===");
    println!("Output directory: {:?}", args.output);

    Ok(())
}
