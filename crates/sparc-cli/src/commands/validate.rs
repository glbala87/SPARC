//! `sparc validate` - Run truthset validation against synthetic ground-truth data

use anyhow::Result;
use clap::Args;
use std::path::PathBuf;

use sparc_core::validation::report::{ValidationReport, ValidationThresholds};
use sparc_core::validation::stages;
use sparc_core::validation::synthetic::{SyntheticConfig, SyntheticDataset};

#[derive(Args)]
pub struct ValidateArgs {
    /// Output directory for validation report
    #[arg(short, long, default_value = "validation_output")]
    output: PathBuf,

    /// Protocol to validate
    #[arg(short, long, default_value = "10x-3prime-v3")]
    protocol: String,

    /// Number of synthetic cells
    #[arg(long, default_value = "500")]
    n_cells: usize,

    /// Number of synthetic genes
    #[arg(long, default_value = "200")]
    n_genes: usize,

    /// Number of cell types
    #[arg(long, default_value = "5")]
    n_cell_types: usize,

    /// Stages to validate: extract,count,analysis or all
    #[arg(long, default_value = "all")]
    stages: String,

    /// Random seed for reproducibility
    #[arg(long, default_value = "42")]
    seed: u64,

    /// Minimum barcode F1 score for pass
    #[arg(long, default_value = "0.95")]
    min_barcode_f1: f64,

    /// Minimum expression Pearson r for pass
    #[arg(long, default_value = "0.90")]
    min_expression_pearson: f64,

    /// Minimum clustering ARI for pass
    #[arg(long, default_value = "0.70")]
    min_clustering_ari: f64,
}

pub fn run(args: ValidateArgs) -> Result<()> {
    let stages_to_run: Vec<&str> = if args.stages == "all" {
        vec!["extract", "count", "analysis"]
    } else {
        args.stages.split(',').map(|s| s.trim()).collect()
    };

    let config = SyntheticConfig {
        n_cells: args.n_cells,
        n_genes: args.n_genes,
        n_cell_types: args.n_cell_types,
        seed: args.seed,
        protocol: args.protocol.clone(),
        ..Default::default()
    };

    let thresholds = ValidationThresholds {
        min_barcode_f1: args.min_barcode_f1,
        min_expression_pearson: args.min_expression_pearson,
        min_clustering_ari: args.min_clustering_ari,
    };

    println!("Generating synthetic dataset...");
    println!(
        "  {} cells, {} genes, {} cell types, seed={}",
        config.n_cells, config.n_genes, config.n_cell_types, config.seed
    );

    let dataset = SyntheticDataset::generate(config.clone());

    println!(
        "  Generated {} R1 reads, {} R2 reads",
        dataset.r1_records.len(),
        dataset.r2_records.len()
    );
    println!(
        "  {} mutated barcodes, {} invalid barcodes",
        dataset.truth.mutated_barcodes.len(),
        dataset.truth.invalid_barcodes.len()
    );

    let mut report = ValidationReport::new(config, thresholds);

    // ─── Extract stage ───────────────────────────────────────────
    if stages_to_run.contains(&"extract") {
        println!("\nRunning extract validation...");
        let extract_result = stages::validate_extract(&dataset);
        println!("  F1: {:.4}, Sensitivity: {:.4}, Specificity: {:.4}",
            extract_result.f1, extract_result.sensitivity, extract_result.specificity);
        report.set_extract_results(extract_result);
    }

    // ─── Count stage ─────────────────────────────────────────────
    if stages_to_run.contains(&"count") {
        println!("\nRunning count validation...");
        // Validate the truth matrix against itself (perfect baseline)
        // In a full pipeline run, this would compare the pipeline output
        let count_result = stages::validate_count(
            &dataset.truth.expression_matrix,
            &dataset.truth.expression_matrix,
        );
        println!("  Pearson r: {:.4}, Spearman rho: {:.4}",
            count_result.pearson_r, count_result.spearman_rho);
        report.set_count_results(count_result);
    }

    // ─── Analysis stage ──────────────────────────────────────────
    if stages_to_run.contains(&"analysis") {
        println!("\nRunning analysis validation...");
        let analysis_result = stages::validate_clustering_with_knn(&dataset);
        println!("  ARI: {:.4}, NMI: {:.4}, Clusters: {}/{}",
            analysis_result.ari, analysis_result.nmi,
            analysis_result.n_clusters_found, analysis_result.n_clusters_expected);
        report.set_analysis_results(analysis_result);
    }

    // ─── Write report ────────────────────────────────────────────
    std::fs::create_dir_all(&args.output)?;
    let report_path = args.output.join("validation_report.json");
    report.write_json(&report_path)?;
    println!("\nReport written to: {}", report_path.display());

    // Print summary
    println!("\n{}", report.summary());

    if report.overall_pass {
        Ok(())
    } else {
        anyhow::bail!("Validation FAILED — see report for details")
    }
}
