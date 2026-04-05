//! SPARC CLI - Single-cell Pipeline Accelerated in Rust Core

mod commands;

use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "sparc")]
#[command(author, version, about = "SPARC: Single-cell Pipeline Accelerated in Rust Core", long_about = None)]
struct Cli {
    /// Enable verbose output
    #[arg(short, long, global = true)]
    verbose: bool,

    /// Number of threads to use
    #[arg(short = 'j', long, global = true, default_value = "0")]
    threads: usize,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract barcodes and UMIs from FASTQ files
    Extract(commands::extract::ExtractArgs),

    /// Generate gene count matrix
    Count(commands::count::CountArgs),

    /// Generate QC report
    Qc(commands::qc::QcArgs),

    /// Run full analysis pipeline
    Pipeline(commands::pipeline::PipelineArgs),

    /// Process multiple samples from a manifest
    Batch(commands::batch::BatchArgs),

    /// Distributed processing (shard/worker/merge)
    Distributed(commands::distributed::DistributedArgs),

    /// Run downstream analysis (normalize, PCA, clustering) on a count matrix
    Analyze(commands::analyze::AnalyzeArgs),

    /// Run truthset validation against synthetic ground-truth data
    Validate(commands::validate::ValidateArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Initialize logger
    if cli.verbose {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    }

    // Set thread count
    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok();
    }

    match cli.command {
        Commands::Extract(args) => commands::extract::run(args),
        Commands::Count(args) => commands::count::run(args),
        Commands::Qc(args) => commands::qc::run(args),
        Commands::Pipeline(args) => commands::pipeline::run(args),
        Commands::Batch(args) => commands::batch::run(args),
        Commands::Distributed(args) => commands::distributed::run(args),
        Commands::Analyze(args) => commands::analyze::run(args),
        Commands::Validate(args) => commands::validate::run(args),
    }
}
