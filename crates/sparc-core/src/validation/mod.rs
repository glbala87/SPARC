//! Truthset validation framework for SPARC
//!
//! Generates synthetic ground-truth datasets and validates each pipeline stage
//! (extract, count, analysis) against expected results with accuracy metrics.

pub mod metrics;
pub mod report;
pub mod stages;
pub mod synthetic;

pub use metrics::*;
pub use report::ValidationReport;
pub use synthetic::{SyntheticConfig, SyntheticDataset, TruthSet};
