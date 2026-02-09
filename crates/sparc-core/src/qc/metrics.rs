//! Quality control metrics calculation

use serde::{Deserialize, Serialize};

/// Quality control metrics for a single-cell dataset
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct QcMetrics {
    /// Total number of reads
    pub total_reads: u64,
    /// Reads with valid barcodes
    pub valid_barcode_reads: u64,
    /// Reads with valid UMIs
    pub valid_umi_reads: u64,
    /// Reads mapped to genome
    pub mapped_reads: u64,
    /// Reads assigned to genes
    pub assigned_reads: u64,
    /// Number of detected cells
    pub num_cells: u64,
    /// Mean reads per cell
    pub mean_reads_per_cell: f64,
    /// Median reads per cell
    pub median_reads_per_cell: f64,
    /// Mean genes per cell
    pub mean_genes_per_cell: f64,
    /// Median genes per cell
    pub median_genes_per_cell: f64,
    /// Total genes detected
    pub total_genes: u64,
    /// Mean UMI counts per cell
    pub mean_umi_per_cell: f64,
    /// Median UMI counts per cell
    pub median_umi_per_cell: f64,
    /// Sequencing saturation
    pub sequencing_saturation: f64,
    /// Fraction of reads in cells
    pub fraction_reads_in_cells: f64,
}

impl QcMetrics {
    pub fn new() -> Self {
        Self::default()
    }

    /// Calculate barcode validity rate
    pub fn barcode_validity_rate(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            self.valid_barcode_reads as f64 / self.total_reads as f64
        }
    }

    /// Calculate mapping rate
    pub fn mapping_rate(&self) -> f64 {
        if self.total_reads == 0 {
            0.0
        } else {
            self.mapped_reads as f64 / self.total_reads as f64
        }
    }

    /// Calculate gene assignment rate
    pub fn assignment_rate(&self) -> f64 {
        if self.mapped_reads == 0 {
            0.0
        } else {
            self.assigned_reads as f64 / self.mapped_reads as f64
        }
    }

    /// Update metrics from cell counts
    pub fn update_from_cells(&mut self, reads_per_cell: &[u64], genes_per_cell: &[u64], umis_per_cell: &[u64]) {
        self.num_cells = reads_per_cell.len() as u64;

        if !reads_per_cell.is_empty() {
            let mut sorted = reads_per_cell.to_vec();
            sorted.sort();
            self.mean_reads_per_cell = sorted.iter().sum::<u64>() as f64 / sorted.len() as f64;
            self.median_reads_per_cell = sorted[sorted.len() / 2] as f64;
        }

        if !genes_per_cell.is_empty() {
            let mut sorted = genes_per_cell.to_vec();
            sorted.sort();
            self.mean_genes_per_cell = sorted.iter().sum::<u64>() as f64 / sorted.len() as f64;
            self.median_genes_per_cell = sorted[sorted.len() / 2] as f64;
        }

        if !umis_per_cell.is_empty() {
            let mut sorted = umis_per_cell.to_vec();
            sorted.sort();
            self.mean_umi_per_cell = sorted.iter().sum::<u64>() as f64 / sorted.len() as f64;
            self.median_umi_per_cell = sorted[sorted.len() / 2] as f64;
        }
    }

    /// Calculate sequencing saturation
    /// saturation = 1 - (unique_reads / total_reads)
    pub fn calculate_saturation(&mut self, unique_reads: u64, total_reads: u64) {
        if total_reads == 0 {
            self.sequencing_saturation = 0.0;
        } else {
            self.sequencing_saturation = 1.0 - (unique_reads as f64 / total_reads as f64);
        }
    }
}

/// QC report containing metrics and summary statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcReport {
    /// Sample name
    pub sample_name: String,
    /// QC metrics
    pub metrics: QcMetrics,
    /// Per-cell metrics (cell_barcode -> (reads, genes, umis))
    pub per_cell_metrics: Vec<CellMetrics>,
    /// Warnings
    pub warnings: Vec<String>,
}

/// Metrics for a single cell
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellMetrics {
    /// Cell barcode
    pub barcode: String,
    /// Number of reads
    pub reads: u64,
    /// Number of genes
    pub genes: u64,
    /// Number of UMIs
    pub umis: u64,
    /// Mitochondrial gene percentage
    pub mito_percent: f64,
}

impl QcReport {
    pub fn new(sample_name: String) -> Self {
        Self {
            sample_name,
            metrics: QcMetrics::new(),
            per_cell_metrics: Vec::new(),
            warnings: Vec::new(),
        }
    }

    /// Add a warning
    pub fn add_warning(&mut self, warning: String) {
        self.warnings.push(warning);
    }

    /// Generate warnings based on metrics thresholds
    pub fn generate_warnings(&mut self) {
        if self.metrics.barcode_validity_rate() < 0.5 {
            self.warnings.push("Low barcode validity rate (<50%)".to_string());
        }
        if self.metrics.mapping_rate() < 0.5 {
            self.warnings.push("Low mapping rate (<50%)".to_string());
        }
        if self.metrics.sequencing_saturation < 0.3 {
            self.warnings.push("Low sequencing saturation (<30%)".to_string());
        }
        if self.metrics.median_genes_per_cell < 200.0 {
            self.warnings.push("Low median genes per cell (<200)".to_string());
        }
    }

    /// Export to JSON
    pub fn to_json(&self) -> serde_json::Result<String> {
        serde_json::to_string_pretty(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metrics_calculation() {
        let mut metrics = QcMetrics::new();
        metrics.total_reads = 1000;
        metrics.valid_barcode_reads = 900;
        metrics.mapped_reads = 800;
        metrics.assigned_reads = 700;

        assert!((metrics.barcode_validity_rate() - 0.9).abs() < 0.001);
        assert!((metrics.mapping_rate() - 0.8).abs() < 0.001);
        assert!((metrics.assignment_rate() - 0.875).abs() < 0.001);
    }

    #[test]
    fn test_update_from_cells() {
        let mut metrics = QcMetrics::new();

        let reads_per_cell = vec![100, 200, 300, 400, 500];
        let genes_per_cell = vec![50, 100, 150, 200, 250];
        let umis_per_cell = vec![80, 160, 240, 320, 400];

        metrics.update_from_cells(&reads_per_cell, &genes_per_cell, &umis_per_cell);

        assert_eq!(metrics.num_cells, 5);
        assert!((metrics.mean_reads_per_cell - 300.0).abs() < 0.001);
        assert!((metrics.median_reads_per_cell - 300.0).abs() < 0.001);
    }
}
