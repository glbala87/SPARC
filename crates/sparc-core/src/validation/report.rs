//! Validation report generation and serialization

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use super::stages::{AnalysisValidationResult, CountValidationResult, ExtractValidationResult};
use super::synthetic::SyntheticConfig;

/// Thresholds for pass/fail determination
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationThresholds {
    /// Minimum barcode detection F1 score
    pub min_barcode_f1: f64,
    /// Minimum expression Pearson correlation
    pub min_expression_pearson: f64,
    /// Minimum clustering ARI
    pub min_clustering_ari: f64,
}

impl Default for ValidationThresholds {
    fn default() -> Self {
        Self {
            min_barcode_f1: 0.95,
            min_expression_pearson: 0.90,
            min_clustering_ari: 0.70,
        }
    }
}

/// Complete validation report
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationReport {
    /// Timestamp of validation run
    pub timestamp: String,
    /// SPARC version
    pub sparc_version: String,
    /// Protocol validated
    pub protocol: String,
    /// Configuration used for synthetic data
    pub synthetic_config: SyntheticConfig,
    /// Extract stage results
    pub extract_results: Option<ExtractValidationResult>,
    /// Count stage results
    pub count_results: Option<CountValidationResult>,
    /// Analysis stage results
    pub analysis_results: Option<AnalysisValidationResult>,
    /// Overall pass/fail
    pub overall_pass: bool,
    /// Thresholds used
    pub thresholds: ValidationThresholds,
    /// Per-stage pass/fail
    pub stage_results: StagePassFail,
}

/// Per-stage pass/fail summary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StagePassFail {
    pub extract: Option<bool>,
    pub count: Option<bool>,
    pub analysis: Option<bool>,
}

impl ValidationReport {
    /// Create a new report
    pub fn new(config: SyntheticConfig, thresholds: ValidationThresholds) -> Self {
        let timestamp = chrono::Utc::now().to_rfc3339();

        Self {
            timestamp,
            sparc_version: env!("CARGO_PKG_VERSION").to_string(),
            protocol: config.protocol.clone(),
            synthetic_config: config,
            extract_results: None,
            count_results: None,
            analysis_results: None,
            overall_pass: false,
            thresholds,
            stage_results: StagePassFail {
                extract: None,
                count: None,
                analysis: None,
            },
        }
    }

    /// Set extract results and evaluate pass/fail
    pub fn set_extract_results(&mut self, results: ExtractValidationResult) {
        let pass = results.f1 >= self.thresholds.min_barcode_f1;
        self.stage_results.extract = Some(pass);
        self.extract_results = Some(results);
        self.evaluate_overall();
    }

    /// Set count results and evaluate pass/fail
    pub fn set_count_results(&mut self, results: CountValidationResult) {
        let pass = results.pearson_r >= self.thresholds.min_expression_pearson;
        self.stage_results.count = Some(pass);
        self.count_results = Some(results);
        self.evaluate_overall();
    }

    /// Set analysis results and evaluate pass/fail
    pub fn set_analysis_results(&mut self, results: AnalysisValidationResult) {
        let pass = results.ari >= self.thresholds.min_clustering_ari;
        self.stage_results.analysis = Some(pass);
        self.analysis_results = Some(results);
        self.evaluate_overall();
    }

    /// Evaluate overall pass: all present stages must pass
    fn evaluate_overall(&mut self) {
        let results = [
            self.stage_results.extract,
            self.stage_results.count,
            self.stage_results.analysis,
        ];
        // Overall passes if at least one stage ran and all ran stages passed
        let ran_any = results.iter().any(|r| r.is_some());
        let all_pass = results.iter().all(|r| r.unwrap_or(true));
        self.overall_pass = ran_any && all_pass;
    }

    /// Serialize to JSON string
    pub fn to_json(&self) -> serde_json::Result<String> {
        serde_json::to_string_pretty(self)
    }

    /// Write JSON report to file
    pub fn write_json<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let file = File::create(path)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
    }

    /// Format a human-readable summary
    pub fn summary(&self) -> String {
        let mut lines = Vec::new();
        lines.push("═══════════════════════════════════════════════════".to_string());
        lines.push("           SPARC Truthset Validation Report        ".to_string());
        lines.push("═══════════════════════════════════════════════════".to_string());
        lines.push(format!("Timestamp:  {}", self.timestamp));
        lines.push(format!("Version:    {}", self.sparc_version));
        lines.push(format!("Protocol:   {}", self.protocol));
        lines.push(format!(
            "Synthetic:  {} cells, {} genes, {} types",
            self.synthetic_config.n_cells,
            self.synthetic_config.n_genes,
            self.synthetic_config.n_cell_types
        ));
        lines.push("───────────────────────────────────────────────────".to_string());

        if let Some(ref ext) = self.extract_results {
            let pass = self.stage_results.extract.unwrap_or(false);
            let icon = if pass { "PASS" } else { "FAIL" };
            lines.push(format!("\n[{}] Extract / Barcode Detection", icon));
            lines.push(format!("  Sensitivity:       {:.4}", ext.sensitivity));
            lines.push(format!("  Specificity:       {:.4}", ext.specificity));
            lines.push(format!("  Precision:         {:.4}", ext.precision));
            lines.push(format!("  Recall:            {:.4}", ext.recall));
            lines.push(format!(
                "  F1 Score:          {:.4}  (threshold: {:.2})",
                ext.f1, self.thresholds.min_barcode_f1
            ));
            lines.push(format!("  Correction Acc:    {:.4}", ext.correction_accuracy));
            lines.push(format!("  Total Reads:       {}", ext.total_reads));
        }

        if let Some(ref cnt) = self.count_results {
            let pass = self.stage_results.count.unwrap_or(false);
            let icon = if pass { "PASS" } else { "FAIL" };
            lines.push(format!("\n[{}] Count Matrix Quantification", icon));
            lines.push(format!(
                "  Pearson r:         {:.4}  (threshold: {:.2})",
                cnt.pearson_r, self.thresholds.min_expression_pearson
            ));
            lines.push(format!("  Spearman rho:      {:.4}", cnt.spearman_rho));
            lines.push(format!("  MAE:               {:.4}", cnt.mae));
            lines.push(format!("  RMSE:              {:.4}", cnt.rmse));
            lines.push(format!("  Cells concordant:  {:.4}", cnt.cells_concordant));
            lines.push(format!("  Genes concordant:  {:.4}", cnt.genes_concordant));
            lines.push(format!(
                "  Cells (truth/obs): {}/{}",
                cnt.truth_cells, cnt.observed_cells
            ));
        }

        if let Some(ref ana) = self.analysis_results {
            let pass = self.stage_results.analysis.unwrap_or(false);
            let icon = if pass { "PASS" } else { "FAIL" };
            lines.push(format!("\n[{}] Analysis / Clustering", icon));
            lines.push(format!(
                "  ARI:               {:.4}  (threshold: {:.2})",
                ana.ari, self.thresholds.min_clustering_ari
            ));
            lines.push(format!("  NMI:               {:.4}", ana.nmi));
            lines.push(format!(
                "  Clusters (exp/found): {}/{}",
                ana.n_clusters_expected, ana.n_clusters_found
            ));
        }

        lines.push("───────────────────────────────────────────────────".to_string());
        let overall_icon = if self.overall_pass { "PASS" } else { "FAIL" };
        lines.push(format!("Overall: [{}]", overall_icon));
        lines.push("═══════════════════════════════════════════════════".to_string());

        lines.join("\n")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_report_json_roundtrip() {
        let config = SyntheticConfig::default();
        let thresholds = ValidationThresholds::default();
        let report = ValidationReport::new(config, thresholds);

        let json = report.to_json().unwrap();
        let parsed: ValidationReport = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.protocol, report.protocol);
    }

    #[test]
    fn test_report_pass_fail() {
        let config = SyntheticConfig::default();
        let thresholds = ValidationThresholds::default();
        let mut report = ValidationReport::new(config, thresholds);

        report.set_extract_results(ExtractValidationResult {
            tp: 90,
            fp: 2,
            tn: 10,
            fn_count: 3,
            sensitivity: 0.97,
            specificity: 0.83,
            precision: 0.978,
            recall: 0.97,
            f1: 0.974,
            correction_accuracy: 0.95,
            total_reads: 105,
        });

        assert_eq!(report.stage_results.extract, Some(true));
        assert!(report.overall_pass);
    }

    #[test]
    fn test_report_summary_format() {
        let config = SyntheticConfig::default();
        let thresholds = ValidationThresholds::default();
        let report = ValidationReport::new(config, thresholds);
        let summary = report.summary();
        assert!(summary.contains("SPARC Truthset Validation Report"));
    }
}
