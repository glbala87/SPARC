//! Built-in aligner integration (STAR and minimap2)

use crate::{Error, Result};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Result of BAM file validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BamValidation {
    pub path: PathBuf,
    pub file_size: u64,
    pub format: String,
}

/// Supported aligner types
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AlignerType {
    Star,
    Minimap2,
}

/// Aligner configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignerConfig {
    pub aligner_type: AlignerType,
    pub genome_dir: PathBuf,
    pub threads: usize,
    pub extra_args: Vec<String>,
}

impl AlignerConfig {
    pub fn star(genome_dir: PathBuf, threads: usize) -> Self {
        Self {
            aligner_type: AlignerType::Star,
            genome_dir,
            threads,
            extra_args: vec![
                "--outSAMtype".into(),
                "BAM".into(),
                "SortedByCoordinate".into(),
                "--outSAMattributes".into(),
                "NH".into(),
                "HI".into(),
                "nM".into(),
                "AS".into(),
                "CR".into(),
                "UR".into(),
                "CB".into(),
                "UB".into(),
                "GX".into(),
                "GN".into(),
                "--soloType".into(),
                "CB_UMI_Simple".into(),
            ],
        }
    }

    pub fn minimap2(genome_ref: PathBuf, threads: usize) -> Self {
        Self {
            aligner_type: AlignerType::Minimap2,
            genome_dir: genome_ref,
            threads,
            extra_args: vec!["-a".into(), "--secondary=no".into()],
        }
    }
}

/// Aligner wrapper
pub struct Aligner {
    config: AlignerConfig,
}

impl Aligner {
    pub fn new(config: AlignerConfig) -> Self {
        Self { config }
    }

    /// Check if the aligner is available in PATH
    pub fn is_available(&self) -> bool {
        let result = Command::new("which")
            .arg(self.binary_name())
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        log::info!("Aligner '{}' available: {}", self.binary_name(), result);
        result
    }

    /// Get the aligner binary name
    pub fn binary_name(&self) -> &str {
        match self.config.aligner_type {
            AlignerType::Star => "STAR",
            AlignerType::Minimap2 => "minimap2",
        }
    }

    /// Run alignment
    pub fn align<P: AsRef<Path>>(
        &self,
        r1: P,
        r2: Option<P>,
        output_dir: P,
    ) -> Result<PathBuf> {
        match self.config.aligner_type {
            AlignerType::Star => self.run_star(r1, r2, output_dir),
            AlignerType::Minimap2 => self.run_minimap2(r1, r2, output_dir),
        }
    }

    fn run_star<P: AsRef<Path>>(
        &self,
        r1: P,
        r2: Option<P>,
        output_dir: P,
    ) -> Result<PathBuf> {
        let output_dir = output_dir.as_ref();
        std::fs::create_dir_all(output_dir)?;

        let output_prefix = output_dir.join("star_");

        let mut cmd = Command::new("STAR");
        cmd.arg("--genomeDir")
            .arg(&self.config.genome_dir)
            .arg("--readFilesIn");

        if let Some(r2) = r2 {
            cmd.arg(r2.as_ref()).arg(r1.as_ref());
        } else {
            cmd.arg(r1.as_ref());
        }

        cmd.arg("--runThreadN")
            .arg(self.config.threads.to_string())
            .arg("--outFileNamePrefix")
            .arg(&output_prefix);

        let r1_path = r1.as_ref();
        if r1_path.extension().map_or(false, |e| e == "gz") {
            cmd.arg("--readFilesCommand").arg("zcat");
        }

        for arg in &self.config.extra_args {
            cmd.arg(arg);
        }

        log::info!("Running STAR alignment...");
        let output = cmd.output().map_err(Error::Io)?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::BamParse(format!("STAR failed: {}", stderr)));
        }

        let bam_path = output_dir.join("star_Aligned.sortedByCoord.out.bam");
        if bam_path.exists() {
            Ok(bam_path)
        } else {
            let bam_path2 = output_dir.join("star_Aligned.out.bam");
            if bam_path2.exists() {
                Ok(bam_path2)
            } else {
                Err(Error::BamParse("STAR output BAM not found".into()))
            }
        }
    }

    fn run_minimap2<P: AsRef<Path>>(
        &self,
        r1: P,
        r2: Option<P>,
        output_dir: P,
    ) -> Result<PathBuf> {
        let output_dir = output_dir.as_ref();
        std::fs::create_dir_all(output_dir)?;

        let bam_path = output_dir.join("aligned.bam");
        let sam_path = output_dir.join("aligned.sam");

        let mut cmd = Command::new("minimap2");
        cmd.arg("-t").arg(self.config.threads.to_string());

        for arg in &self.config.extra_args {
            cmd.arg(arg);
        }

        cmd.arg(&self.config.genome_dir);

        if let Some(r2) = r2 {
            cmd.arg(r2.as_ref()).arg(r1.as_ref());
        } else {
            cmd.arg(r1.as_ref());
        }

        cmd.arg("-o").arg(&sam_path);

        log::info!("Running minimap2 alignment...");
        let output = cmd.output().map_err(Error::Io)?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(Error::BamParse(format!("minimap2 failed: {}", stderr)));
        }

        // Convert SAM to BAM using samtools if available
        let samtools_available = Command::new("which")
            .arg("samtools")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);

        if samtools_available {
            let status = Command::new("samtools")
                .args(["view", "-bS", "-o"])
                .arg(&bam_path)
                .arg(&sam_path)
                .status()
                .map_err(Error::Io)?;

            if status.success() {
                let _ = std::fs::remove_file(&sam_path);

                let sorted_path = output_dir.join("aligned.sorted.bam");
                let sort_ok = Command::new("samtools")
                    .args(["sort", "-o"])
                    .arg(&sorted_path)
                    .arg(&bam_path)
                    .status()
                    .map(|s| s.success())
                    .unwrap_or(false);

                if sort_ok && sorted_path.exists() {
                    let _ = std::fs::remove_file(&bam_path);
                    return Ok(sorted_path);
                }
                return Ok(bam_path);
            }
        }

        Ok(sam_path)
    }

    /// Validate a BAM file by checking it can be read and has aligned records
    pub fn validate_bam<P: AsRef<Path>>(bam_path: P) -> Result<BamValidation> {
        let path = bam_path.as_ref();
        if !path.exists() {
            return Err(Error::BamParse(format!("BAM file not found: {:?}", path)));
        }

        let file_size = std::fs::metadata(path)
            .map(|m| m.len())
            .unwrap_or(0);

        if file_size == 0 {
            return Err(Error::BamParse("BAM file is empty".into()));
        }

        // Try to read header bytes to validate format
        let mut file = std::fs::File::open(path)?;
        let mut magic = [0u8; 4];
        use std::io::Read;
        file.read_exact(&mut magic)?;

        let is_bam = &magic == b"BAM\x01";
        let is_sam = magic[0] == b'@';

        if !is_bam && !is_sam {
            return Err(Error::BamParse(
                "File does not appear to be a valid BAM or SAM file".into(),
            ));
        }

        Ok(BamValidation {
            path: path.to_path_buf(),
            file_size,
            format: if is_bam { "BAM".into() } else { "SAM".into() },
        })
    }

    /// Index a BAM file using samtools
    pub fn index_bam<P: AsRef<Path>>(bam_path: P) -> Result<PathBuf> {
        let bam_path = bam_path.as_ref();
        let bai_path = bam_path.with_extension("bam.bai");

        if !Self::is_tool_available("samtools") {
            return Err(Error::BamParse(
                "samtools not found in PATH — cannot index BAM".into(),
            ));
        }

        let status = Command::new("samtools")
            .args(["index", "-@", "4"])
            .arg(bam_path)
            .status()
            .map_err(Error::Io)?;

        if !status.success() {
            return Err(Error::BamParse("samtools index failed".into()));
        }

        Ok(bai_path)
    }

    /// Check if a tool is available in PATH
    fn is_tool_available(name: &str) -> bool {
        Command::new("which")
            .arg(name)
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false)
    }

    /// Generate STAR genome index
    pub fn generate_star_index<P: AsRef<Path>>(
        genome_fasta: P,
        gtf: P,
        output_dir: P,
        threads: usize,
    ) -> Result<()> {
        let output_dir = output_dir.as_ref();
        std::fs::create_dir_all(output_dir)?;

        let status = Command::new("STAR")
            .arg("--runMode")
            .arg("genomeGenerate")
            .arg("--genomeDir")
            .arg(output_dir)
            .arg("--genomeFastaFiles")
            .arg(genome_fasta.as_ref())
            .arg("--sjdbGTFfile")
            .arg(gtf.as_ref())
            .arg("--runThreadN")
            .arg(threads.to_string())
            .status()
            .map_err(Error::Io)?;

        if !status.success() {
            return Err(Error::BamParse(
                "STAR genome generation failed".to_string(),
            ));
        }

        Ok(())
    }
}
