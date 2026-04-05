//! Synthetic ground-truth dataset generation
//!
//! Generates deterministic synthetic scRNA-seq data with known cell types,
//! barcode sequences, UMI counts, and expression profiles for validation.

use std::collections::HashMap;

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::Poisson;
use serde::{Deserialize, Serialize};

use crate::barcode::Whitelist;
use crate::count::CountMatrix;
use crate::fastq::FastqRecord;
use crate::ReadStructure;

/// Configuration for synthetic data generation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyntheticConfig {
    /// Number of cells to simulate
    pub n_cells: usize,
    /// Number of genes to simulate
    pub n_genes: usize,
    /// Number of distinct cell types
    pub n_cell_types: usize,
    /// Number of marker genes per cell type (elevated expression)
    pub n_markers_per_type: usize,
    /// Fold-change for marker genes
    pub marker_fold_change: f64,
    /// Base mean expression for non-marker genes
    pub base_expression: f64,
    /// Fraction of barcodes to mutate (for correction testing)
    pub mutation_rate: f64,
    /// Fraction of completely invalid barcodes (for specificity testing)
    pub invalid_barcode_rate: f64,
    /// Barcode length
    pub barcode_len: usize,
    /// UMI length
    pub umi_len: usize,
    /// Random seed for reproducibility
    pub seed: u64,
    /// Protocol name
    pub protocol: String,
}

impl Default for SyntheticConfig {
    fn default() -> Self {
        Self {
            n_cells: 500,
            n_genes: 200,
            n_cell_types: 5,
            n_markers_per_type: 10,
            marker_fold_change: 8.0,
            base_expression: 2.0,
            mutation_rate: 0.1,
            invalid_barcode_rate: 0.05,
            barcode_len: 16,
            umi_len: 12,
            seed: 42,
            protocol: "10x-3prime-v3".to_string(),
        }
    }
}

/// Ground-truth labels and expected outputs
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TruthSet {
    /// Valid barcode sequences
    pub barcodes: Vec<String>,
    /// Gene names
    pub genes: Vec<String>,
    /// Barcode -> cell type index
    pub cell_types: HashMap<String, usize>,
    /// Cell type names
    pub cell_type_names: Vec<String>,
    /// Ground truth count matrix (genes x cells)
    pub expression_matrix: CountMatrix,
    /// UMI counts per (barcode, gene) pair
    pub umi_counts: HashMap<(String, String), u32>,
    /// Mutated barcodes: (original, mutated, hamming_distance)
    pub mutated_barcodes: Vec<(String, String, u32)>,
    /// Invalid barcodes (not in whitelist, not correctable)
    pub invalid_barcodes: Vec<String>,
    /// Mean expression profile per cell type (cell_type_idx -> gene_idx -> mean)
    pub expression_profiles: Vec<Vec<f64>>,
}

/// A complete synthetic dataset for validation
pub struct SyntheticDataset {
    /// Barcode whitelist
    pub whitelist: Whitelist,
    /// R1 FASTQ records (barcode + UMI)
    pub r1_records: Vec<FastqRecord>,
    /// R2 FASTQ records (cDNA / gene tag)
    pub r2_records: Vec<FastqRecord>,
    /// Ground truth
    pub truth: TruthSet,
    /// Read structure
    pub read_structure: ReadStructure,
    /// Configuration used
    pub config: SyntheticConfig,
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

impl SyntheticDataset {
    /// Generate a complete synthetic dataset
    pub fn generate(config: SyntheticConfig) -> Self {
        let mut rng = StdRng::seed_from_u64(config.seed);

        // 1. Generate whitelist barcodes
        let barcodes = generate_barcodes(&mut rng, config.n_cells, config.barcode_len);
        let whitelist = Whitelist::from_vec(barcodes.clone()).expect("Valid barcodes");

        // 2. Generate gene names
        let genes: Vec<String> = (0..config.n_genes)
            .map(|i| format!("GENE_{:04}", i))
            .collect();

        // 3. Assign cell types
        let cell_type_names: Vec<String> = (0..config.n_cell_types)
            .map(|i| format!("CellType_{}", i))
            .collect();

        let mut cell_types = HashMap::new();
        for (i, bc) in barcodes.iter().enumerate() {
            let ct = i % config.n_cell_types;
            cell_types.insert(bc.clone(), ct);
        }

        // 4. Build expression profiles per cell type
        let expression_profiles = build_expression_profiles(
            &config,
            &mut rng,
        );

        // 5. Sample counts from Poisson distribution
        let (count_matrix, umi_counts) = sample_expression(
            &barcodes,
            &genes,
            &cell_types,
            &expression_profiles,
            &mut rng,
        );

        // 6. Generate mutated barcodes
        let mutated_barcodes = generate_mutations(&barcodes, &config, &mut rng);

        // 7. Generate invalid barcodes
        let n_invalid = (config.n_cells as f64 * config.invalid_barcode_rate) as usize;
        let invalid_barcodes = generate_invalid_barcodes(
            &mut rng,
            n_invalid,
            config.barcode_len,
            &whitelist,
        );

        // 8. Build FASTQ records
        let (r1_records, r2_records) = build_fastq_records(
            &barcodes,
            &genes,
            &umi_counts,
            &mutated_barcodes,
            &invalid_barcodes,
            &config,
            &mut rng,
        );

        let read_structure = ReadStructure::new(
            0,
            config.barcode_len,
            config.barcode_len,
            config.umi_len,
            0,
        );

        let truth = TruthSet {
            barcodes: barcodes.clone(),
            genes,
            cell_types,
            cell_type_names,
            expression_matrix: count_matrix,
            umi_counts,
            mutated_barcodes,
            invalid_barcodes,
            expression_profiles,
        };

        Self {
            whitelist,
            r1_records,
            r2_records,
            truth,
            read_structure,
            config,
        }
    }
}

/// Generate random DNA barcodes
fn generate_barcodes(rng: &mut StdRng, n: usize, len: usize) -> Vec<String> {
    let mut barcodes = Vec::with_capacity(n);
    let mut seen = ahash::AHashSet::new();

    while barcodes.len() < n {
        let bc: String = (0..len)
            .map(|_| BASES[rng.gen_range(0..4)] as char)
            .collect();
        if seen.insert(bc.clone()) {
            barcodes.push(bc);
        }
    }
    barcodes
}

/// Build mean expression profiles for each cell type
fn build_expression_profiles(
    config: &SyntheticConfig,
    rng: &mut StdRng,
) -> Vec<Vec<f64>> {
    let mut profiles = Vec::with_capacity(config.n_cell_types);

    for ct in 0..config.n_cell_types {
        let mut profile = vec![config.base_expression; config.n_genes];

        // Assign marker genes for this cell type
        let marker_start = ct * config.n_markers_per_type;
        let marker_end = (marker_start + config.n_markers_per_type).min(config.n_genes);

        for gene_idx in marker_start..marker_end {
            profile[gene_idx] = config.base_expression * config.marker_fold_change;
        }

        // Add some random variation to non-marker genes
        for val in profile.iter_mut() {
            let noise: f64 = 1.0 + rng.gen_range(-0.3..0.3);
            *val = (*val * noise).max(0.1);
        }

        profiles.push(profile);
    }
    profiles
}

/// Sample expression counts from Poisson distribution
fn sample_expression(
    barcodes: &[String],
    genes: &[String],
    cell_types: &HashMap<String, usize>,
    profiles: &[Vec<f64>],
    rng: &mut StdRng,
) -> (CountMatrix, HashMap<(String, String), u32>) {
    let n_genes = genes.len();
    let n_cells = barcodes.len();

    // Dense matrix: genes x cells
    let mut dense = vec![vec![0u32; n_cells]; n_genes];
    let mut umi_counts = HashMap::new();

    for (cell_idx, bc) in barcodes.iter().enumerate() {
        let ct = cell_types[bc];
        let profile = &profiles[ct];

        for (gene_idx, &mean) in profile.iter().enumerate() {
            if mean < 0.01 {
                continue;
            }
            let poisson = Poisson::new(mean).unwrap_or_else(|_| Poisson::new(0.1).unwrap());
            let count: u32 = rng.sample::<f64, _>(poisson) as u32;
            if count > 0 {
                dense[gene_idx][cell_idx] = count;
                umi_counts.insert(
                    (bc.clone(), genes[gene_idx].clone()),
                    count,
                );
            }
        }
    }

    let matrix = CountMatrix::from_dense(barcodes.to_vec(), genes.to_vec(), dense);
    (matrix, umi_counts)
}

/// Generate 1-bp mutated versions of selected barcodes
fn generate_mutations(
    barcodes: &[String],
    config: &SyntheticConfig,
    rng: &mut StdRng,
) -> Vec<(String, String, u32)> {
    let n_mutate = (barcodes.len() as f64 * config.mutation_rate) as usize;
    let mut mutations = Vec::with_capacity(n_mutate);

    for i in 0..n_mutate.min(barcodes.len()) {
        let original = &barcodes[i];
        let mut chars: Vec<u8> = original.bytes().collect();

        // Mutate one random position
        let pos = rng.gen_range(0..chars.len());
        let original_base = chars[pos];
        loop {
            let new_base = BASES[rng.gen_range(0..4)];
            if new_base != original_base {
                chars[pos] = new_base;
                break;
            }
        }

        let mutated = String::from_utf8(chars).unwrap();
        mutations.push((original.clone(), mutated, 1));
    }
    mutations
}

/// Generate barcodes that don't match anything in the whitelist
fn generate_invalid_barcodes(
    rng: &mut StdRng,
    n: usize,
    len: usize,
    whitelist: &Whitelist,
) -> Vec<String> {
    let mut invalid = Vec::with_capacity(n);

    while invalid.len() < n {
        // Use 'N' characters to ensure no correction is possible
        let bc: String = (0..len)
            .map(|_| {
                if rng.gen_bool(0.3) {
                    'N'
                } else {
                    BASES[rng.gen_range(0..4)] as char
                }
            })
            .collect();
        if !whitelist.contains(&bc) {
            invalid.push(bc);
        }
    }
    invalid
}

/// Generate random UMI sequence
fn generate_umi(rng: &mut StdRng, len: usize) -> Vec<u8> {
    (0..len).map(|_| BASES[rng.gen_range(0..4)]).collect()
}

/// Build FASTQ records from the truth data
fn build_fastq_records(
    barcodes: &[String],
    genes: &[String],
    umi_counts: &HashMap<(String, String), u32>,
    mutated_barcodes: &[(String, String, u32)],
    invalid_barcodes: &[String],
    config: &SyntheticConfig,
    rng: &mut StdRng,
) -> (Vec<FastqRecord>, Vec<FastqRecord>) {
    let mut r1_records = Vec::new();
    let mut r2_records = Vec::new();
    let mut read_idx: u64 = 0;

    let high_qual = vec![b'I'; config.barcode_len + config.umi_len]; // Q40
    let r2_len = 100;
    let r2_qual = vec![b'I'; r2_len];

    // Reads from valid barcodes (exact matches)
    for bc in barcodes.iter() {
        for gene in genes.iter() {
            let key = (bc.clone(), gene.clone());
            if let Some(&count) = umi_counts.get(&key) {
                for _ in 0..count {
                    let umi = generate_umi(rng, config.umi_len);
                    let mut r1_seq = bc.as_bytes().to_vec();
                    r1_seq.extend_from_slice(&umi);

                    let read_id = format!("READ_{:08}:{}:{}", read_idx, bc, gene);

                    r1_records.push(FastqRecord::new(
                        read_id.clone(),
                        r1_seq,
                        high_qual.clone(),
                    ));

                    // R2: random cDNA tagged with gene in read name
                    let r2_seq: Vec<u8> = (0..r2_len).map(|_| BASES[rng.gen_range(0..4)]).collect();
                    r2_records.push(FastqRecord::new(read_id, r2_seq, r2_qual.clone()));

                    read_idx += 1;
                }
            }
        }
    }

    // Reads from mutated barcodes (should be corrected)
    for (original, mutated, _dist) in mutated_barcodes {
        // Pick a gene that this cell expresses
        for gene in genes.iter().take(3) {
            let key = (original.clone(), gene.clone());
            if umi_counts.contains_key(&key) {
                let umi = generate_umi(rng, config.umi_len);
                let mut r1_seq = mutated.as_bytes().to_vec();
                r1_seq.extend_from_slice(&umi);

                let read_id = format!("READ_{:08}:{}:{}", read_idx, mutated, gene);

                r1_records.push(FastqRecord::new(
                    read_id.clone(),
                    r1_seq,
                    high_qual.clone(),
                ));

                let r2_seq: Vec<u8> = (0..r2_len).map(|_| BASES[rng.gen_range(0..4)]).collect();
                r2_records.push(FastqRecord::new(read_id, r2_seq, r2_qual.clone()));

                read_idx += 1;
            }
        }
    }

    // Reads from invalid barcodes (should be rejected)
    for invalid_bc in invalid_barcodes {
        let umi = generate_umi(rng, config.umi_len);
        let mut r1_seq = invalid_bc.as_bytes().to_vec();
        r1_seq.extend_from_slice(&umi);

        let read_id = format!("READ_{:08}:{}:NONE", read_idx, invalid_bc);

        r1_records.push(FastqRecord::new(
            read_id.clone(),
            r1_seq,
            high_qual.clone(),
        ));

        let r2_seq: Vec<u8> = (0..r2_len).map(|_| BASES[rng.gen_range(0..4)]).collect();
        r2_records.push(FastqRecord::new(read_id, r2_seq, r2_qual.clone()));

        read_idx += 1;
    }

    (r1_records, r2_records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_default() {
        let config = SyntheticConfig {
            n_cells: 50,
            n_genes: 20,
            n_cell_types: 3,
            n_markers_per_type: 5,
            seed: 42,
            ..Default::default()
        };
        let dataset = SyntheticDataset::generate(config);

        assert_eq!(dataset.truth.barcodes.len(), 50);
        assert_eq!(dataset.truth.genes.len(), 20);
        assert_eq!(dataset.truth.cell_type_names.len(), 3);
        assert_eq!(dataset.truth.expression_matrix.n_rows, 20);
        assert_eq!(dataset.truth.expression_matrix.n_cols, 50);
        assert!(!dataset.r1_records.is_empty());
        assert_eq!(dataset.r1_records.len(), dataset.r2_records.len());
    }

    #[test]
    fn test_deterministic() {
        let config = SyntheticConfig {
            n_cells: 10,
            n_genes: 5,
            seed: 123,
            ..Default::default()
        };
        let d1 = SyntheticDataset::generate(config.clone());
        let d2 = SyntheticDataset::generate(config);

        assert_eq!(d1.truth.barcodes, d2.truth.barcodes);
        assert_eq!(d1.truth.expression_matrix.values, d2.truth.expression_matrix.values);
    }

    #[test]
    fn test_marker_genes_enriched() {
        let config = SyntheticConfig {
            n_cells: 100,
            n_genes: 50,
            n_cell_types: 3,
            n_markers_per_type: 5,
            marker_fold_change: 10.0,
            base_expression: 2.0,
            seed: 42,
            ..Default::default()
        };
        let dataset = SyntheticDataset::generate(config);

        // Check that marker genes have higher mean expression for their cell type
        for ct in 0..3 {
            let profile = &dataset.truth.expression_profiles[ct];
            let marker_start = ct * 5;
            let marker_mean: f64 = profile[marker_start..marker_start + 5].iter().sum::<f64>() / 5.0;
            let nonmarker_mean: f64 = profile.iter()
                .enumerate()
                .filter(|(i, _)| *i < marker_start || *i >= marker_start + 5)
                .map(|(_, v)| v)
                .sum::<f64>() / (50 - 5) as f64;
            assert!(
                marker_mean > nonmarker_mean * 3.0,
                "Marker genes should be enriched: marker_mean={}, nonmarker_mean={}",
                marker_mean,
                nonmarker_mean
            );
        }
    }
}
