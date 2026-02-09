//! Count matrix generation

use ahash::AHashMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::{Error, Result};

/// Sparse count matrix in COO format
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CountMatrix {
    /// Cell barcodes (column names)
    pub barcodes: Vec<String>,
    /// Gene names/IDs (row names)
    pub genes: Vec<String>,
    /// Row indices (gene indices)
    pub rows: Vec<usize>,
    /// Column indices (cell indices)
    pub cols: Vec<usize>,
    /// Count values
    pub values: Vec<u32>,
    /// Number of rows (genes)
    pub n_rows: usize,
    /// Number of columns (cells)
    pub n_cols: usize,
}

impl CountMatrix {
    pub fn new() -> Self {
        Self {
            barcodes: Vec::new(),
            genes: Vec::new(),
            rows: Vec::new(),
            cols: Vec::new(),
            values: Vec::new(),
            n_rows: 0,
            n_cols: 0,
        }
    }

    /// Create from dense matrix (for small matrices)
    pub fn from_dense(barcodes: Vec<String>, genes: Vec<String>, data: Vec<Vec<u32>>) -> Self {
        let n_rows = genes.len();
        let n_cols = barcodes.len();

        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut values = Vec::new();

        for (i, row) in data.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                if val > 0 {
                    rows.push(i);
                    cols.push(j);
                    values.push(val);
                }
            }
        }

        Self {
            barcodes,
            genes,
            rows,
            cols,
            values,
            n_rows,
            n_cols,
        }
    }

    /// Get count for a specific gene and cell
    pub fn get(&self, gene_idx: usize, cell_idx: usize) -> u32 {
        for (i, (&r, &c)) in self.rows.iter().zip(self.cols.iter()).enumerate() {
            if r == gene_idx && c == cell_idx {
                return self.values[i];
            }
        }
        0
    }

    /// Get total counts per cell
    pub fn counts_per_cell(&self) -> Vec<u64> {
        let mut counts = vec![0u64; self.n_cols];
        for (i, &c) in self.cols.iter().enumerate() {
            counts[c] += self.values[i] as u64;
        }
        counts
    }

    /// Get total counts per gene
    pub fn counts_per_gene(&self) -> Vec<u64> {
        let mut counts = vec![0u64; self.n_rows];
        for (i, &r) in self.rows.iter().enumerate() {
            counts[r] += self.values[i] as u64;
        }
        counts
    }

    /// Get number of genes detected per cell
    pub fn genes_per_cell(&self) -> Vec<u64> {
        let mut genes: Vec<ahash::AHashSet<usize>> =
            (0..self.n_cols).map(|_| ahash::AHashSet::new()).collect();
        for (&r, &c) in self.rows.iter().zip(self.cols.iter()) {
            genes[c].insert(r);
        }
        genes.iter().map(|s| s.len() as u64).collect()
    }

    /// Get number of cells expressing each gene
    pub fn cells_per_gene(&self) -> Vec<u64> {
        let mut cells: Vec<ahash::AHashSet<usize>> =
            (0..self.n_rows).map(|_| ahash::AHashSet::new()).collect();
        for (&r, &c) in self.rows.iter().zip(self.cols.iter()) {
            cells[r].insert(c);
        }
        cells.iter().map(|s| s.len() as u64).collect()
    }

    /// Write to Matrix Market format
    pub fn write_mtx<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Header
        writeln!(writer, "%%MatrixMarket matrix coordinate integer general")?;
        writeln!(writer, "%")?;
        writeln!(
            writer,
            "{} {} {}",
            self.n_rows,
            self.n_cols,
            self.values.len()
        )?;

        // Data (1-indexed)
        for (i, ((&r, &c), &v)) in self
            .rows
            .iter()
            .zip(self.cols.iter())
            .zip(self.values.iter())
            .enumerate()
        {
            if i > 0 {
                writeln!(writer)?;
            }
            write!(writer, "{} {} {}", r + 1, c + 1, v)?;
        }

        Ok(())
    }

    /// Write barcodes to file
    pub fn write_barcodes<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        for barcode in &self.barcodes {
            writeln!(writer, "{}", barcode)?;
        }
        Ok(())
    }

    /// Write genes to file
    pub fn write_genes<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        for gene in &self.genes {
            writeln!(writer, "{}\t{}", gene, gene)?; // gene_id, gene_name
        }
        Ok(())
    }
}

impl Default for CountMatrix {
    fn default() -> Self {
        Self::new()
    }
}

/// Gene counter for building count matrix
pub struct GeneCounter {
    /// Barcode -> index mapping
    barcode_index: AHashMap<String, usize>,
    /// Gene -> index mapping
    gene_index: AHashMap<String, usize>,
    /// Counts: (gene_idx, cell_idx) -> count
    counts: AHashMap<(usize, usize), u32>,
    /// Barcodes in order
    barcodes: Vec<String>,
    /// Genes in order
    genes: Vec<String>,
}

impl GeneCounter {
    pub fn new() -> Self {
        Self {
            barcode_index: AHashMap::new(),
            gene_index: AHashMap::new(),
            counts: AHashMap::new(),
            barcodes: Vec::new(),
            genes: Vec::new(),
        }
    }

    /// Add a count for a barcode-gene pair
    pub fn add_count(&mut self, barcode: &str, gene: &str, count: u32) {
        let cell_idx = *self.barcode_index.entry(barcode.to_string()).or_insert_with(|| {
            let idx = self.barcodes.len();
            self.barcodes.push(barcode.to_string());
            idx
        });

        let gene_idx = *self.gene_index.entry(gene.to_string()).or_insert_with(|| {
            let idx = self.genes.len();
            self.genes.push(gene.to_string());
            idx
        });

        *self.counts.entry((gene_idx, cell_idx)).or_insert(0) += count;
    }

    /// Increment count by 1
    pub fn increment(&mut self, barcode: &str, gene: &str) {
        self.add_count(barcode, gene, 1);
    }

    /// Build the count matrix
    pub fn build(self) -> CountMatrix {
        let n_rows = self.genes.len();
        let n_cols = self.barcodes.len();

        let mut rows = Vec::with_capacity(self.counts.len());
        let mut cols = Vec::with_capacity(self.counts.len());
        let mut values = Vec::with_capacity(self.counts.len());

        for ((gene_idx, cell_idx), count) in self.counts {
            rows.push(gene_idx);
            cols.push(cell_idx);
            values.push(count);
        }

        CountMatrix {
            barcodes: self.barcodes,
            genes: self.genes,
            rows,
            cols,
            values,
            n_rows,
            n_cols,
        }
    }

    /// Get number of cells
    pub fn num_cells(&self) -> usize {
        self.barcodes.len()
    }

    /// Get number of genes
    pub fn num_genes(&self) -> usize {
        self.genes.len()
    }
}

impl Default for GeneCounter {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gene_counter() {
        let mut counter = GeneCounter::new();

        counter.increment("CELL1", "GENE1");
        counter.increment("CELL1", "GENE1");
        counter.increment("CELL1", "GENE2");
        counter.increment("CELL2", "GENE1");

        let matrix = counter.build();

        assert_eq!(matrix.n_rows, 2);
        assert_eq!(matrix.n_cols, 2);
        assert_eq!(matrix.values.len(), 3);
    }

    #[test]
    fn test_count_matrix_stats() {
        let barcodes = vec!["CELL1".to_string(), "CELL2".to_string()];
        let genes = vec!["GENE1".to_string(), "GENE2".to_string()];
        let data = vec![vec![10, 5], vec![3, 8]];

        let matrix = CountMatrix::from_dense(barcodes, genes, data);

        let counts_per_cell = matrix.counts_per_cell();
        assert_eq!(counts_per_cell, vec![13, 13]);

        let counts_per_gene = matrix.counts_per_gene();
        assert_eq!(counts_per_gene, vec![15, 11]);
    }
}
