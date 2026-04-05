//! Integration tests for sparc-core

use sparc_core::{
    barcode::{BarcodeCorrector, BarcodeMatch, Whitelist},
    count::{CountMatrix, GeneCounter},
    protocols::{DropSeq, InDrop, Protocol, SciRNA, SmartSeq2, TenX3Prime, TenX5Prime},
    streaming::{StreamConfig, StreamingProcessor},
    ReadStructure,
};

// ===== Protocol Tests =====

#[test]
fn test_tenx_3prime_v3_extraction() {
    let protocol = TenX3Prime::v3();
    assert_eq!(protocol.name(), "10x Genomics 3' Gene Expression");

    // v3: 16bp barcode + 12bp UMI = 28bp minimum
    let seq = b"ACGTACGTACGTACGTAAAAAAAAAAAA_extra_cdna";
    let qual = vec![30u8; seq.len()];

    let result = protocol.extract_r1(seq, &qual);
    assert!(result.is_ok());

    let components = result.unwrap();
    assert_eq!(components.barcode.len(), 16);
    assert_eq!(components.umi.len(), 12);
}

#[test]
fn test_tenx_3prime_v2_extraction() {
    let protocol = TenX3Prime::v2();
    assert_eq!(protocol.version(), "v2");

    // v2: 16bp barcode + 10bp UMI = 26bp minimum
    let seq = b"ACGTACGTACGTACGTAAAAAAAAAA_extra";
    let qual = vec![30u8; seq.len()];

    let result = protocol.extract_r1(seq, &qual);
    assert!(result.is_ok());

    let components = result.unwrap();
    assert_eq!(components.barcode.len(), 16);
    assert_eq!(components.umi.len(), 10);
}

#[test]
fn test_tenx_5prime_v2_extraction() {
    let protocol = TenX5Prime::v2();
    assert_eq!(protocol.name(), "10x Genomics 5' Gene Expression");

    let seq = b"ACGTACGTACGTACGTAAAAAAAAAA_extra";
    let qual = vec![30u8; seq.len()];

    let result = protocol.extract_r1(seq, &qual);
    assert!(result.is_ok());
}

#[test]
fn test_dropseq_extraction() {
    let protocol = DropSeq::new();
    assert_eq!(protocol.name(), "Drop-seq");

    // Drop-seq: 12bp barcode + 8bp UMI = 20bp minimum
    let seq = b"ACGTACGTACGTAAAAAAAA_extra";
    let qual = vec![30u8; seq.len()];

    let result = protocol.extract_r1(seq, &qual);
    assert!(result.is_ok());

    let components = result.unwrap();
    assert_eq!(components.barcode.len(), 12);
    assert_eq!(components.umi.len(), 8);
}

#[test]
fn test_indrop_extraction() {
    let protocol = InDrop::new();
    assert_eq!(protocol.name(), "inDrop");

    // inDrop: 16bp barcode + 6bp UMI = 22bp minimum
    let seq = b"ACGTACGTACGTACGTAAAAAA_extra";
    let qual = vec![30u8; seq.len()];

    let result = protocol.extract_r1(seq, &qual);
    assert!(result.is_ok());

    let components = result.unwrap();
    assert_eq!(components.barcode.len(), 16);
    assert_eq!(components.umi.len(), 6);
}

#[test]
fn test_scirna_extraction() {
    let protocol = SciRNA::new();
    assert_eq!(protocol.name(), "sci-RNA-seq");

    // sci-RNA: 10bp barcode + 8bp UMI = 18bp minimum
    let seq = b"ACGTACGTACAAAAAAAA_extra";
    let qual = vec![30u8; seq.len()];

    let result = protocol.extract_r1(seq, &qual);
    assert!(result.is_ok());

    let components = result.unwrap();
    assert_eq!(components.barcode.len(), 10);
    assert_eq!(components.umi.len(), 8);
}

#[test]
fn test_smartseq2_extraction() {
    let protocol = SmartSeq2::new("MySample".to_string());
    assert_eq!(protocol.name(), "SMART-seq2");

    let seq = b"ACGTACGTACGT";
    let qual = vec![30u8; seq.len()];

    let result = protocol.extract_r1(seq, &qual);
    assert!(result.is_ok());

    let components = result.unwrap();
    assert_eq!(components.barcode_str(), "MySample");
}

#[test]
fn test_extraction_too_short() {
    let protocol = TenX3Prime::v3();
    let seq = b"SHORT";
    let qual = vec![30u8; seq.len()];

    let result = protocol.extract_r1(seq, &qual);
    assert!(result.is_err());
}

// ===== Barcode Tests =====

#[test]
fn test_whitelist_contains() {
    let barcodes = vec!["AAAA".to_string(), "CCCC".to_string(), "GGGG".to_string()];
    let whitelist = Whitelist::from_vec(barcodes).unwrap();

    assert!(whitelist.contains("AAAA"));
    assert!(whitelist.contains("CCCC"));
    assert!(!whitelist.contains("TTTT"));
    assert_eq!(whitelist.len(), 3);
}

#[test]
fn test_barcode_corrector_exact() {
    let barcodes = vec!["AAAA".to_string(), "CCCC".to_string()];
    let whitelist = Whitelist::from_vec(barcodes).unwrap();
    let corrector = BarcodeCorrector::new(whitelist, 1);

    match corrector.match_barcode("AAAA") {
        BarcodeMatch::Exact(bc) => assert_eq!(bc, "AAAA"),
        _ => panic!("Expected exact match"),
    }
}

#[test]
fn test_barcode_corrector_corrected() {
    let barcodes = vec!["AAAA".to_string(), "CCCC".to_string()];
    let whitelist = Whitelist::from_vec(barcodes).unwrap();
    let corrector = BarcodeCorrector::new(whitelist, 1);

    // One mismatch from AAAA
    match corrector.match_barcode("AAAT") {
        BarcodeMatch::Corrected(original, corrected, dist) => {
            assert_eq!(original, "AAAT");
            assert_eq!(corrected, "AAAA");
            assert_eq!(dist, 1);
        }
        _ => panic!("Expected corrected match"),
    }
}

#[test]
fn test_barcode_corrector_no_match() {
    let barcodes = vec!["AAAA".to_string(), "CCCC".to_string()];
    let whitelist = Whitelist::from_vec(barcodes).unwrap();
    let corrector = BarcodeCorrector::new(whitelist, 1);

    // Two mismatches from everything (max_mismatch=1)
    match corrector.match_barcode("TTTT") {
        BarcodeMatch::NoMatch(_) => {}
        _ => panic!("Expected no match"),
    }
}

// ===== Count Matrix Tests =====

#[test]
fn test_gene_counter_basic() {
    let mut counter = GeneCounter::new();

    counter.increment("CELL1", "GENE1");
    counter.increment("CELL1", "GENE1");
    counter.increment("CELL1", "GENE2");
    counter.increment("CELL2", "GENE1");
    counter.increment("CELL2", "GENE3");

    assert_eq!(counter.num_cells(), 2);
    assert_eq!(counter.num_genes(), 3);

    let matrix = counter.build();
    assert_eq!(matrix.n_rows, 3); // 3 genes
    assert_eq!(matrix.n_cols, 2); // 2 cells
    assert_eq!(matrix.values.len(), 4); // 4 non-zero entries
}

#[test]
fn test_gene_counter_add_count() {
    let mut counter = GeneCounter::new();

    counter.add_count("CELL1", "GENE1", 5);
    counter.add_count("CELL1", "GENE1", 3);

    let matrix = counter.build();
    assert_eq!(matrix.get(0, 0), 8); // 5 + 3
}

#[test]
fn test_count_matrix_from_dense() {
    let barcodes = vec!["C1".to_string(), "C2".to_string(), "C3".to_string()];
    let genes = vec!["G1".to_string(), "G2".to_string()];
    let data = vec![vec![10, 0, 5], vec![0, 8, 3]];

    let matrix = CountMatrix::from_dense(barcodes, genes, data);

    assert_eq!(matrix.n_rows, 2);
    assert_eq!(matrix.n_cols, 3);
    assert_eq!(matrix.get(0, 0), 10);
    assert_eq!(matrix.get(0, 1), 0);
    assert_eq!(matrix.get(0, 2), 5);
    assert_eq!(matrix.get(1, 1), 8);
}

#[test]
fn test_count_matrix_stats() {
    let barcodes = vec!["C1".to_string(), "C2".to_string()];
    let genes = vec!["G1".to_string(), "G2".to_string(), "G3".to_string()];
    let data = vec![vec![10, 5], vec![3, 8], vec![0, 2]];

    let matrix = CountMatrix::from_dense(barcodes, genes, data);

    let counts_per_cell = matrix.counts_per_cell();
    assert_eq!(counts_per_cell, vec![13, 15]);

    let counts_per_gene = matrix.counts_per_gene();
    assert_eq!(counts_per_gene, vec![15, 11, 2]);

    let genes_per_cell = matrix.genes_per_cell();
    assert_eq!(genes_per_cell, vec![2, 3]);

    let cells_per_gene = matrix.cells_per_gene();
    assert_eq!(cells_per_gene, vec![2, 2, 1]);
}

#[test]
fn test_count_matrix_to_csr() {
    let barcodes = vec!["C1".to_string(), "C2".to_string()];
    let genes = vec!["G1".to_string(), "G2".to_string()];
    let data = vec![vec![10, 5], vec![3, 8]];

    let matrix = CountMatrix::from_dense(barcodes, genes, data);
    let csr = matrix.to_csr();

    assert_eq!(csr.n_rows, 2);
    assert_eq!(csr.n_cols, 2);
    assert_eq!(csr.indptr.len(), 3); // n_rows + 1
}

#[test]
fn test_count_matrix_empty() {
    let matrix = CountMatrix::new();
    assert_eq!(matrix.n_rows, 0);
    assert_eq!(matrix.n_cols, 0);
    assert_eq!(matrix.values.len(), 0);
}

// ===== Streaming Tests =====

#[test]
fn test_stream_config_defaults() {
    let config = StreamConfig::default();
    assert_eq!(config.chunk_size, 1_000_000);
    assert_eq!(config.max_memory, 4 * 1024 * 1024 * 1024);
}

#[test]
fn test_streaming_processor_creation() {
    let config = StreamConfig {
        chunk_size: 50_000,
        max_memory: 2 * 1024 * 1024 * 1024,
    };
    let _processor = StreamingProcessor::new(config);
}

// ===== Analysis Tests =====

#[test]
fn test_normalize_total() {
    use sparc_core::analysis::normalize::normalize_total;

    let mut data = vec![vec![10.0, 20.0, 30.0], vec![5.0, 15.0, 25.0]];
    normalize_total(&mut data, 100.0);

    let sum0: f64 = data[0].iter().sum();
    let sum1: f64 = data[1].iter().sum();

    assert!((sum0 - 100.0).abs() < 1e-6);
    assert!((sum1 - 100.0).abs() < 1e-6);
}

#[test]
fn test_log1p_transform() {
    use sparc_core::analysis::normalize::log1p_transform;

    let mut data = vec![vec![0.0, 1.0, 10.0]];
    log1p_transform(&mut data);

    assert!((data[0][0] - 0.0f64.ln_1p()).abs() < 1e-6);
    assert!((data[0][1] - 1.0f64.ln_1p()).abs() < 1e-6);
    assert!((data[0][2] - 10.0f64.ln_1p()).abs() < 1e-6);
}

#[test]
fn test_scale() {
    use sparc_core::analysis::normalize::scale;

    let mut data = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
    scale(&mut data, None);

    // After scaling, each feature (column) should have zero mean
    let col0_mean: f64 = data.iter().map(|r| r[0]).sum::<f64>() / 3.0;
    assert!(col0_mean.abs() < 1e-6);
}

#[test]
fn test_pca() {
    use sparc_core::analysis::pca::pca;

    // Simple 2D data with clear first PC
    let data = vec![
        vec![1.0, 2.0],
        vec![2.0, 4.0],
        vec![3.0, 6.0],
        vec![4.0, 8.0],
    ];

    let (coords, variances) = pca(&data, 1);
    assert_eq!(coords.len(), 4);
    assert_eq!(coords[0].len(), 1);
    assert_eq!(variances.len(), 1);
    assert!(variances[0] > 0.0);
}

#[test]
fn test_knn_graph() {
    use sparc_core::analysis::neighbors::build_knn_graph;

    let coords = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![10.0, 10.0],
    ];

    let graph = build_knn_graph(&coords, 2);
    assert_eq!(graph.len(), 4);

    // Each point should have 2 neighbors
    for neighbors in &graph {
        assert_eq!(neighbors.len(), 2);
    }
}

#[test]
fn test_label_propagation() {
    use sparc_core::analysis::cluster::label_propagation;
    use sparc_core::analysis::neighbors::build_knn_graph;

    // Two clear clusters
    let coords = vec![
        vec![0.0, 0.0],
        vec![0.1, 0.1],
        vec![0.2, 0.0],
        vec![10.0, 10.0],
        vec![10.1, 10.1],
        vec![10.2, 10.0],
    ];

    let graph = build_knn_graph(&coords, 2);
    let labels = label_propagation(&graph, 1.0, 100);

    assert_eq!(labels.len(), 6);
    // Points 0,1,2 should be in same cluster
    assert_eq!(labels[0], labels[1]);
    assert_eq!(labels[1], labels[2]);
    // Points 3,4,5 should be in same cluster
    assert_eq!(labels[3], labels[4]);
    assert_eq!(labels[4], labels[5]);
    // The two clusters should be different
    assert_ne!(labels[0], labels[3]);
}

// ===== Read Structure Tests =====

#[test]
fn test_read_structure() {
    let rs = ReadStructure::new(0, 16, 16, 12, 0);
    assert_eq!(rs.barcode_start, 0);
    assert_eq!(rs.barcode_len, 16);
    assert_eq!(rs.umi_start, 16);
    assert_eq!(rs.umi_len, 12);
}

// ===== MTX Write/Read Roundtrip =====

#[test]
fn test_mtx_write() {
    use std::io::Read;
    let dir = tempfile::tempdir().unwrap();

    let barcodes = vec!["C1".to_string(), "C2".to_string()];
    let genes = vec!["G1".to_string(), "G2".to_string()];
    let data = vec![vec![10, 5], vec![0, 8]];

    let matrix = CountMatrix::from_dense(barcodes, genes, data);

    matrix.write_mtx(dir.path().join("matrix.mtx")).unwrap();
    matrix
        .write_barcodes(dir.path().join("barcodes.tsv"))
        .unwrap();
    matrix.write_genes(dir.path().join("genes.tsv")).unwrap();

    // Verify files exist
    assert!(dir.path().join("matrix.mtx").exists());
    assert!(dir.path().join("barcodes.tsv").exists());
    assert!(dir.path().join("genes.tsv").exists());

    // Verify MTX header
    let mut mtx_content = String::new();
    std::fs::File::open(dir.path().join("matrix.mtx"))
        .unwrap()
        .read_to_string(&mut mtx_content)
        .unwrap();
    assert!(mtx_content.starts_with("%%MatrixMarket"));
}
