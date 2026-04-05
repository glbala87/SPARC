//! Community detection via label propagation clustering

use ahash::AHashMap;

/// Simple community detection using label propagation
/// Takes a KNN graph and returns cluster assignments
pub fn label_propagation(
    graph: &[Vec<(usize, f64)>],
    resolution: f64,
    max_iterations: usize,
) -> Vec<usize> {
    let n = graph.len();
    let mut labels: Vec<usize> = (0..n).collect();

    for _ in 0..max_iterations {
        let mut changed = false;

        for i in 0..n {
            if graph[i].is_empty() {
                continue;
            }

            let mut label_weights: AHashMap<usize, f64> = AHashMap::new();
            for &(j, dist) in &graph[i] {
                let weight = 1.0 / (dist + 1e-10) * resolution;
                *label_weights.entry(labels[j]).or_insert(0.0) += weight;
            }

            if let Some((&best_label, _)) = label_weights
                .iter()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            {
                if best_label != labels[i] {
                    labels[i] = best_label;
                    changed = true;
                }
            }
        }

        if !changed {
            break;
        }
    }

    // Renumber labels to be contiguous
    let mut label_map: AHashMap<usize, usize> = AHashMap::new();
    let mut next_label = 0;
    for label in &mut labels {
        let new_label = *label_map.entry(*label).or_insert_with(|| {
            let l = next_label;
            next_label += 1;
            l
        });
        *label = new_label;
    }

    labels
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_label_propagation() {
        // Two clear clusters
        let graph = vec![
            vec![(1, 0.1), (2, 0.2)],  // cluster A
            vec![(0, 0.1), (2, 0.15)],  // cluster A
            vec![(0, 0.2), (1, 0.15)],  // cluster A
            vec![(4, 0.1), (5, 0.2)],   // cluster B
            vec![(3, 0.1), (5, 0.15)],  // cluster B
            vec![(3, 0.2), (4, 0.15)],  // cluster B
        ];

        let labels = label_propagation(&graph, 1.0, 100);
        assert_eq!(labels.len(), 6);

        // First 3 should be same cluster, last 3 should be same cluster
        assert_eq!(labels[0], labels[1]);
        assert_eq!(labels[1], labels[2]);
        assert_eq!(labels[3], labels[4]);
        assert_eq!(labels[4], labels[5]);
        assert_ne!(labels[0], labels[3]);
    }
}
