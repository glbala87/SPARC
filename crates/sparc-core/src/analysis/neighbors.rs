//! K-nearest neighbors graph construction

/// Build k-nearest neighbor graph from coordinates
/// Returns adjacency list: cell_index -> Vec<(neighbor_index, distance)>
pub fn build_knn_graph(coords: &[Vec<f64>], k: usize) -> Vec<Vec<(usize, f64)>> {
    let n = coords.len();
    let mut graph = Vec::with_capacity(n);

    for i in 0..n {
        let mut distances: Vec<(usize, f64)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (j, euclidean_distance(&coords[i], &coords[j])))
            .collect();

        distances.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        distances.truncate(k);

        graph.push(distances);
    }

    graph
}

fn euclidean_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Convert KNN graph to symmetric connectivity matrix
pub fn knn_to_connectivity(graph: &[Vec<(usize, f64)>], n: usize) -> Vec<Vec<bool>> {
    let mut conn = vec![vec![false; n]; n];
    for (i, neighbors) in graph.iter().enumerate() {
        for &(j, _) in neighbors {
            conn[i][j] = true;
            conn[j][i] = true;
        }
    }
    conn
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_knn_graph() {
        let coords = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![10.0, 10.0],
        ];

        let graph = build_knn_graph(&coords, 1);
        assert_eq!(graph.len(), 3);
        // Point 0's nearest neighbor should be point 1
        assert_eq!(graph[0][0].0, 1);
        // Point 1's nearest neighbor should be point 0
        assert_eq!(graph[1][0].0, 0);
    }

    #[test]
    fn test_connectivity() {
        let coords = vec![vec![0.0], vec![1.0], vec![100.0]];
        let graph = build_knn_graph(&coords, 1);
        let conn = knn_to_connectivity(&graph, 3);

        assert!(conn[0][1]);
        assert!(conn[1][0]);
    }
}
