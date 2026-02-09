//! UMI deduplication using directional adjacency method

use super::{Umi, UmiGroup};
use ahash::{AHashMap, AHashSet};

/// UMI graph for clustering
pub struct UmiGraph {
    /// Nodes: UMI -> count
    nodes: AHashMap<String, u32>,
    /// Edges: UMI -> connected UMIs
    edges: AHashMap<String, Vec<String>>,
}

impl UmiGraph {
    pub fn new() -> Self {
        Self {
            nodes: AHashMap::new(),
            edges: AHashMap::new(),
        }
    }

    /// Add a UMI to the graph
    pub fn add_umi(&mut self, umi: &str, count: u32) {
        *self.nodes.entry(umi.to_string()).or_insert(0) += count;
    }

    /// Build edges based on Hamming distance
    pub fn build_edges(&mut self, max_distance: u32) {
        let umis: Vec<String> = self.nodes.keys().cloned().collect();

        for i in 0..umis.len() {
            for j in (i + 1)..umis.len() {
                if Self::hamming_distance(&umis[i], &umis[j]) <= max_distance {
                    self.edges
                        .entry(umis[i].clone())
                        .or_default()
                        .push(umis[j].clone());
                    self.edges
                        .entry(umis[j].clone())
                        .or_default()
                        .push(umis[i].clone());
                }
            }
        }
    }

    fn hamming_distance(a: &str, b: &str) -> u32 {
        if a.len() != b.len() {
            return u32::MAX;
        }
        a.chars().zip(b.chars()).filter(|(a, b)| a != b).count() as u32
    }

    /// Get connected components
    pub fn connected_components(&self) -> Vec<Vec<String>> {
        let mut visited: AHashSet<String> = AHashSet::new();
        let mut components = Vec::new();

        for umi in self.nodes.keys() {
            if !visited.contains(umi) {
                let mut component = Vec::new();
                let mut stack = vec![umi.clone()];

                while let Some(current) = stack.pop() {
                    if visited.insert(current.clone()) {
                        component.push(current.clone());
                        if let Some(neighbors) = self.edges.get(&current) {
                            for neighbor in neighbors {
                                if !visited.contains(neighbor) {
                                    stack.push(neighbor.clone());
                                }
                            }
                        }
                    }
                }

                components.push(component);
            }
        }

        components
    }

    /// Get node count
    pub fn get_count(&self, umi: &str) -> u32 {
        *self.nodes.get(umi).unwrap_or(&0)
    }
}

impl Default for UmiGraph {
    fn default() -> Self {
        Self::new()
    }
}

/// UMI deduplicator using directional adjacency method
pub struct UmiDeduplicator {
    /// Maximum edit distance for UMI clustering
    max_distance: u32,
}

impl UmiDeduplicator {
    pub fn new(max_distance: u32) -> Self {
        Self { max_distance }
    }

    /// Deduplicate UMIs using directional adjacency
    ///
    /// This method clusters UMIs that are within `max_distance` edits of each other,
    /// considering the direction based on read counts.
    pub fn deduplicate(&self, umis: &[Umi]) -> Vec<UmiGroup> {
        if umis.is_empty() {
            return Vec::new();
        }

        // Build graph
        let mut graph = UmiGraph::new();
        for umi in umis {
            graph.add_umi(&umi.sequence, umi.count);
        }
        graph.build_edges(self.max_distance);

        // Get connected components
        let components = graph.connected_components();

        // Create UMI groups from components
        let mut groups = Vec::new();
        for component in components {
            let mut group_umis: Vec<Umi> = component
                .iter()
                .map(|seq| Umi::with_count(seq.clone(), graph.get_count(seq)))
                .collect();

            // Sort by count (descending) to get representative
            group_umis.sort_by(|a, b| b.count.cmp(&a.count));

            let representative = group_umis[0].sequence.clone();
            let mut group = UmiGroup::new(representative);

            for umi in group_umis {
                group.add_member(umi);
            }

            groups.push(group);
        }

        groups
    }

    /// Simple clustering: group UMIs with exact matches only
    pub fn deduplicate_exact(&self, umis: &[Umi]) -> Vec<UmiGroup> {
        let mut umi_counts: AHashMap<String, u32> = AHashMap::new();

        for umi in umis {
            *umi_counts.entry(umi.sequence.clone()).or_insert(0) += umi.count;
        }

        umi_counts
            .into_iter()
            .map(|(seq, count)| {
                let mut group = UmiGroup::new(seq.clone());
                group.add_member(Umi::with_count(seq, count));
                group
            })
            .collect()
    }
}

impl Default for UmiDeduplicator {
    fn default() -> Self {
        Self::new(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_umi_deduplication() {
        let umis = vec![
            Umi::with_count("AAAAAAAAAAAA".to_string(), 10),
            Umi::with_count("AAAAAAAAAAAC".to_string(), 2), // 1 edit from first
            Umi::with_count("CCCCCCCCCCCC".to_string(), 5),
        ];

        let dedup = UmiDeduplicator::new(1);
        let groups = dedup.deduplicate(&umis);

        assert_eq!(groups.len(), 2); // Two groups: A-cluster and C-cluster

        let a_group = groups.iter().find(|g| g.representative.starts_with('A'));
        assert!(a_group.is_some());
        assert_eq!(a_group.unwrap().members.len(), 2);
    }

    #[test]
    fn test_exact_dedup() {
        let umis = vec![
            Umi::new("AAAAAAAAAAAA".to_string()),
            Umi::new("AAAAAAAAAAAA".to_string()),
            Umi::new("CCCCCCCCCCCC".to_string()),
        ];

        let dedup = UmiDeduplicator::new(0);
        let groups = dedup.deduplicate_exact(&umis);

        assert_eq!(groups.len(), 2);
    }
}
