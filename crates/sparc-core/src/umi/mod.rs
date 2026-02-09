//! UMI processing module

mod dedup;

pub use dedup::{UmiDeduplicator, UmiGraph};

/// A UMI with associated data
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Umi {
    /// UMI sequence
    pub sequence: String,
    /// Read count for this UMI
    pub count: u32,
}

impl Umi {
    pub fn new(sequence: String) -> Self {
        Self { sequence, count: 1 }
    }

    pub fn with_count(sequence: String, count: u32) -> Self {
        Self { sequence, count }
    }
}

/// UMI group after deduplication
#[derive(Debug, Clone)]
pub struct UmiGroup {
    /// Representative UMI (highest count)
    pub representative: String,
    /// All UMIs in this group
    pub members: Vec<Umi>,
    /// Total deduplicated count
    pub total_count: u32,
}

impl UmiGroup {
    pub fn new(representative: String) -> Self {
        Self {
            representative,
            members: Vec::new(),
            total_count: 0,
        }
    }

    pub fn add_member(&mut self, umi: Umi) {
        self.total_count += umi.count;
        self.members.push(umi);
    }
}
